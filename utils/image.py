import ee
import os 
import time
import geemap

os.system('clear')

ee.Authenticate()

ee.Initialize(project='ee-terraventis')
# Test to fetch and print details of a sample image

# Define the coordinates of the vertices of a polygon (agricultural field)

field_coords = [[30.125561339140923,-1.5585790911711437],
[30.1255935256491,-1.55911533444342],
[30.12586174655059,-1.559147509035414],
[30.12639818835357,-1.5589759112057722],
[30.126548392058403,-1.5583431440880757],
[30.125818831206352,-1.5582144456679896],
[30.125561339140923,-1.5585790911711437]]

# Create an ee.Geometry.Polygon object
field_polygon = ee.Geometry.Polygon([field_coords])


def mask_s2_clouds(image):
  """Masks clouds in a Sentinel-2 image using the QA band.

  Args:
      image (ee.Image): A Sentinel-2 image.

  Returns:
      ee.Image: A cloud-masked Sentinel-2 image.
  """
  qa = image.select('QA60')

  # Bits 10 and 11 are clouds and cirrus, respectively.
  cloud_bit_mask = 1 << 10
  cirrus_bit_mask = 1 << 11

  # Both flags should be set to zero, indicating clear conditions.
  mask = (
      qa.bitwiseAnd(cloud_bit_mask)
      .eq(0)
      .And(qa.bitwiseAnd(cirrus_bit_mask).eq(0))
  )

  return image.updateMask(mask).divide(10000)

dataset = (
    ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
    .filterDate('2020-01-01', '2020-01-30')
    # Pre-filter to get less cloudy granules.
    .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 20))
    .map(mask_s2_clouds)
)

# out_dir = os.path.join(os.path.expanduser("~"), "Downloads")

# print(dataset.aggregate_array("system:index").getInfo())

# geemap.ee_export_image_collection(dataset, out_dir=out_dir)

# geemap.ee_export_image_collection_to_drive(dataset, folder="export", scale=10)

# Get the first least cloudy image and clip it to the field polygon
cropped_image = dataset.first().clip(field_polygon)


# def get_band_types(image):
#     band_types = {}
#     bands = image.bandNames().getInfo()
#     for band in bands:
#         band_types[band] = image.select(band).dataType().getInfo()
#     return band_types

# # Print band data types
# print(get_band_types(cropped_image))


def convert_band_types(image):
    bands = image.bandNames().getInfo()
    converted_image = ee.Image()
    for band in bands:
        converted_image = converted_image.addBands(image.select(band).toUint16())
    return converted_image

# Apply the conversion
converted_image = convert_band_types(cropped_image)

# Continue with your export process


print(cropped_image.getInfo())

# Export the image to Google Drive
export_task = ee.batch.Export.image.toDrive(**{
    'image': dataset.toBands(),
    'description': 'agricultural_field_image',
    'folder': 'GEE_Images',
    'fileNamePrefix': 'field_image_2',
    'scale': 10,  # Adjust based on the desired resolution
    'region': field_polygon,
    'fileFormat': 'GeoTIFF'
})

# Start the export task
export_task.start()

# Print the status of the export
while export_task.active():
    print('Exporting image... Please wait.')
    time.sleep(30)  # Sleep for 30 seconds; adjust as necessary

print('Export task completed:', export_task.status())
