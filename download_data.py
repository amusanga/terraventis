import os
import ee
import geemap
import datetime
from tqdm import tqdm

os.system('clear')
# Authenticate and initialize the Earth Engine
ee.Authenticate()

ee.Initialize(project='ee-terraventis')

# Applies scaling factors.
def apply_scale_factors(image):
  optical_bands = image.select('SR_B.').multiply(0.0000275).add(-0.2)
  return image.addBands(optical_bands, None, True)


field_coords = [[30.125561339140923,-1.5585790911711437],
[30.1255935256491,-1.55911533444342],
[30.12586174655059,-1.559147509035414],
[30.12639818835357,-1.5589759112057722],
[30.126548392058403,-1.5583431440880757],
[30.125818831206352,-1.5582144456679896],
[30.125561339140923,-1.5585790911711437]]

# Define the Area of Interest (AOI)
# Replace with your AOI or use a known location
aoi = ee.Geometry.Polygon([field_coords])
feature = ee.Feature(aoi, {})
roi = feature.geometry()

# Define the date range
end_date = datetime.datetime.now()
start_date = end_date - datetime.timedelta(days=10)

# Create a list of datetime objects for each hour between 10 AM and 2 PM for each day in the date range
date_range = []
current_date = start_date
while current_date <= end_date:
    for hour in range(10, 15):  # 10 AM to 2 PM (14 PM)
        date_range.append(current_date.replace(hour=hour, minute=0, second=0))
    current_date += datetime.timedelta(days=1)

# Function to filter the collection for a specific hour
def get_hourly_landsat_image(datetime_obj):
    start_hour = datetime_obj
    end_hour = datetime_obj + datetime.timedelta(hours=1)

    # Filter the Sentinel-2 dataset
    collection = (ee.ImageCollection('LANDSAT/LC09/C02/T1_L2')
                .filterDate(start_hour.strftime('%Y-%m-%dT%H:%M:%S'), end_hour.strftime('%Y-%m-%dT%H:%M:%S'))
                .select(['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7']))
    
    # Get the first image in the collection for the hour
    collection = collection.map(apply_scale_factors)
    image = collection.first()

    return image, start_hour.strftime('%Y-%m-%d_%H-%M-%S')


# Function to filter the collection for a specific hour
def get_hourly_sentinel_image(datetime_obj):
    start_hour = datetime_obj
    end_hour = datetime_obj + datetime.timedelta(hours=1)

    # Filter the Sentinel-2 dataset
    collection = (ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
                .filterDate(start_hour.strftime('%Y-%m-%dT%H:%M:%S'), end_hour.strftime('%Y-%m-%dT%H:%M:%S'))
                .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 20))
                .select(['B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B8A', 'B9', 'B11', 'B12']))

    # Get the first image in the collection for the hour
    image = collection.first()

    return image, start_hour.strftime('%Y-%m-%d_%H-%M-%S')


# Create a directory to save images locally
import os

data_type = 'sentinel'

if data_type == 'landsat':
    output_dir = 'LANDSAT_IMAGES'
elif data_type =='sentinel':
    output_dir = 'SENTINEL_IMAGES'

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# Loop through each hour and download the image
for dt in tqdm(date_range, total=len(date_range)):

    if data_type == 'landsat':
         image, date_str = get_hourly_landsat_image(dt)
         filename = f'{output_dir}/Landsat_{date_str}.tif'
    elif data_type =='sentinel':
         image, date_str = get_hourly_sentinel_image(dt)
         filename = f'{output_dir}/Sentinel_{date_str}.tif'

    if image.getInfo():
        image = image.clip(roi).unmask()
        # Export the image
        geemap.ee_export_image(
            image,
            filename=filename,
            scale=90,
            region=roi,
            file_per_band=False
        )
    else:
        print(f"No image found for {date_str}")
print(f"Images downloaded to {output_dir}")



