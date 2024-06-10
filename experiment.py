import ee
import json
import datetime
import geopandas as gpd

# Authenticate and initialize the Earth Engine
ee.Authenticate()
ee.Initialize(project='ee-terraventis')


# # Function to read the GeoJSON file and extract the coordinates
# def read_geojson(file_path):
#     with open(file_path) as f:
#         geojson = json.load(f)
#     return geojson['features'][0]['geometry']['coordinates']

# # Read the polygon coordinates from the GeoJSON file
# geojson_file = '/home/cesaire/Projects/terraventis/farmland/3.geojson'  # Replace with the path to your GeoJSON file
# polygon_coords = read_geojson(geojson_file)

# # Convert the coordinates to a format suitable for Earth Engine
# aoi = ee.Geometry.Polygon(polygon_coords)

# # Define the date range
# start_date = '2022-01-01'  # Replace with your start date
# end_date = '2023-01-01'  # Replace with your end date

# # Filter the Sentinel-2 collection
# collection = (ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
#               .filterBounds(aoi)
#               .filterDate(start_date, end_date))
#               #.filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 20)))  # Filtering images with less than 20% cloud cover

# # Function to compute monthly collection
# def get_monthly_collection(date):
#     start = ee.Date(date)
#     end = start.advance(1, 'month')
#     return collection.filterDate(start, end)

# # Function to compute monthly median
# def get_monthly_median(date):
#     monthly_collection = get_monthly_collection(date)
#     return monthly_collection.median().set('system:time_start', ee.Date(date).millis())

# # Create a list of dates for each month in the date range
# def generate_monthly_dates(start_date, end_date):
#     start = ee.Date(start_date)
#     end = ee.Date(end_date)
#     months = ee.List.sequence(0, end.difference(start, 'month').subtract(1)).map(lambda n: start.advance(n, 'month'))
#     return months

# # Get list of monthly dates
# months = generate_monthly_dates(start_date, end_date)

# # Compute the monthly median for each month
# monthly_medians = ee.ImageCollection.fromImages(months.map(lambda date: get_monthly_median(date)))

# # Function to compute NDVI
# def compute_ndvi(image):
#     ndvi = image.normalizedDifference(['B8', 'B4']).rename('NDVI')
#     return image.addBands(ndvi)

# # Function to compute SAVI
# def compute_savi(image):
#     L = 0.5
#     savi = image.expression(
#         '((NIR - RED) / (NIR + RED + L)) * (1 + L)', {
#             'NIR': image.select('B8'),
#             'RED': image.select('B4'),
#             'L': L
#         }).rename('SAVI')
#     return image.addBands(savi)

# # Function to compute EVI
# def compute_evi(image):
#     evi = image.expression(
#         '2.5 * ((NIR - RED) / (NIR + 6 * RED - 7.5 * BLUE + 1))', {
#             'NIR': image.select('B8'),
#             'RED': image.select('B4'),
#             'BLUE': image.select('B2')
#         }).rename('EVI')
#     return image.addBands(evi)

# # Add all computed indices to the image collection
# def add_indices(image):
#     image = compute_ndvi(image)
#     image = compute_savi(image)
#     image = compute_evi(image)
#     return image

# # Apply the function to each monthly median image
# monthly_indices = monthly_medians.map(add_indices)

# # Function to extract index values and save to JSON
# def extract_and_save_indices(image, aoi, index_name):
#     stats = image.reduceRegion(
#         reducer=ee.Reducer.minMax().combine(
#             reducer2=ee.Reducer.mean(),
#             sharedInputs=True
#         ),
#         geometry=aoi,
#         scale=10,
#         maxPixels=1e9
#     ).getInfo()


#     return {
#         'date': ee.Date(image.get('system:time_start')).format('YYYY-MM').getInfo(),
#         f'{index_name}_min': stats.get(f'{index_name}_min'),
#         f'{index_name}_max': stats.get(f'{index_name}_max'),
#         f'{index_name}_mean': stats.get(f'{index_name}_mean')
#     }

# # Extract the index values for each month and save to a list
# index_data = []
# index_names = ['NDVI', 'SAVI', 'EVI']

# for index_name in index_names:
#     monthly_indices_list = monthly_indices.toList(monthly_indices.size())
#     for i in range(monthly_indices.size().getInfo()):
#         index_data.append(extract_and_save_indices(ee.Image(monthly_indices_list.get(i)), aoi, index_name))

# # Save the list to a JSON file
# with open('monthly_indices.json', 'w') as json_file:
#     json.dump(index_data, json_file, indent=4)

# print('Index data saved to monthly_indices.json')


# Define the Area of Interest (AOI) - replace with your coordinates
latitude = 37.7749  # Example latitude
longitude = -122.4194  # Example longitude
aoi = ee.Geometry.Point([longitude, latitude]).buffer(10000)  # 10 km buffer around the point

# Define the date range
start_date = '2022-01-01'  # Replace with your start date
end_date = '2022-12-31'  # Replace with your end date

# Function to extract and print data for a given dataset and variable
def extract_climate_data(collection, variable, reducer, scale=1000):
    data = collection.filterDate(start_date, end_date).select(variable).mean().reduceRegion(
        reducer=reducer,
        geometry=aoi,
        scale=scale,
        maxPixels=1e9
    ).getInfo()
    print(f'{variable} ({collection.get("system:id").getInfo()}):', data)

# Extract precipitation data from CHIRPS
chirps_collection = ee.ImageCollection('UCSB-CHG/CHIRPS/DAILY')
extract_climate_data(chirps_collection, 'precipitation', ee.Reducer.mean())

# Extract evapotranspiration data from MODIS
modis_et_collection = ee.ImageCollection('MODIS/061/MOD16A2GF')
extract_climate_data(modis_et_collection, 'ET', ee.Reducer.mean())

# Extract temperature data from MODIS
modis_temp_collection = ee.ImageCollection('MODIS/061/MOD11A2')
extract_climate_data(modis_temp_collection, 'LST_Day_1km', ee.Reducer.mean())

# Extract temperature data from ERA5
era5_collection = ee.ImageCollection('ECMWF/ERA5_LAND/DAILY_AGGR')
extract_climate_data(era5_collection, 'mean_2m_air_temperature', ee.Reducer.mean())