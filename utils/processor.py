import ee
import os
import json
import geopandas as gpd

from collections import defaultdict

# Authenticate and initialize the Earth Engine
ee.Authenticate()
ee.Initialize(project='ee-terraventis')


class ComputeIndices:
    def __init__(self, output_dir='./output'):
        self.indices = ['NDVI', 'SAVI', 'LAI', 'TCARI_OSAVI', 'WDRVI', 'GNDVI', 'NDMI', 'NMDI', 'RECl', 'NDRE', 'NDWI', 'ARVI', 'EVI', 'VARI', 'GCI']
        self.out_dir = output_dir
        
        if not os.path.exists(self.out_dir):
            os.makedirs(self.out_dir)
    
    def get_sentinel_collection(self, start_date, end_date):
        # Filter the Sentinel-2 dataset
        collection = (ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
                      .filterBounds(self.roi)
                      .filterDate(start_date, end_date))
                      #.select(['B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B8A', 'B9', 'B11', 'B12']))
        
        return collection

    def get_weather_collections(self, start_date, end_date):
        collections = {
            'production': self.get_production_collection(start_date, end_date),
            'precipitation': self.get_chirps_collection(start_date, end_date),
            'evapotranspiration': self.get_evapotranspiration_collection(start_date, end_date),
            'temperature': self.get_temperature_collection(start_date, end_date),
            'rainfall': self.get_rainfall_collection(start_date, end_date),
            'wind_speed': self.get_wind_speed_collection(start_date, end_date),
        }
        return collections
    

    def get_chirps_collection(self, start_date, end_date):
        # Extract precipitation data from CHIRPS
        collection = (ee.ImageCollection('UCSB-CHG/CHIRPS/DAILY')
                      .filterBounds(self.roi)
                      .filterDate(start_date, end_date)
                      .select('precipitation'))

        return collection
    
    def get_production_collection(self, start_date, end_date):
        # Extract production data from MODIS
        collection = (ee.ImageCollection('MODIS/061/MOD17A2H')
                      .filterBounds(self.roi)
                      .filterDate(start_date, end_date)
                      .select('Gpp'))

        return collection
    
    def get_evapotranspiration_collection(self, start_date, end_date):
        # Extract production data from MODIS
        collection = (ee.ImageCollection('MODIS/061/MOD16A2GF')
                      .filterBounds(self.roi)
                      .filterDate(start_date, end_date)
                      .select('ET'))

        return collection
    
    def get_temperature_collection(self, start_date, end_date):
        # Extract production data from MODIS
        collection = (ee.ImageCollection('MODIS/061/MOD11A2')
                      .filterBounds(self.roi)
                      .filterDate(start_date, end_date)
                      .select('LST_Day_1km'))

        return collection
    
    def get_rainfall_collection(self, start_date, end_date):
        # Extract production data from MODIS
        collection = (ee.ImageCollection('TRMM/3B42')
                      .filterBounds(self.roi)
                      .filterDate(start_date, end_date)
                      .select('precipitation', 'HQprecipitation', 'IRprecipitation'))

        return collection
    
    def get_wind_speed_collection(self, start_date, end_date):

        # Extract production data from MODIS
        collection = (ee.ImageCollection('ECMWF/ERA5_LAND/DAILY_AGGR')
                      .filterBounds(self.roi)
                      .filterDate(start_date, end_date)
                      .select('u_component_of_wind_10m', 'v_component_of_wind_10m'))

        return collection


        # Function to compute monthly collection
    def get_monthly_collection(self, date, collection):
        start = ee.Date(date)
        end = start.advance(1, 'month')
        return collection.filterDate(start, end)

    # Function to compute monthly median
    def get_monthly_median(self, date, collection):
        monthly_collection = self.get_monthly_collection(date, collection)
        return monthly_collection.median().set('system:time_start', ee.Date(date).millis())

    # Create a list of dates for each month in the date range
    def generate_monthly_dates(self, start_date, end_date):
        start = ee.Date(start_date)
        end = ee.Date(end_date)
        months = ee.List.sequence(0, end.difference(start, 'month').subtract(1)).map(lambda n: start.advance(n, 'month'))
        return months

    # Function to compute NDVI
    def compute_ndvi(self, image):
        ndvi = image.normalizedDifference(['B8', 'B4']).rename('NDVI')
        return image.addBands(ndvi)

    # Function to compute SAVI
    def compute_savi(self, image):
        L = 0.5
        savi = image.expression(
            '((NIR - RED) / (NIR + RED + L)) * (1 + L)', {
                'NIR': image.select('B8'),
                'RED': image.select('B4'),
                'L': L
            }).rename('SAVI')
        return image.addBands(savi)

    # Function to compute LAI
    def compute_lai(self, image):
        lai = image.expression(
            '3.618 * EVI - 0.118', {
                'EVI': image.expression(
                    '2.5 * ((NIR - RED) / (NIR + 6 * RED - 7.5 * BLUE + 1))', {
                        'NIR': image.select('B8'),
                        'RED': image.select('B4'),
                        'BLUE': image.select('B2')
                    })
            }).rename('LAI')
        return image.addBands(lai)

    # Function to compute TCARI/OSAVI
    def compute_tcari_osavi(self, image):
        tcari = image.expression(
            '3 * ((RE - RED) - 0.2 * (RE - GREEN) * (RE / RED))', {
                'RE': image.select('B5'),
                'RED': image.select('B4'),
                'GREEN': image.select('B3')
            })
        osavi = image.expression(
            '(NIR - RED) / (NIR + RED + 0.16)', {
                'NIR': image.select('B8'),
                'RED': image.select('B4')
            })
        tcari_osavi = tcari.divide(osavi).rename('TCARI_OSAVI')
        return image.addBands(tcari_osavi)

    # Function to compute WDRVI
    def compute_wdrvi(self, image):
        alpha = 0.1
        wdrvi = image.expression(
            '(alpha * NIR - RED) / (alpha * NIR + RED)', {
                'NIR': image.select('B8'),
                'RED': image.select('B4'),
                'alpha': alpha
            }).rename('WDRVI')
        return image.addBands(wdrvi)

    # Function to compute GNDVI
    def compute_gndvi(self, image):
        gndvi = image.normalizedDifference(['B8', 'B3']).rename('GNDVI')
        return image.addBands(gndvi)

    # Function to compute NDMI
    def compute_ndmi(self, image):
        ndmi = image.normalizedDifference(['B8', 'B11']).rename('NDMI')
        return image.addBands(ndmi)

    # Function to compute NMDI
    def compute_nmdi(self, image):
        nmdi = image.expression(
            '(NIR - (SWIR1 + SWIR2)) / (NIR + (SWIR1 + SWIR2))', {
                'NIR': image.select('B8'),
                'SWIR1': image.select('B11'),
                'SWIR2': image.select('B12')
            }).rename('NMDI')
        return image.addBands(nmdi)

    # Function to compute RECl
    def compute_recl(self, image):
        recl = image.expression(
            'RE / RED - 1', {
                'RE': image.select('B5'),
                'RED': image.select('B4')
            }).rename('RECl')
        return image.addBands(recl)

    # Function to compute NDRE
    def compute_ndre(self, image):
        ndre = image.normalizedDifference(['B8', 'B5']).rename('NDRE')
        return image.addBands(ndre)

    # Function to compute NDWI
    def compute_ndwi(self, image):
        ndwi = image.normalizedDifference(['B3', 'B8']).rename('NDWI')
        return image.addBands(ndwi)

    # Function to compute ARVI
    def compute_arvi(self, image):
        arvi = image.expression(
            '(NIR - 2 * RED + BLUE) / (NIR + 2 * RED + BLUE)', {
                'NIR': image.select('B8'),
                'RED': image.select('B4'),
                'BLUE': image.select('B2')
            }).rename('ARVI')
        return image.addBands(arvi)

    # Function to compute EVI
    def compute_evi(self, image):
        evi = image.expression(
            '2.5 * ((NIR - RED) / (NIR + 6 * RED - 7.5 * BLUE + 1))', {
                'NIR': image.select('B8'),
                'RED': image.select('B4'),
                'BLUE': image.select('B2')
            }).rename('EVI')
        return image.addBands(evi)

    # Function to compute VARI
    def compute_vari(self, image):
        vari = image.expression(
            '(GREEN - RED) / (GREEN + RED - BLUE)', {
                'GREEN': image.select('B3'),
                'RED': image.select('B4'),
                'BLUE': image.select('B2')
            }).rename('VARI')
        return image.addBands(vari)

    # Function to compute GCI
    def compute_gci(self, image):
        gci = image.expression(
            '(NIR / GREEN) - 1', {
                'NIR': image.select('B8'),
                'GREEN': image.select('B3')
            }).rename('GCI')
        return image.addBands(gci)

    # Add all computed indices to the image collection
    def add_vegetation_indices(self, image):
        image = self.compute_ndvi(image)
        image = self.compute_savi(image)
        image = self.compute_lai(image)
        image = self.compute_tcari_osavi(image)
        image = self.compute_wdrvi(image)
        image = self.compute_gndvi(image)
        image = self.compute_ndmi(image)
        image = self.compute_nmdi(image)
        image = self.compute_recl(image)
        image = self.compute_ndre(image)
        image = self.compute_ndwi(image)
        image = self.compute_arvi(image)
        image = self.compute_evi(image)
        image = self.compute_vari(image)
        image = self.compute_gci(image)
        return image
    
        # Function to extract index values and save to JSON
    def extract_and_save_indices(self,image, aoi, index_name):
        stats = image.reduceRegion(
            reducer=ee.Reducer.minMax().combine(
                reducer2=ee.Reducer.mean(),
                sharedInputs=True
            ),
            geometry=aoi,
            scale=10,
            maxPixels=1e9
        ).getInfo()

        return {
            'date': ee.Date(image.get('system:time_start')).format('YYYY-MM').getInfo(),
            f'{index_name}_min': stats.get(f'{index_name}_min'),
            f'{index_name}_max': stats.get(f'{index_name}_max'),
            f'{index_name}_mean': stats.get(f'{index_name}_mean')
        }
    
            # Function to extract index values and save to JSON
    def extract_and_save_weather_indices(self,image, aoi, index_name):
        stats = image.reduceRegion(
            reducer=ee.Reducer.mean(),
            geometry=aoi,
            scale=10,
            maxPixels=1e9
        ).getInfo()


        return {
            'date': ee.Date(image.get('system:time_start')).format('YYYY-MM').getInfo(),
            f'{index_name}_mean': stats.get(f'{index_name}'),
        }
    
    # Function to extract weather indices
    def get_weather_indices(self, collections):
        weather_data = []

        for collection_name, monthly_indices in collections.items():
            monthly_indices_list = monthly_indices.toList(monthly_indices.size())
            for i in range(monthly_indices.size().getInfo()):
                weather_data.append(self.extract_and_save_weather_indices(ee.Image(monthly_indices_list.get(i)), self.roi, collection_name))

        return weather_data

    
    def get_vegetation_indices(self, collection):
        # Apply the function to each monthly median image
        monthly_indices = collection.map(self.add_vegetation_indices)
        # Extract the index values for each month and save to a list
        index_data = []
        index_names = ['NDVI', 'SAVI', 'EVI']

        for index_name in index_names:
            monthly_indices_list = monthly_indices.toList(monthly_indices.size())
            for i in range(monthly_indices.size().getInfo()):
                index_data.append(self.extract_and_save_indices(ee.Image(monthly_indices_list.get(i)), self.roi, index_name))
        return index_data
    
    def get_vegetation_data(self, months, sentil_collection):
        monthly_medians = ee.ImageCollection.fromImages(months.map(lambda date: self.get_monthly_median(date, sentil_collection)))
        aggregate_vegetation_data = self.aggregate_data(self.get_vegetation_indices(monthly_medians))
        return aggregate_vegetation_data
    
    def get_weather_data(self, months, weather_collection):
        monthly_medians = {}
        for collection_name, collection in weather_collection.items():
            monthly_medians[collection_name] = ee.ImageCollection.fromImages(months.map(lambda date: self.get_monthly_median(date, collection)))
        aggregate_weather_data = self.aggregate_data(self.get_weather_indices(monthly_medians))
        return aggregate_weather_data
    
    def aggregate_data(self, index_data):
        aggregated_data = defaultdict(dict)
        for entry in index_data:
            date = entry.pop('date')
            for key, value in entry.items():
                index_name = key.split('_')[0]
                metric = key.split('_')[1]
                if index_name not in aggregated_data[date]:
                    aggregated_data[date][index_name] = {}
                aggregated_data[date][index_name][metric] = value

        # Convert defaultdict to regular dict
        aggregated_data = {date: dict(indices) for date, indices in aggregated_data.items()}
        return aggregated_data

    def computeIndexes(self, start_date, end_date):
        #
        sentil_collection = self.get_sentinel_collection(start_date, end_date)
        weather_collection = self.get_weather_collections(start_date, end_date)

        # Get list of monthly dates
        months = self.generate_monthly_dates(start_date, end_date)
        # Compute the monthly median for each month
        aggregate_vegetation_data = self.get_vegetation_data(months, sentil_collection)
        # aggregate_vegetation_data = None
        aggregate_weather_data = self.get_weather_data(months, weather_collection)

        # Merge the weather and vegetation data
        aggregated_data = {"weather":aggregate_weather_data, "vegetation":aggregate_vegetation_data}
        # Save to a JSON file
        output_file = 'aggregated_output.json'
        with open(output_file, 'w') as json_file:
            json.dump(aggregated_data, json_file, indent=4)

        print(f'Aggregated data saved to {output_file}')

        print('Index data saved to monthly_indices.json')





class Field(ComputeIndices):
    def __init__(self, json_file):
        super().__init__()
        self.roi = self.read_roi(json_file)

    # Function to read the GeoJSON file and extract the coordinates
    def read_roi(self, file_path):
        with open(file_path) as f:
            geojson = json.load(f)
        polygon_coords =  geojson['features'][0]['geometry']['coordinates']

        roi = ee.Geometry.Polygon(polygon_coords)
        return roi

