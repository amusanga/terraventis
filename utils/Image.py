import ee
import os
import json
import datetime
import numpy as np
from json import JSONEncoder
from tqdm import tqdm

# Authenticate and initialize the Earth Engine
ee.Authenticate()
ee.Initialize(project='ee-terraventis')

class ComputeIndices:
    def __init__(self, output_dir='./output'):
        self.indices = ['NDVI', 'SAVI', 'LAI', 'TCARI_OSAVI', 'WDRVI', 'GNDVI', 'NDMI', 'NMDI', 'RECl', 'NDRE', 'NDWI', 'ARVI', 'EVI', 'VARI', 'GCI']
        self.out_dir = output_dir
        
        if not os.path.exists(self.out_dir):
            os.makedirs(self.out_dir)
    
    def get_date_range(self, days):
        # Define the date range
        self.end_date = datetime.datetime.now()
        self.start_date = self.end_date - datetime.timedelta(days=days)

        # Create a list of datetime objects for each hour between 10 AM and 2 PM for each day in the date range
        self.date_range = []
        current_date =  self.start_date
        while current_date <= self.end_date:
            for hour in [7, 12, 17]:  # 10 AM to 2 PM (14 PM)
                self.date_range.append(current_date.replace(hour=hour, minute=0, second=0))
            current_date += datetime.timedelta(days=1)
        return self.date_range

    # Split the AOI into north, central, and south parts
    def split_aoi(self, aoi):
        bounds = aoi.bounds().coordinates().get(0).getInfo()
        ymin = bounds[0][1]
        ymax = bounds[2][1]
        ymid1 = ymin + (ymax - ymin) / 3
        ymid2 = ymin + 2 * (ymax - ymin) / 3

        north = ee.Geometry.Polygon([[[bounds[0][0], ymid2], [bounds[1][0], ymid2], [bounds[2][0], ymax], [bounds[3][0], ymax], [bounds[0][0], ymid2]]])
        central = ee.Geometry.Polygon([[[bounds[0][0], ymid1], [bounds[1][0], ymid1], [bounds[2][0], ymid2], [bounds[3][0], ymid2], [bounds[0][0], ymid1]]])
        south = ee.Geometry.Polygon([[[bounds[0][0], ymin], [bounds[1][0], ymin], [bounds[2][0], ymid1], [bounds[3][0], ymid1], [bounds[0][0], ymin]]])
        
        return north, central, south
    
    def get_sentinel_collection(self, start_date, end_date):
        # Filter the Sentinel-2 dataset
        self.collection = (ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
                      .filterBounds(self.roi)
                      .filterDate(start_date, end_date)
                      .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 20))
                      .select(['B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B8A', 'B9', 'B11', 'B12']))
        
        return self.collection

        # Function to compute monthly median
    def get_monthly_median(self, year, month):
        start = ee.Date.fromYMD(year, month, 1)
        end = start.advance(1, 'month')
        monthly_collection = self.collection.filterDate(start, end)
        return monthly_collection.median().set('system:time_start', start.millis())

    def get_monthly_sentinel_collections(self, start_date, end_date):

        # Create a list of dates for each month in the date range
        months = ee.List.sequence(start_date.get('year'), end_date.get('year')).map(lambda year: ee.List.sequence(1, 12).map(lambda month: ee.Date.fromYMD(year, month, 1))).flatten()

        print("months : ", months)
        # Compute the monthly median for each month
        monthly_medians = ee.ImageCollection.fromImages(months.map(lambda date: self.get_monthly_median(ee.Date(date).get('year'), ee.Date(date).get('month'))))

        return monthly_medians
                
        # Function to filter the collection for a specific hour
    def get_hourly_landsat_image(self, datetime_obj):
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
    def get_hourly_sentinel_image(self, datetime_obj):
        start_hour = datetime_obj
        end_hour = datetime_obj + datetime.timedelta(hours=1)

        # Filter the Sentinel-2 dataset
        collection = (ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
                    .filterBounds(self.roi)
                    .filterDate(start_hour.strftime('%Y-%m-%dT%H:%M:%S'), end_hour.strftime('%Y-%m-%dT%H:%M:%S'))
                    #.filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 20))
                    .select(['B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B8A', 'B9', 'B11', 'B12']))

        # Get the first image in the collection for the hour
        image = collection.first()

        return image, start_hour.strftime('%Y-%m-%d_%H-%M-%S')

    # Function to compute NDVI
    def compute_ndvi(self, image):
        return image.normalizedDifference(['B8', 'B4']).rename('NDVI')

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
    def add_indices(self, image):

        # if image.getInfo() is None:
        #      print('Image is empty')
        #      return None
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

    # Function to compute statistics for a given region
    def compute_stats(self, image, region, region_name):
        stats = image.reduceRegion(
            reducer=ee.Reducer.minMax().combine(
                reducer2=ee.Reducer.mean(),
                sharedInputs=True
            ),
            geometry=region,
            scale=10,
            maxPixels=1e9
        )
        stats = stats.set('region', region_name)
        return stats

    # Compute statistics for each index in each region
    def get_index_stats(self, index_name,  north, central, south):
        image = self.image_indexes.select(index_name)
        north_stats = self.compute_stats(image, north, 'North').getInfo()
        central_stats = self.compute_stats(image, central, 'Central').getInfo()
        south_stats = self.compute_stats(image, south, 'South').getInfo()
        return {
            'North': north_stats,
            'Central': central_stats,
            'South': south_stats
        }
    
    def computeIndexes_old(self, days, dataset='sentinel'):

        # Get the image collection
        date_range = self.get_date_range(days)
        # Split the AOI into 3 regions
        north, central, south = self.split_aoi(self.roi)
        # Loop through each hour

        for dt in tqdm(date_range, total=len(date_range), desc='Computing indices'):
            if dataset == 'landsat':
                image, date_str = self.get_hourly_landsat_image(dt)
            elif dataset =='sentinel':
                image, date_str = self.get_hourly_sentinel_image(dt)
          
            if image.getInfo():
                self.image_indexes = self.add_indices(image)
            else:
                continue

            out_metadata_path = os.path.join(self.out_dir, f'{dataset}_{date_str}.json')
            filename = f'{dataset}_{date_str}.tif'
            
            
            try:
                index_stats = {index:  self.get_index_stats(index, north, central, south) for index in self.indices}
            except Exception as e:
                print(e)
                continue

            output_json = {
                'filename': filename,
                'index_stats': index_stats
            }
            print(output_json)

            with open(out_metadata_path, "w") as json_file:
                json.dump(
                    output_json,
                    json_file,
                    indent=4,
                    separators=(",", ":"),
                    cls=NumpyArrayEncoder,
                )


    def computeIndexes(self, start_date, end_date):

        collection = self.get_sentinel_collection(start_date, end_date)
        # Split the AOI into 3 regions
        north, central, south = self.split_aoi(self.roi)
        # Loop through each hour
        monthly_medians = self.get_monthly_sentinel_collections(start_date, end_date)
        
        # Apply the function to each monthly median image
        self.image_indexes  = monthly_medians.map(self.add_indices)

        try:
            # Apply the function to each monthly median image
            index_stats = {index:  self.get_index_stats(index, north, central, south) for index in self.indices}
        except Exception as e:
            index_stats = []
            print(e)

        
        output_json = {
            'index_stats': index_stats
        }
        print(output_json)

        # with open(out_metadata_path, "w") as json_file:
        #     json.dump(
        #         output_json,
        #         json_file,
        #         indent=4,
        #         separators=(",", ":"),
        #         cls=NumpyArrayEncoder,
        #     )




class Field(ComputeIndices):
    def __init__(self, json_file):
        super().__init__()
        self.roi = self.read_roi(json_file)
    
    def read_roi(self, json_file):
        '''
        Read the json file containing the coordinates of the field
        '''
        with open(json_file) as f:
            data = json.load(f)

        coords = data['features'][0]['geometry']['coordinates']
        aoi = ee.Geometry.Polygon(coords)
        feature = ee.Feature(aoi, {})
        roi = feature.geometry()
        
        return roi
    
# Applies scaling factors.
def apply_scale_factors(image):
  optical_bands = image.select('SR_B.').multiply(0.0000275).add(-0.2)
  return image.addBands(optical_bands, None, True)

class NumpyArrayEncoder(JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return JSONEncoder.default(self, obj)