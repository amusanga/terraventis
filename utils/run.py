import ee
import datetime

# Authenticate and initialize the Earth Engine
ee.Authenticate()
ee.Initialize()

# Define the Area of Interest (AOI) - replace with your coordinates
latitude = 37.7749  # Example latitude
longitude = -122.4194  # Example longitude
aoi = ee.Geometry.Point([longitude, latitude]).buffer(10000)  # 10 km buffer around the point

# Split the AOI into north, central, and south parts
def split_aoi(aoi):
    bounds = aoi.bounds().coordinates().get(0).getInfo()
    ymin = bounds[0][1]
    ymax = bounds[2][1]
    ymid1 = ymin + (ymax - ymin) / 3
    ymid2 = ymin + 2 * (ymax - ymin) / 3

    north = ee.Geometry.Polygon([[[bounds[0][0], ymid2], [bounds[1][0], ymid2], [bounds[2][0], ymax], [bounds[3][0], ymax], [bounds[0][0], ymid2]]])
    central = ee.Geometry.Polygon([[[bounds[0][0], ymid1], [bounds[1][0], ymid1], [bounds[2][0], ymid2], [bounds[3][0], ymid2], [bounds[0][0], ymid1]]])
    south = ee.Geometry.Polygon([[[bounds[0][0], ymin], [bounds[1][0], ymin], [bounds[2][0], ymid1], [bounds[3][0], ymid1], [bounds[0][0], ymin]]])
    
    return north, central, south

north, central, south = split_aoi(aoi)

# Define the date range
end_date = ee.Date(datetime.datetime.now())
start_date = end_date.advance(-2, 'year')

# Filter the Sentinel-2 collection
collection = (ee.ImageCollection('COPERNICUS/S2')
              .filterBounds(aoi)
              .filterDate(start_date, end_date))

# Function to compute NDVI
def compute_ndvi(image):
    return image.normalizedDifference(['B8', 'B4']).rename('NDVI')

# Function to compute SAVI
def compute_savi(image):
    L = 0.5
    savi = image.expression(
        '((NIR - RED) / (NIR + RED + L)) * (1 + L)', {
            'NIR': image.select('B8'),
            'RED': image.select('B4'),
            'L': L
        }).rename('SAVI')
    return image.addBands(savi)

# Function to compute LAI
def compute_lai(image):
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
def compute_tcari_osavi(image):
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
def compute_wdrvi(image):
    alpha = 0.1
    wdrvi = image.expression(
        '(alpha * NIR - RED) / (alpha * NIR + RED)', {
            'NIR': image.select('B8'),
            'RED': image.select('B4'),
            'alpha': alpha
        }).rename('WDRVI')
    return image.addBands(wdrvi)

# Function to compute GNDVI
def compute_gndvi(image):
    gndvi = image.normalizedDifference(['B8', 'B3']).rename('GNDVI')
    return image.addBands(gndvi)

# Function to compute NDMI
def compute_ndmi(image):
    ndmi = image.normalizedDifference(['B8', 'B11']).rename('NDMI')
    return image.addBands(ndmi)

# Function to compute NMDI
def compute_nmdi(image):
    nmdi = image.expression(
        '(NIR - (SWIR1 + SWIR2)) / (NIR + (SWIR1 + SWIR2))', {
            'NIR': image.select('B8'),
            'SWIR1': image.select('B11'),
            'SWIR2': image.select('B12')
        }).rename('NMDI')
    return image.addBands(nmdi)

# Function to compute RECl
def compute_recl(image):
    recl = image.expression(
        'RE / RED - 1', {
            'RE': image.select('B5'),
            'RED': image.select('B4')
        }).rename('RECl')
    return image.addBands(recl)

# Function to compute NDRE
def compute_ndre(image):
    ndre = image.normalizedDifference(['B8', 'B5']).rename('NDRE')
    return image.addBands(ndre)

# Function to compute NDWI
def compute_ndwi(image):
    ndwi = image.normalizedDifference(['B3', 'B8']).rename('NDWI')
    return image.addBands(ndwi)

# Function to compute ARVI
def compute_arvi(image):
    arvi = image.expression(
        '(NIR - 2 * RED + BLUE) / (NIR + 2 * RED + BLUE)', {
            'NIR': image.select('B8'),
            'RED': image.select('B4'),
            'BLUE': image.select('B2')
        }).rename('ARVI')
    return image.addBands(arvi)

# Function to compute EVI
def compute_evi(image):
    evi = image.expression(
        '2.5 * ((NIR - RED) / (NIR + 6 * RED - 7.5 * BLUE + 1))', {
            'NIR': image.select('B8'),
            'RED': image.select('B4'),
            'BLUE': image.select('B2')
        }).rename('EVI')
    return image.addBands(evi)

# Function to compute VARI
def compute_vari(image):
    vari = image.expression(
        '(GREEN - RED) / (GREEN + RED - BLUE)', {
            'GREEN': image.select('B3'),
            'RED': image.select('B4'),
            'BLUE': image.select('B2')
        }).rename('VARI')
    return image.addBands(vari)

# Function to compute GCI
def compute_gci(image):
    gci = image.expression(
        '(NIR / GREEN) - 1', {
            'NIR': image.select('B8'),
            'GREEN': image.select('B3')
        }).rename('GCI')
    return image.addBands(gci)

# Add all computed indices to the image collection
def add_indices(image):
    image = compute_ndvi(image)
    image = compute_savi(image)
    image = compute_lai(image)
    image = compute_tcari_osavi(image)
    image = compute_wdrvi(image)
    image = compute_gndvi(image)
    image = compute_ndmi(image)
    image = compute_nmdi(image)
    image = compute_recl(image)
    image = compute_ndre(image)
    image = compute_ndwi(image)
    image = compute_arvi(image)
    image = compute_evi(image)
    image = compute_vari(image)
    image = compute_gci(image)
    return image

# Apply the function to each image in the collection
indexed_collection = collection.map(add_indices)

# Function to compute statistics for a given region
def compute_stats(image, region, region_name):
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
def get_index_stats(index_name):
    index_collection = indexed_collection.select(index_name)
    north_stats = index_collection.map(lambda img: compute_stats(img, north, 'North')).getInfo()
    central_stats = index_collection.map(lambda img: compute_stats(img, central, 'Central')).getInfo()
    south_stats = index_collection.map(lambda img: compute_stats(img, south, 'South')).getInfo()
    return {
        'North': north_stats,
        'Central': central_stats,
        'South': south_stats
    }

indices = ['NDVI', 'SAVI', 'LAI', 'TCARI_OSAVI', 'WDRVI', 'GNDVI', 'NDMI', 'NMDI', 'RECl', 'NDRE', 'NDWI', 'ARVI', 'EVI', 'VARI', 'GCI']
index_stats = {index: get_index_stats(index) for index in indices}

# Print stats for each region and index
for index, stats in index_stats.items():
    print(f"Index: {index}")
    for region, stat in stats.items():
        print(f"Region: {region}")
        print(stat)
        print("\n")

