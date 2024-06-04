import rasterio
import numpy as np

# Open the downloaded GeoTIFF file
with rasterio.open('/home/cesaire/Projects/terraventis/LANDSAT_IMAGES/Landsat_2024-05-16_11-00-00.tif') as src:
    # Read the image as a NumPy array
    array = src.read(1)  # Reads the first band

print(array)