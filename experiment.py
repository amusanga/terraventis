import rasterio
import numpy as np

# Open the downloaded GeoTIFF file
with rasterio.open('home/cesaire/Downloads/field_image.tif') as src:
    # Read the image as a NumPy array
    array = src.read(1)  # Reads the first band

print(array)