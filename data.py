import ee
import geemap
import os

ee.Authenticate()

ee.Initialize(project='ee-terraventis')

image = ee.Image("LANDSAT/LE7_TOA_5YEAR/1999_2003")

Map = geemap.Map()
Map

landsat_vis = {"bands": ["B4", "B3", "B2"], "gamma": 1.4}
Map.addLayer(image, landsat_vis, "LE7_TOA_5YEAR/1999_2003", True, 0.7)

# Draw any shapes on the map using the Drawing tools before executing this code block
feature = Map.draw_last_feature

if feature is None:
    geom = ee.Geometry.Polygon(
        [
            [
                [-115.413031, 35.889467],
                [-115.413031, 36.543157],
                [-114.034328, 36.543157],
                [-114.034328, 35.889467],
                [-115.413031, 35.889467],
            ]
        ]
    )
    feature = ee.Feature(geom, {})

roi = feature.geometry()

filename = os.path.join("/home/cesaire/Geo-Data", "landsat.tif")

image = image.clip(roi).unmask()
geemap.ee_export_image(
    image, filename=filename, scale=90, region=roi, file_per_band=False
)
geemap.ee_export_image(
    image, filename=filename, scale=90, region=roi, file_per_band=True
)