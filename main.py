import os
import ee
import datetime
# from utils.Image import Field
from utils.processor import Field

os.system('clear')
# Define the date range

# Define the date range
start_date = '2018-01-01'  # Replace with your start date
end_date = '2024-04-02'  # Replace with your end date
# start = ee.Date(start_date)
# end = ee.Date(end_date)

field = Field("/home/cesaire/Projects/terraventis/farmland/6.geojson")
aggregated_data = field.computeIndexes(start_date, end_date)

print("hbhhbh")