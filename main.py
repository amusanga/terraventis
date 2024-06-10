import os
import ee
import datetime
# from utils.Image import Field
from utils.processor import Field

os.system('clear')
# Define the date range

# Define the date range
start_date = '2022-01-01'  # Replace with your start date
end_date = '2023-01-01'  # Replace with your end date
# start = ee.Date(start_date)
# end = ee.Date(end_date)

# months = ee.List.sequence(start.get('year'), end.get('year')).map(lambda year: ee.List.sequence(1, 12).map(lambda month: ee.Date.fromYMD(year, month, 1))).flatten()
# # Filter the months to be within the start and end date
# months = months.filter(ee.Filter.calendarRange(start.get('year'), end.get('year'), 'year'))


# # Compute the monthly median for each month
# monthly_medians = ee.ImageCollection.fromImages(months.map(lambda date: get_monthly_median(get_monthly_collection(ee.Date(date).get('year'), ee.Date(date).get('month')), ee.Date(date).get('year'), ee.Date(date).get('month'))))


field_1 = Field("/home/cesaire/Projects/terraventis/farmland/3.geojson")
field_1.computeIndexes(start_date, end_date)

print("hbhhbh")