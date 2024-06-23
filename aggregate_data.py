import json
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.dates as mdates
from sklearn.preprocessing import MinMaxScaler


# Load JSON data
file_path = 'aggregated_output.json'
with open(file_path) as f:
    data = json.load(f)

# Initialize an empty list to store the extracted data
weather_data = []
vegetation_data = []

# Loop through the data and extract relevant information
for date, info in data.get('weather', {}).items():
    entry = {'date': date}
    if 'precipitation' in info and 'mean' in info['precipitation']:
        entry['precipitation'] = info['precipitation']['mean'].get('precipitation', None)
    if 'evapotranspiration' in info and 'mean' in info['evapotranspiration']:
        entry['evapotranspiration'] = info['evapotranspiration']['mean'].get('ET', None)
    if 'temperature' in info and 'mean' in info['temperature']:
        entry['temperature'] = info['temperature']['mean'].get('LST_Day_1km', None)
    if 'wind' in info and 'speed' in info['wind']:
        entry['wind_u'] = info['wind']['speed'].get('u_component_of_wind_10m', None)
        entry['wind_v'] = info['wind']['speed'].get('v_component_of_wind_10m', None)
    if 'production' in info and 'mean' in info['production']:
        try:
            entry['production'] = info['production']['mean']['Gpp']
        except:
            entry['production'] = None
    weather_data.append(entry)

# Loop through the data and extract relevant information
for date, info in data.get('vegetation', {}).items():
    entry = {'date': date}
    if 'NDVI' in info:
        entry['NDVI_mean'] = info['NDVI']['mean']
        entry['NDVI_min'] = info['NDVI']['min']
        entry['NDVI_max'] = info['NDVI']['max']
        
    if 'SAVI' in info:
        entry['SAVI_mean'] = info['SAVI']['mean']
        entry['SAVI_min'] = info['SAVI']['min']
        entry['SAVI_max'] = info['SAVI']['max']

    if 'EVI' in info:
        entry['EVI_mean'] = info['EVI']['mean']
        entry['EVI_min'] = info['EVI']['min']
        entry['EVI_max'] = info['EVI']['max']

    if 'NDRE' in info:
        entry['NDRE_mean'] = info['NDRE']['mean']
        entry['NDRE_min'] = info['NDRE']['min']
        entry['NDRE_max'] = info['NDRE']['max']

    if 'NDWI' in info:
        entry['NDWI_mean'] = info['NDWI']['mean']
        entry['NDWI_min'] = info['NDWI']['min']
        entry['NDWI_max'] = info['NDWI']['max']

    if 'ARVI' in info:
        entry['ARVI_mean'] = info['ARVI']['mean']
        entry['ARVI_min'] = info['ARVI']['min']
        entry['ARVI_max'] = info['ARVI']['max']

    if 'VARI' in info:
        entry['VARI_mean'] = info['VARI']['mean']
        entry['VARI_min'] = info['VARI']['min']
        entry['VARI_max'] = info['VARI']['max']

    if 'GCI' in info:
        entry['GCI_mean'] = info['GCI']['mean']
        entry['GCI_min'] = info['GCI']['min']
        entry['GCI_max'] = info['GCI']['max']
    
    vegetation_data.append(entry)

# Convert to DataFrame
df_weather_data = pd.DataFrame(weather_data)
df_vegetation_data = pd.DataFrame(vegetation_data)
df_combined_data = pd.merge(df_vegetation_data, df_weather_data, on="date")

# Convert 'date' column to datetime format
df_combined_data['date'] = pd.to_datetime(df_combined_data['date'])

# Extract the year from the date
df_combined_data['year'] = df_combined_data['date'].dt.year

df_combined_data = df_combined_data[df_combined_data['year'] != 2018]


# Normalize metrics
def normalize_metrics(df, metrics):
    scaler = MinMaxScaler()
    df[metrics] = scaler.fit_transform(df[metrics])
    return df

# Plot time series for each metric by year
metrics = [ 'SAVI_mean', 'EVI_mean', 'precipitation','NDVI_mean'] #, 'NDRE_mean'] #, 'NDWI_mean', 'ARVI_mean', 'VARI_mean', 'GCI_mean'] #, 'precipitation', 'evapotranspiration', 'temperature', 'wind_u', 'wind_v', 'production']



# Function to plot multiple metrics in the same subplot
def plot_metrics(metrics):

    if len(metrics) > 1:
        normalized_data = normalize_metrics(df_combined_data.copy(), metrics)
    else:
        normalized_data = df_combined_data.copy()
    years = normalized_data['year'].unique()
    num_years = len(years)
    
    fig, axes = plt.subplots(num_years, 1, figsize=(9, 6 * num_years), sharex=False)
    fig.suptitle(f'Time Series of {", ".join(metrics)} by Year', fontsize=14)
    
    for i, year in enumerate(years):
        yearly_data = normalized_data[normalized_data['year'] == year]
        for metric in metrics:
            axes[i].plot(yearly_data['date'], yearly_data[metric], label=f'{metric} - {year}')
        axes[i].set_title(f'{year}')
        axes[i].xaxis.set_major_locator(mdates.MonthLocator())
        axes[i].xaxis.set_major_formatter(mdates.DateFormatter('%m'))
        axes[i].legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
        axes[i].grid(True)
        plt.setp(axes[i].xaxis.get_majorticklabels())
    
    fig.text(0.04, 0.5, 'Values', va='center', rotation='vertical', fontsize=12)
    plt.xlabel('Dates in months')
    plt.tight_layout(rect=[0.05, 0.03, 1, 0.95], h_pad=12.0)
    plt.show()

plot_metrics(metrics)

# Correlation analysis
metrics = [ 'NDVI_mean','VARI_mean', 'GCI_mean','EVI_mean','precipitation', 'evapotranspiration', 'temperature', 'wind_u', 'wind_v', 'production']
# Drop the 'date' and 'year' columns for correlation analysis
corr_df = df_combined_data[metrics]

# Calculate the correlation matrix
corr_matrix = corr_df.corr()

# Plot heatmap of the correlation matrix
plt.figure(figsize=(10, 8))
sns.heatmap(corr_matrix, annot=True, cmap='coolwarm', fmt='.2f', linewidths=0.5)
plt.title('Correlation Matrix of Weather Metrics')
plt.show()
