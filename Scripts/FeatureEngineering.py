# Checking if KOI's are in habitable zone and calculating against ESI(earth similarity index)

# Using SQL to filter earth-like size, temperature and insolation flux 
import pandas as pd
import os

# Set up relative paths
script_dir = os.path.dirname(os.path.abspath(__file__))
project_dir = os.path.dirname(script_dir)
data_dir = os.path.join(project_dir, "Data")

# Load the cleaned data
cleaned_data = pd.read_csv(os.path.join(data_dir, 'cleaned_data.csv'))

'''
Initially we standardized the data before the first filter and no results were returned,
so we applied the first filter before standardization and we have 5 potentials returned. 
I then realized we need to apply all feature engineering before standardization (which we only need for training models later on)
'''

# Apply the filtering criteria on the cleaned data
filtered_data = cleaned_data[
    (cleaned_data['koi_prad'] >= 0.5) & (cleaned_data['koi_prad'] <= 1.5) &
    (cleaned_data['koi_teq'] >= 273) & (cleaned_data['koi_teq'] <= 373) &
    (cleaned_data['koi_insol'] >= 0.35) & (cleaned_data['koi_insol'] <= 1.5)
]

# Count the number of exoplanets that match the criteria
count_exoplanets = filtered_data.shape[0]

# Print the result
print(f"Number of exoplanets that fall within the criteria: {count_exoplanets}")

# Save the filtered data
filtered_data.to_csv(os.path.join(data_dir, 'FirstFilterCanditates.csv'), index=False)

# Apply the second habitability filter which involves referencing the stellar temperatures, each star type has its own 'goldilocks zone'
# This is also important to check if the star can provide enough solar energy for the habitats 

# Check star type and its relative goldilocks zone 
df1 = pd.read_csv(os.path.join(data_dir, 'FirstFilterCanditates.csv'))

# You can classify the star type based on the effective temperature (koi_steff):
'''
M-type: 2400 K - 3700 K
K-type: 3700 K - 5200 K
G-type: 5200 K - 6000 K
F-type: 6000 K - 7500 K
'''
# Calc stellar luminosity (Stefan-Boltzmann law)
sigma = 5.67e-8  # Stefan-Boltzmann constant in W/m^2/K^4
solar_radius = 6.96e8  # Solar radius in meters
solar_luminosity = 3.828e26  # Solar luminosity in watts

# drop any NA values
filtered_data = df1.dropna(subset=['koi_steff', 'koi_srad'])

# Classify star type based on effective temperature
def classify_star_type(teff):
    if 2400 <= teff < 3700:
        return 'M-type'
    elif 3700 <= teff < 5200:
        return 'K-type'
    elif 5200 <= teff < 6000:
        return 'G-type'
    elif 6000 <= teff < 7500:
        return 'F-type'
    else:
        return 'Other'

# Apply the classify_star_type function to the DataFrame
df1['star_type'] = df1['koi_steff'].apply(classify_star_type)

# Print the kepid and star_type columns
print(df1[['kepid', 'star_type']])

df2 = pd.read_csv(os.path.join(data_dir, 'FirstFilterCanditates.csv'))

# Drop any NA values
filtered_data = df2.dropna(subset=['koi_steff', 'koi_srad', 'koi_period', 'koi_insol'])

# Calc generic day cycle for each planet
filtered_data['generic_day_length'] = filtered_data['koi_period'] / 365.25

# Convert generic day length into Earth hours
filtered_data['generic_day_length_hours'] = filtered_data['generic_day_length'] * 24 * 365.25

# Estimate hours per day of sunlight (assuming 50/50 day/night cycle)
filtered_data['hours_of_sunlight_per_year'] = filtered_data['generic_day_length_hours'] / 2 

# Calc the total hours of sunlight per year 
filtered_data['hours_of_sunlight_per_year'] = filtered_data['hours_of_sunlight_per_year'] * 365.25

print(filtered_data[['kepid', 'generic_day_length_hours']])

# Save the final filtered data with all calculations
filtered_data.to_csv(os.path.join(data_dir, 'filtered_exoplanet_data.csv'), index=False)

# Count and print remaining records 
count_exoplanets = filtered_data.shape[0]
print(f"Number of exoplanets that fall within the criteria: {count_exoplanets}")
print("Habitability analysis completed. Results saved to 'filtered_exoplanet_data.csv'")
