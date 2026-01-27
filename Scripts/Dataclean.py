import pandas as pd
import os

# Set up relative paths
script_dir = os.path.dirname(os.path.abspath(__file__))
project_dir = os.path.dirname(script_dir)
data_dir = os.path.join(project_dir, "Data")

# Load the data
data = pd.read_csv(os.path.join(data_dir, 'cumulative_2024.09.30_07.27.05.csv'))
print("\nData loaded successfully")


# Define relevant columns for habitability analysis based on your dataset key
relevant_columns = [
    'kepid', 'kepoi_name', 'koi_disposition', 'koi_prad', 'koi_teq', 'koi_insol',
    'koi_period', 'koi_duration', 'koi_depth', 'koi_steff', 'koi_slogg', 'koi_srad',
    'ra', 'dec', 'koi_fpflag_nt', 'koi_fpflag_ss', 'koi_fpflag_co', 'koi_fpflag_ec'
]

# Filter the data to include only the relevant columns
filtered_data = data[relevant_columns]

# Check for missing values in the relevent columns
missing_values = filtered_data.isnull().sum()
print("\nMissing values before clean:")
print(missing_values)
print()  # This adds an empty line

# Drop rows with missing data
cleaned_data = filtered_data.dropna()

# Check for missing values in the cleaned data
missing_values_cleaned = cleaned_data.isnull().sum()
print("Missing values after cleaning:")
print(missing_values_cleaned)
print()  # This adds an empty line

# Save the cleaned data to a new CSV file
cleaned_data.to_csv(os.path.join(data_dir, 'cleaned_data.csv'), index=False)

print()  # This adds an empty line
print("Data cleaning completed. Cleaned data saved to 'cleaned_data.csv'")
