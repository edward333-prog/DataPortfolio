# Script to calculate habitable score - A measure of how closely an exoplanet's physical properties resemble Earth's. 
# ...It's based on the assumption that Earth-like conditions are most likely to support life as we know it.

# ESI - Earth similarity index. Radius, Density, Escape velocity, Surface temp  

import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt 
import os

# Set up relative paths
script_dir = os.path.dirname(os.path.abspath(__file__))
project_dir = os.path.dirname(script_dir)
data_dir = os.path.join(project_dir, "Data")
viz_dir = os.path.join(project_dir, "Visualisations")

# Ensure visualization directory exists
os.makedirs(viz_dir, exist_ok=True)

# Load the filtered data
filtered_data = pd.read_csv(os.path.join(data_dir, 'filtered_exoplanet_data.csv'))
overall_data = pd.read_csv(os.path.join(data_dir, 'cumulative_2024.09.30_07.27.05.csv'))
r''
print("\n") # Create a gap for nice readability 
print(f"Loaded {len(filtered_data)} exoplanets for Habitability score")

# Filter out False Positives 
filtered_data = filtered_data[filtered_data['koi_disposition'] != 'FALSE POSITIVE']
print(f"After removing FALSE POSITIVES: {len(filtered_data)} exoplanets remain")


# Define Earth reference values
earth_radius = 1.0 # Earth radii
earth_temp = 288 # Earth temperature in Kelvin
earth_insolation = 1.0 # Earth insolation flux
earth_period = 365.25 # Earth orbital period in days
earth_escape_velocity = 11.186 # km/s

# Constraints for escape velocity calculation 
G = 6.67430e-11 # Gravitational constant in m^3 kg^-1 s^-2
earth_radius_m = 6.371e6 # Earth radius in meters
earth_mass_kg = 5.972e24 # Earth mass in kg

# Define tolerances for different factors (to be assessed for importance)
radius_tolerance = 0.25
temp_tolerance = 0.30
insolation_tolerance = 0.25
period_tolerance = 0.10
star_type_tolerance = 0.10
escape_velocity_tolerance = 0.20

# Adjust tolerances to ensure they sum to 1.0
total_tolerance = radius_tolerance + temp_tolerance + insolation_tolerance + period_tolerance + star_type_tolerance + escape_velocity_tolerance
radius_tolerance /= total_tolerance
temp_tolerance /= total_tolerance
insolation_tolerance /= total_tolerance
period_tolerance /= total_tolerance
star_type_tolerance /= total_tolerance
escape_velocity_tolerance /= total_tolerance

print("\nAdjusted tolerance weights (normalized to sum to 1.0)")
print(f"Radius tolerance: {radius_tolerance:.3f}")
print(f"Temperature tolerance: {temp_tolerance:.3f}")
print(f"Insolation tolerance: {insolation_tolerance:.3f}")
print(f"Period tolerance: {period_tolerance:.3f}")
print(f"Star type tolerance: {star_type_tolerance:.3f}")
print(f"Escape velocity tolerance: {escape_velocity_tolerance:.3f}")

# Function to estimate mass from radius using an empirical relationship for rocky exoplanets 

def estimate_mass_from_radius(radius_earth_units):
    """
    Estimate planet mass from radius using an empirical mass-radius relationship. 
    input: radius in Earth radii
    output: mass in Earth masses 

    For rocky planets (R < 1.5 R_Earth), we use M ∝ R^3.7 relationship
    based on observant exoplanet populations. 
    """
    if radius_earth_units < 1.5:
        return radius_earth_units ** 3.7 # This accounts for the observation that slighly larger rocky planets tend to be denser due to compression 
    else:
        return radius_earth_units ** 2.0 # For larger planets we'll use a different relationship but this is less important for habitability 


# Calculate Earth Similarity Index (ESI) components 

# Radius similarity (closer to Earth's radius is better)
filtered_data['radius_similarity'] = (1 - abs(filtered_data['koi_prad'] - earth_radius) / 
                                     (filtered_data['koi_prad'] + earth_radius)) ** radius_tolerance


# Temperature similarity (closer to Earth's temperature is better)
filtered_data['temp_similarity'] = (1 - abs(filtered_data['koi_teq'] - earth_temp) / 
                                     (filtered_data['koi_teq'] + earth_temp)) ** temp_tolerance

# Insolation similarity (closer to Earth's insolation is better)
filtered_data['insolation_similarity'] = (1 - abs(filtered_data['koi_insol'] - earth_insolation) / 
                                     (filtered_data['koi_insol'] + earth_insolation)) ** insolation_tolerance

# Period similarity (logarithmic scale since periods vary widely)
filtered_data['period_similarity'] = (1 - abs(np.log10(filtered_data['koi_period']) - np.log10(earth_period)) / 
                                     (np.log10(filtered_data['koi_period']) + np.log10(earth_period))) ** period_tolerance

# Star type scoring (G-type stars like our own sun are considered optimal for life)
def score_star_type(temp):
    # G-type stars (5200-6000k) get highest scores 
    if 5200 <= temp <= 6000:
        return 1.0
    # k-type stars (3700-5200k) get second highest scores
    elif 3700 <= temp < 5200:
        return 0.8 
    # F-type starts (6000-7500k) get third highest scores 
    elif 6000 <= temp < 7500:
        return 0.6 
    # M-type stars (2400-3700k) get lowest scores 
    elif 2400 <= temp < 3700:
        return 0.4
    # Other star types get lowest scores  
    else: 
        return 0.2

filtered_data['star_type_score'] = filtered_data['koi_steff'].apply(score_star_type) ** star_type_tolerance

# Calculate escape velocity using empirical mass-radius relationship 
# Note: This is an approximation, with significant uncertainty 
filtered_data['estimated_mass'] = filtered_data['koi_prad'].apply(estimate_mass_from_radius) 

# Convert Earth masses to kg and Earth radii to meters 
filtered_data['mass_kg'] = filtered_data['estimated_mass'] * earth_mass_kg
filtered_data['radius_m'] = filtered_data['koi_prad'] * earth_radius_m

# Calculate escape velocity in km/s 
filtered_data['escape_velocity'] = np.sqrt(2 * G * filtered_data['mass_kg'] / filtered_data['radius_m']) / 1000 

# Print estimated masses and escape velocity 
print("\nEstimated planet Properties:")
print(filtered_data[['kepoi_name', 'koi_prad', 'estimated_mass', 'escape_velocity']].to_string(index=False))
print("\nNote: These mass estimates are appoximately based on empirical relationships and have significant uncertainty.")
print("      (potentially ±50% or more).") 

# Calculate escape velocity similarity (compared to Earth)
filtered_data['escape_velocity_similarity'] = (1 - abs  (filtered_data['escape_velocity'] - earth_escape_velocity) / 
                                     (filtered_data['escape_velocity'] + earth_escape_velocity)) ** escape_velocity_tolerance

# Calculate overall habitability score (geometric mean of all components)
filtered_data['habitability_score'] = (filtered_data['radius_similarity'] * 
                                    filtered_data['temp_similarity'] * 
                                    filtered_data['insolation_similarity'] * 
                                    filtered_data['period_similarity'] * 
                                    filtered_data['star_type_score'] * 
                                    filtered_data['escape_velocity_similarity']) ** (1/6)

"""
Why Use Geometric Mean Instead of Arithmetic Mean?
The geometric mean has several important implications for our habitability score:

Zero Effect: If any single factor is zero (completely unsuitable), the entire habitability score becomes zero. This makes sense because a planet with, for example, a completely unsuitable temperature couldn't support life regardless of how perfect its other characteristics are.

Balanced Approach: Planets need to score reasonably well across all factors to get a high overall score. A planet can't compensate for a very poor score in one area by having excellent scores in others.

Sensitivity to Low Values: The geometric mean is more sensitive to low values than the arithmetic mean. This means that significant deviations from Earth-like conditions in any one factor will substantially reduce the overall score.

Multiplicative Effects: In real planetary systems, habitability factors often interact multiplicatively rather than additively. For example, temperature and atmospheric pressure together determine whether water can exist as a liquid.

Scientific Precedent
The geometric mean is commonly used in the Earth Similarity Index (ESI) and other habitability metrics precisely because of these properties. It better reflects the reality that habitability requires meeting multiple conditions simultaneously, not just on average.
"""


# Scale to 0-110 for easier interpretation 
filtered_data['habitability_score'] = filtered_data['habitability_score'] * 100

# Create a binary habitability classification (For ML classification tasks)
#  Planets with score > 70 are considered potentially habitable  
filtered_data['is_habitable'] = (filtered_data['habitability_score'] > 70).astype(int) 

# Print results 
print("\nHabitability Scores:")
print(filtered_data[['kepoi_name', 'koi_disposition','habitability_score', 'is_habitable']].sort_values('habitability_score', ascending=False))


# Count potentially habitable planets with a percentage

habitability_count = filtered_data['is_habitable'].sum()
overall_count = len(overall_data) 
percentage = (habitability_count / overall_count) * 100
habitability_count = filtered_data['is_habitable'].sum()
print(f"Number of potentially habitable planets: {habitability_count}, out of {overall_count} ({percentage:.2f}%)")



# Create a detailed table showing all similarity scores
print("\nDetailed Similarity Scores (0-1 scale before tolerance):")
detailed_scores = filtered_data[['kepoi_name', 'radius_similarity', 'temp_similarity',
                                'insolation_similarity', 'period_similarity',
                                'star_type_score', 'escape_velocity_similarity',
                                'habitability_score']] 
print(detailed_scores.sort_values('habitability_score', ascending=False). to_string(index=False))

# Save the updated data with habitability scores
filtered_data.to_csv(os.path.join(data_dir, 'exoplanets_score.csv'), index=False)

print("Habitability scores saved to 'exoplanets_score.csv'")

## Refer to 'ExoplanetHabitabilityAnalysis.txt' for more details ##

# Create radar chart for visualizing planet characteristics
def create_radar_chart(planet_data):
    # Radar chart metrics
    categories = ['Radius', 'Temperature', 'Insolation', 'Period', 'Star Type', 'Escape Velocity']
    
    # Get the similarity scores
    values = [
        planet_data['radius_similarity'],
        planet_data['temp_similarity'],
        planet_data['insolation_similarity'],
        planet_data['period_similarity'],
        planet_data['star_type_score'],
        planet_data['escape_velocity_similarity']
    ]
    
    # Number of variables
    N = len(categories)
    
    # Create angles for each category
    angles = [n / float(N) * 2 * np.pi for n in range(N)]
    angles += angles[:1]  # Close the loop
    
    # Add the values for each category (close the loop)
    values += values[:1]
    
    # Set up the plot
    fig, ax = plt.subplots(figsize=(10, 10), subplot_kw=dict(polar=True))
    
    # Create a radial grid for the background with higher resolution
    r_grid = np.linspace(0.8, 1.0, 200)  
    theta_grid = np.linspace(0, 2*np.pi, 200)  
    
    # Create meshgrid
    r, theta = np.meshgrid(r_grid, theta_grid)
    
    # Create a simple color value based on the radius (distance from center)
    # 0 at inner edge (r=0.8) and 1 at outer edge (r=1.0)
    z = (r - 0.8) / 0.2  
    
    # Create a colormap - use the same one for all elements
    cmap = plt.cm.RdYlGn  
    
    # Draw the colormap as a polar pcolormesh with smoother rendering
    ax.pcolormesh(theta, r, z, cmap=cmap, alpha=0.8, shading='gouraud')
    
    # Define explicit colors for each grid line from red to green
    grid_values = [0.8, 0.85, 0.9, 0.95, 1.0]
    
    # Red to green colors matching the RdYlGn colormap
    grid_colors = [
        '#d7191c',  
        '#fdae61',  
        '#ffffbf',  
        '#a6d96a',  
        '#1a9641'   
    ]
    
    # Add colored grid lines with thick, very visible lines
    for r, color in zip(grid_values, grid_colors):
        circle = plt.Circle((0, 0), r, transform=ax.transData._b, 
                           fill=False, color=color, linewidth=3.5, alpha=1.0)
        ax.add_artist(circle)
        
        ax.text(
            np.deg2rad(0),  
            r,
            f"{r:.2f}",
            color='black',
            fontsize=10,
            fontweight='bold',
            ha='center',
            va='center',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8)
        )
    
    ax.set_facecolor('none')
    
    ax.plot(angles, values, linewidth=1.5, color='black', zorder=10)
    
    ax.fill(angles, values, alpha=0.2, color='blue', zorder=9)
    
    ax.set_theta_offset(np.pi / 2)
    ax.set_theta_direction(-1)
    
    for angle in angles[:-1]:  
        ax.plot([angle, angle], [0.8, 1.1], color='black', linewidth=1.5, alpha=0.7, zorder=5)
    
    plt.xticks(angles[:-1], categories, fontsize=12, fontweight='bold', color='black')
    
    ax.set_rticks(grid_values)
    ax.set_yticklabels([])
    
    ax.set_ylim(0.8, 1.0)
    
    plt.title(f"Planet {planet_data['kepoi_name']} Similarity Metrics\nHabitability Score: {planet_data['habitability_score']:.2f}", y=1.1)
    
    return fig

# Create radar charts for each potential candidate and save 
for idx, row in filtered_data.iterrows():
    fig = create_radar_chart(row)
    plt.savefig(os.path.join(viz_dir, f"{row['kepoi_name']}_radar.png"))
    plt.close()

print(f"Radar charts saved to {viz_dir}")