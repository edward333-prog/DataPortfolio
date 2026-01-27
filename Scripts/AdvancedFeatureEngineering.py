"""

Temporal Stability Analysis and Stella Evolution Impact as part of advanced feature engineering

Temporal Stability Analysis
Here we are testing for Tidal Locked planets and this affect of habitability 
-A tidal locked planet is a planet with matching rotational and orbital periods

Stella Evolution Impact 
How the stars evolution will affect the habitability of the planet over time, and at what stage the planet is in its habitable period

"""

import pandas as pd 
import numpy as np 

import os


# Set up relative paths
script_dir = os.path.dirname(os.path.abspath(__file__))
project_dir = os.path.dirname(script_dir)
data_dir = os.path.join(project_dir, "Data")
viz_dir = os.path.join(project_dir, "Visualisations")

# Ensure visualization directory exists
os.makedirs(viz_dir, exist_ok=True)

# Load data 
scored_data = pd.read_csv(os.path.join(data_dir, 'exoplanets_score.csv'))
print(f"Loaded {len(scored_data)} exoplanets for temporal stability analysis")

# Filter to only include potentially habitable planets (habitability score > 70)
habitable_planets = scored_data[scored_data['habitability_score'] > 70].copy()
print(f"Analyzing {len(habitable_planets)} potentially habitable planets")

# Constants
G = 6.67430e-11  # Gravitational constant in m^3 kg^-1 s^-2
SECONDS_PER_DAY = 86400
EARTH_MASS_KG = 5.972e24  # Earth mass in kg
EARTH_RADIUS_M = 6.371e6  # Earth radius in meters
SOLAR_MASS_KG = 1.989e30  # Solar mass in kg
SOLAR_LUMINOSITY = 3.828e26  # Solar luminosity in watts
EARTH_YEAR_DAYS = 365.25  # Earth year in days
AU_TO_METERS = 149597870700  # 1 AU in meters (AU- Astronomical Unit, The AVG distance from Earth to our Sun)
STEFAN_BOLTZMANN = 5.670374419e-8  # Stefan-Boltzmann constant in W/(m^2·K^4)

# Define functions for tidal locking analysis 
def estimate_star_mass(stellar_temp):
    """
    Estimate star mass from effective temperature using an approximate relation
    Returns mass in solar masses
    """
    if stellar_temp < 3700:  # M-dwarf
        return 0.1 + 0.2 * (stellar_temp / 3700)
    elif stellar_temp < 5200:  # K-dwarf
        return 0.5 + 0.3 * ((stellar_temp - 3700) / 1500)
    elif stellar_temp < 6000:  # G-dwarf
        return 0.8 + 0.2 * ((stellar_temp - 5200) / 800)
    elif stellar_temp < 7500:  # F-dwarf
        return 1.0 + 0.6 * ((stellar_temp - 6000) / 1500)
    else:  # A-dwarf or hotter
        return 1.6 + 0.4 * min(1, (stellar_temp - 7500) / 2500)

def estimate_star_luminosity(stellar_temp, stellar_radius=None):
    """
    Estimate star luminosity from effective temperature using Stefan-Boltzmann law
    Returns luminosity in solar units
    """
    if stellar_radius is None:
        # Estimate radius from temperature if not provided
        # This is a very rough approximation
        if stellar_temp < 3700:  # M-dwarf
            stellar_radius = 0.5  # Solar radii
        elif stellar_temp < 5200:  # K-dwarf
            stellar_radius = 0.7  # Solar radii
        elif stellar_temp < 6000:  # G-dwarf
            stellar_radius = 1.0  # Solar radii
        elif stellar_temp < 7500:  # F-dwarf
            stellar_radius = 1.3  # Solar radii
        else:  # A-dwarf or hotter
            stellar_radius = 1.7  # Solar radii
    
    # Convert radius to meters
    stellar_radius_m = stellar_radius * 6.957e8  # Solar radius in meters
    
    # Calculate luminosity using Stefan-Boltzmann law: L = 4πR²σT⁴
    luminosity = 4 * np.pi * (stellar_radius_m**2) * STEFAN_BOLTZMANN * (stellar_temp**4)
    
    # Convert to solar units
    luminosity_solar = luminosity / SOLAR_LUMINOSITY
    
    return luminosity_solar

"""

The Keplar space telescope prioritised dwarf stars, so the majority of our data is based around this,

M-dwarf (M V): 2400 K - 3700 K
K-dwarf (K V): 3700 K - 5200 K
G-dwarf (G V): 5200 K - 6000 K (like our Sun)
F-dwarf (F V): 6000 K - 7500 K

"""

def calculate_semi_major_axis(period_days, star_mass_solar, planet_mass_earth):
    """
    Calculate semi-major axis using Keplar's Third Law 
    Returns semi-major axis in AU
    """
    # Convert masses to kg
    star_mass_kg = star_mass_solar * SOLAR_MASS_KG
    planet_mass_kg = planet_mass_earth * EARTH_MASS_KG

    # Convert period to seconds
    period_seconds = period_days * SECONDS_PER_DAY

    # Calculate semi-major axis in meters
    semi_major_axis_m = ((G * (star_mass_kg + planet_mass_kg) * period_seconds**2) / (4 * np.pi**2))**(1/3)

    # Convert to AU 
    semi_major_axis_au = semi_major_axis_m / AU_TO_METERS

    return semi_major_axis_au


def estimate_tidal_locking_probability(period_days, stellar_temp, planet_radius, planet_mass):
    """
    Estimate the likelihood of tidal locking based on orbital period and star type
    Returns a value between 0-1 where 1 means highly likely to be tidally locked
    """
    # Simplified model based on orbital period and star type
    # More sophisticated models could include more detailed tidal calculations
    
    # Classify star type based on temperature
    if stellar_temp < 3700:  # M-dwarf
        # M-dwarfs have stronger tidal effects
        if period_days < 50:
            return min(1.0, 2.0 * (1 - (period_days / 50)))
        else:
            return max(0.0, 1.0 - (period_days / 100))
    elif stellar_temp < 5200:  # K-dwarf
        if period_days < 30:
            return min(1.0, 1.5 * (1 - (period_days / 30)))
        else:
            return max(0.0, 0.8 - (period_days / 150))
    else:  # G-dwarf or hotter
        if period_days < 10:
            return min(1.0, (1 - (period_days / 10)))
        else:
            return max(0.0, 0.5 - (period_days / 200))

def estimate_seasonal_variation(eccentricity, semi_major_axis):
    """
    Estimate seasonal temperature variation based on orbital eccentricity
    Returns a percentage variation in insolation between perihelion and aphelion
    Perihelion and aphelion are terms that describe the extreme points in an object's elliptical orbit around,
    the Sun or another star
    """
    # Calculate perihelion and aphelion distances 
   
    perihelion = semi_major_axis * (1 - eccentricity)
    aphelion = semi_major_axis * (1 + eccentricity)
    
    # Calculate insolation at perihelion and aphelion (inverse square law)
    insolation_perihelion = 1 / (perihelion**2)
    insolation_aphelion = 1 / (aphelion**2)
    
    # Calculate percentage variation
    variation_percent = 100 * (insolation_perihelion - insolation_aphelion) / insolation_aphelion
    
    return variation_percent        

def calculate_habitable_zone_boundaries(star_mass, star_luminosity, time_gyr=0):
    """
    Calculate the inner and outer boundaries of the habitable zone
    Returns inner and outer boundaries in AU
    
    time_gyr: time in billions of years from now (0 = present)
    """
    # Stellar evolution model: luminosity increases over time
    # For Sun-like stars, luminosity increases ~10% per billion years
    # For other stars, scale by mass
    luminosity_factor = 1.0 + (time_gyr * 0.1 * (star_mass**0.4))
    evolved_luminosity = star_luminosity * luminosity_factor
    
    # Conservative habitable zone boundaries (Kopparapu et al. 2013)
    # Inner edge: runaway greenhouse
    # Outer edge: maximum greenhouse
    inner_edge = 0.99 * np.sqrt(evolved_luminosity)  # AU
    outer_edge = 1.70 * np.sqrt(evolved_luminosity)  # AU
    
    return inner_edge, outer_edge

def estimate_habitable_lifetime(star_mass, star_luminosity, planet_semi_major_axis, planet_temp, star_type):
    """
    Estimate how long a planet will remain in the habitable temperature range
    Returns an estimated habitable lifetime in billions of years and current phase
    """
    # Define temperature range for habitability
    habitable_temp_min = 273  # Kelvin (0°C)
    habitable_temp_max = 323  # Kelvin (50°C)
    
    # Check if planet currently has habitable temperature
    currently_habitable_temp = habitable_temp_min <= planet_temp <= habitable_temp_max
    
    # Estimate main sequence lifetime in billions of years
    main_sequence_lifetime = 10 * (star_mass**-2.5)  # Sun's lifetime ~10 billion years
    
    # Project temperature evolution over time
    future_times = np.linspace(0, main_sequence_lifetime, 100)  # Sample points throughout star's lifetime
    future_temps = []
    
    for t in future_times:
        # Calculate luminosity increase factor based on star type
        if star_type == "M-dwarf":
            lum_factor = 1.0 + (t * 0.05)  # M-dwarfs evolve more slowly
        elif star_type == "K-dwarf":
            lum_factor = 1.0 + (t * 0.08)
        else:  # G-dwarf or hotter
            lum_factor = 1.0 + (t * 0.1 * (star_mass**0.4))
            
        # Estimate temperature change (T ∝ L^0.25 in basic radiative equilibrium)
        temp_factor = lum_factor**0.25
        future_temp = planet_temp * temp_factor
        future_temps.append(future_temp)
    
    # Find when planet enters and exits habitable temperature range
    time_enter = 0
    time_exit = main_sequence_lifetime
    entered_habitable_zone = False
    
    for i, temp in enumerate(future_temps):
        is_habitable = habitable_temp_min <= temp <= habitable_temp_max
        
        # If planet is in habitable temperature range at this time
        if is_habitable:
            # If we haven't set an entry time yet, this is when it enters
            if not entered_habitable_zone:
                time_enter = future_times[i]
                entered_habitable_zone = True
        # If planet is outside habitable range and we've previously found it inside
        elif entered_habitable_zone:
            time_exit = future_times[i]
            break
    
    # If never entered habitable zone
    if not entered_habitable_zone:
        time_enter = 0
        time_exit = 0
    
    # Calculate habitable lifetime
    habitable_lifetime = time_exit - time_enter
    
    # Determine current phase in habitable lifetime
    if not currently_habitable_temp and entered_habitable_zone and time_enter > 0:
        phase = "Pre-habitable (temperature-based)"
    elif currently_habitable_temp:
        # Calculate what percentage through its habitable period the planet is
        if habitable_lifetime > 0:
            current_time = 0  # Current time is 0 (now)
            percent_through = (current_time - time_enter) / habitable_lifetime * 100
            if percent_through < 33:
                phase = "Early habitable phase (temperature-based)"
            elif percent_through < 66:
                phase = "Middle habitable phase (temperature-based)"
            else:
                phase = "Late habitable phase (temperature-based)"
        else:
            phase = "Briefly habitable (temperature-based)"
    else:
        if entered_habitable_zone and time_enter <= 0:
            phase = "Post-habitable (temperature-based)"
        else:
            phase = "Never habitable (temperature-based)"
    
    return habitable_lifetime, phase

def calculate_stellar_evolution_impact(star_mass, star_luminosity, planet_semi_major_axis, planet_temp, star_type):
    """
    Calculate how stellar evolution will impact the planet's habitability over time
    Returns a score from 0-1 where 1 means stable habitability over long timescales
    """
    # Define tempurate range for habitability 
    habitable_temp_min = 273 # Kelvin (0°C)
    habitable_temp_max = 323  # Kelvin (50°C)

    # Check if planet currently has habitable tempurature 
    currently_habitable_temp = habitable_temp_min <= planet_temp <= habitable_temp_max

    # Estimate tempurature evolution based on star type and stella evolution 
    future_times = [1,2,5]
    future_temps = []

    for t in future_times:
        # Calculate luminocity increase factor based on star type 
        if star_type == "M-dwarf":
            lum_factor = 1.0 + (t * 0.05) # M-dwarfs evolve more slowly
        elif star_type == "K-dwarf":
            lum_factor = 1.0 + (t * 0.08)
        else: # G-dwarf or hotter 
            lum_factor = 1.0 + (t * 0.1 * (star_mass**0.4))

        # Estimate tempurature change (T ∝ L^0.25 in basic radiative equilibrium)
        temp_factor = lum_factor**0.25
        future_temp = planet_temp * temp_factor
        future_temps.append(future_temp)

 # Count how many future periods maintain habitable temperatures
    habitable_duration = sum(1 for temp in future_temps 
                           if habitable_temp_min <= temp <= habitable_temp_max)
    
    # Calculate stability score
    if not currently_habitable_temp:
        stability_score = 0.0
    else:
        stability_score = 0.5 + (habitable_duration / (2 * len(future_times)))
    
    return stability_score


# Apply calculations to our dataset in new columns 
habitable_planets['star_mass_solar'] = habitable_planets['koi_steff'].apply(estimate_star_mass)
habitable_planets['star_luminosity_solar'] = habitable_planets['koi_steff'].apply(estimate_star_luminosity)
habitable_planets['semi_major_axis_au'] = habitable_planets.apply(
    lambda row: calculate_semi_major_axis(
        row['koi_period'], 
        row['star_mass_solar'], 
        row['estimated_mass']
    ), 
    axis=1
)


# Calculate habitable zone boundaries
habitable_planets[['hz_inner_edge', 'hz_outer_edge']] = habitable_planets.apply(
    lambda row: pd.Series(calculate_habitable_zone_boundaries(
        row['star_mass_solar'], 
        row['star_luminosity_solar']
    )), 
    axis=1
)

# Calculate if planet is in habitable zone
habitable_planets['in_habitable_zone'] = habitable_planets.apply(
    lambda row: row['hz_inner_edge'] <= row['semi_major_axis_au'] <= row['hz_outer_edge'],
    axis=1
)

# Calculate tidal locking probability
habitable_planets['tidal_locking_probability'] = habitable_planets.apply(
    lambda row: estimate_tidal_locking_probability(
        row['koi_period'], 
        row['koi_steff'], 
        row['koi_prad'],
        row['estimated_mass']
    ), 
    axis=1
)

# Calculate seasonal variation (using eccentricity if available, otherwise assume circular orbit)
if 'koi_eccen' in habitable_planets.columns:
    habitable_planets['seasonal_variation_percent'] = habitable_planets.apply(
        lambda row: estimate_seasonal_variation(row['koi_eccen'], row['semi_major_axis_au']), 
        axis=1
    )
else:
    # Assume small eccentricity of 0.02 (similar to Earth) if not available
    habitable_planets['seasonal_variation_percent'] = habitable_planets.apply(
        lambda row: estimate_seasonal_variation(0.02, row['semi_major_axis_au']), 
        axis=1
    )

# Calculate stellar evolution impact
habitable_planets['stellar_evolution_impact'] = habitable_planets.apply(
    lambda row: calculate_stellar_evolution_impact(
        row['star_mass_solar'],
        row['star_luminosity_solar'],
        row['semi_major_axis_au'],
        row['koi_teq'],  # Planet equilibrium temperature
        # Determine star type based on temperature
        "M-dwarf" if row['koi_steff'] < 3700 else 
        "K-dwarf" if row['koi_steff'] < 5200 else
        "G-dwarf" if row['koi_steff'] < 6000 else
        "F-dwarf" if row['koi_steff'] < 7500 else "A-dwarf"
    ),
    axis=1
)

# Calculate estimated habitable lifetime
habitable_planets[['habitable_lifetime_gyr', 'habitable_phase']] = habitable_planets.apply(
    lambda row: pd.Series(estimate_habitable_lifetime(
        row['star_mass_solar'], 
        row['star_luminosity_solar'],
        row['semi_major_axis_au'],
        row['koi_teq'],  # Planet equilibrium temperature
        "M-dwarf" if row['koi_steff'] < 3700 else 
        "K-dwarf" if row['koi_steff'] < 5200 else
        "G-dwarf" if row['koi_steff'] < 6000 else
        "F-dwarf" if row['koi_steff'] < 7500 else "A-dwarf"
    )), 
    axis=1
)

# Adjust habitability score based on temporal and stellar factors
# Moderate tidal locking penalty - being tidally locked reduces habitability
tidal_locking_penalty = 0.2  # 20% penalty for fully tidally locked planets
seasonal_variation_penalty = 0.01  # 1% penalty per 10% seasonal variation
stellar_evolution_weight = 0.15  # 15% weight for stellar evolution impact

habitable_planets['habitability_with_temporal_factors'] = habitable_planets['habitability_score'] * (
    1 - (habitable_planets['tidal_locking_probability'] * tidal_locking_penalty) - 
    (habitable_planets['seasonal_variation_percent'] / 100 * seasonal_variation_penalty) +
    (habitable_planets['stellar_evolution_impact'] * stellar_evolution_weight)
)

# Print results
print("\nTemporal Stability and Stellar Evolution Analysis Results:")
print("=" * 80)
for idx, planet in habitable_planets.iterrows():
    print(f"\nPlanet: {planet['kepoi_name']} ({planet['koi_disposition']})")
    print(f"  Original Habitability Score: {planet['habitability_score']:.1f}%")
    print(f"  Adjusted Score with Temporal & Stellar Factors: {planet['habitability_with_temporal_factors']:.1f}%")
    print(f"  Tidal Locking Probability: {planet['tidal_locking_probability']:.1f} ({planet['tidal_locking_probability']*100:.0f}%)")
    print(f"  Seasonal Variation: {planet['seasonal_variation_percent']:.1f}%")
    print(f"  Semi-major Axis: {planet['semi_major_axis_au']:.2f} AU")
    print(f"  Habitable Zone: {planet['hz_inner_edge']:.2f} AU - {planet['hz_outer_edge']:.2f} AU")
    print(f"  Currently in 'Earth-like' Habitable Zone: {'Yes' if planet['in_habitable_zone'] else 'No'}")
    
    # Check if planet has high habitability score but isn't in habitable zone
    if planet['habitability_score'] > 70 and not planet['in_habitable_zone']:
        print("\nNOTE: Despite being outside the traditional distance-based habitable zone, this planet has")
        print("Earth-like temperature conditions that may support habitability. The terminator zone")
        print("between day and night sides (if tidally locked) could provide a temperate region.\n")
    
    print(f"  Planet Temperature rating: {planet['stellar_evolution_impact']:.2f} (1 = Suitable)")
    print(f"  Estimated Habitable Lifetime: {planet['habitable_lifetime_gyr']:.1f} billion years")
    print(f"  Habitable Phase: {planet['habitable_phase']}")
    print(f"  Orbital Period: {planet['koi_period']:.1f} days")

# Save the results to a CSV file
habitable_planets.to_csv(os.path.join(data_dir, 'exoplanets_temporal_analysis.csv'), index=False)
print("\nResults saved to exoplanets_temporal_analysis.csv")