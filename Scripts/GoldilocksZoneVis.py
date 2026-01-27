import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import plotly.graph_objects as go
import plotly.io as pio

# Set up relative paths
script_dir = os.path.dirname(os.path.abspath(__file__))
project_dir = os.path.dirname(script_dir)
data_dir = os.path.join(project_dir, "Data")
viz_dir = os.path.join(project_dir, "Visualisations")

# Ensure visualization directory exists
os.makedirs(viz_dir, exist_ok=True)

# Define star type temperature ranges
M_dwarf = "2400K-3700K"
K_dwarf = "3700K-5200K" 
G_dwarf = "5200K-6000K"  # Sun-like
F_dwarf = "6000K-7500K"
A_dwarf = "7500K-10000K"

# Extract min and max temperatures for each star type
def extract_temp_range(temp_range_str):
    temp_min, temp_max = temp_range_str.replace('K', '').split('-')
    return int(temp_min), int(temp_max)

# Star type representative temperatures (middle of ranges)
star_types = {
    "M-dwarf": extract_temp_range(M_dwarf),
    "K-dwarf": extract_temp_range(K_dwarf),
    "G-dwarf": extract_temp_range(G_dwarf),  # Sun-like
    "F-dwarf": extract_temp_range(F_dwarf),
    "A-dwarf": extract_temp_range(A_dwarf)
}

# Function to estimate star mass based on temperature (simplified version)
def estimate_star_mass(star_temp):
    """Estimate stellar mass in solar masses based on effective temperature"""
    if star_temp < 3700:  # M-dwarf
        return 0.1 + (star_temp - 2400) / (3700 - 2400) * 0.4
    elif star_temp < 5200:  # K-dwarf
        return 0.5 + (star_temp - 3700) / (5200 - 3700) * 0.3
    elif star_temp < 6000:  # G-dwarf
        return 0.8 + (star_temp - 5200) / (6000 - 5200) * 0.2
    elif star_temp < 7500:  # F-dwarf
        return 1.0 + (star_temp - 6000) / (7500 - 6000) * 0.5
    else:  # A-dwarf
        return 1.5 + (star_temp - 7500) / (10000 - 7500) * 1.0

# Function to estimate star luminosity based on temperature
def estimate_star_luminosity(star_temp):
    """Estimate stellar luminosity in solar luminosities based on effective temperature"""
    # Stefan-Boltzmann Law: L ∝ R² × T⁴
    # Approximating R based on mass-radius relationships for main sequence stars
    star_mass = estimate_star_mass(star_temp)
    
    # Simplified mass-radius-luminosity relationship
    if star_temp < 3700:  # M-dwarf
        return (star_temp / 5778)**4 * (star_mass**0.8)
    elif star_temp < 5200:  # K-dwarf
        return (star_temp / 5778)**4 * (star_mass**0.9)
    else:  # G, F, A dwarfs
        return (star_temp / 5778)**4 * (star_mass**1.0)

# Function to calculate habitable zone boundaries
def calculate_habitable_zone_boundaries(star_temp):
    """
    Calculate the inner and outer boundaries of the habitable zone
    Returns boundaries in AU
    """
    # Estimate star properties
    star_luminosity = estimate_star_luminosity(star_temp)
    
    # Calculate habitable zone boundaries
    # Inner edge (runaway greenhouse) - conservative estimate
    inner_edge = 0.99 * np.sqrt(star_luminosity)
    
    # Outer edge (maximum greenhouse) - conservative estimate
    outer_edge = 1.70 * np.sqrt(star_luminosity)
    
    return inner_edge, outer_edge

# Function to calculate temperature at a given distance from star
def calculate_temperature_at_distance(star_temp, distance_au):
    """
    Calculate equilibrium temperature of planet at given distance from star
    Based on radiative equilibrium: T ∝ L^0.25 / sqrt(d)
    """
    star_luminosity = estimate_star_luminosity(star_temp)
    
    # Earth's equilibrium temperature (~255K without greenhouse effect)
    earth_temp = 255
    
    # Calculate temperature using radiative equilibrium principle
    # T ∝ L^0.25 / sqrt(d)
    temperature = earth_temp * (star_luminosity**0.25) / np.sqrt(distance_au)
    
    return temperature

# Create arrays for visualization
def generate_visualization_data():
    """Generate data for habitable zone visualization"""
    results = []
    
    # Generate data points for each star type
    for star_type, (min_temp, max_temp) in star_types.items():
        # Use average temperature for representative calculations
        star_temp = (min_temp + max_temp) / 2
        
        # Calculate habitable zone
        inner_edge, outer_edge = calculate_habitable_zone_boundaries(star_temp)
        
        # Star properties
        star_mass = estimate_star_mass(star_temp)
        star_luminosity = estimate_star_luminosity(star_temp)
        
        # Generate distances from 0.01 to 5 AU
        distances = np.linspace(0.01, 5, 500)
        
        # Calculate temperatures at each distance
        temperatures = [calculate_temperature_at_distance(star_temp, d) for d in distances]
        
        # Determine habitable temperature range (273K-323K for complex life)
        habitable_temp_min = 273  # 0°C
        habitable_temp_max = 323  # 50°C
        
        # Create temperature-based habitable zone
        temp_habitable_mask = [(temp >= habitable_temp_min and temp <= habitable_temp_max) for temp in temperatures]
        
        # Find temperature-based boundaries
        temp_based_inner_edge = None
        temp_based_outer_edge = None
        
        for i, is_habitable in enumerate(temp_habitable_mask):
            if is_habitable and temp_based_inner_edge is None:
                temp_based_inner_edge = distances[i]
            elif not is_habitable and temp_based_inner_edge is not None and temp_based_outer_edge is None:
                temp_based_outer_edge = distances[i-1]
                break
        
        # If outer edge wasn't found, it extends beyond our distance range
        if temp_based_outer_edge is None and temp_based_inner_edge is not None:
            temp_based_outer_edge = distances[-1]
        
        results.append({
            'star_type': star_type,
            'star_temp': star_temp,
            'star_mass': star_mass,
            'star_luminosity': star_luminosity,
            'distance_based_inner_edge': inner_edge,
            'distance_based_outer_edge': outer_edge,
            'temp_based_inner_edge': temp_based_inner_edge,
            'temp_based_outer_edge': temp_based_outer_edge,
            'distances': distances,
            'temperatures': temperatures,
            'temp_habitable_mask': temp_habitable_mask
        })
    
    return results

# Create visualization
def create_habitability_visualization(show_both_methods=True):
    """
    Create visualization of habitable zones for different star types
    
    Parameters:
    - show_both_methods: If True, show both distance-based and temperature-based habitability
    """
    # Generate data
    viz_data = generate_visualization_data()
    
    # Set up plot
    fig, axs = plt.subplots(len(viz_data), 1, figsize=(12, 15), sharex=True)
    fig.suptitle('Habitable Zones by Star Type', fontsize=16)
    
    # Color scheme
    distance_based_color = 'lightblue'
    temp_based_color = 'lightgreen'
    
    # Plot each star type
    for i, data in enumerate(viz_data):
        ax = axs[i]
        
        # Plot temperature curve
        ax.plot(data['distances'], data['temperatures'], 'r-', linewidth=2)
        
        # Add horizontal lines for habitable temperature range
        ax.axhline(y=273, color='b', linestyle='--', alpha=0.5, label='Freezing (273K)')
        ax.axhline(y=323, color='r', linestyle='--', alpha=0.5, label='Too Hot (323K)')
        
        # Plot distance-based habitable zone
        ax.axvspan(data['distance_based_inner_edge'], data['distance_based_outer_edge'], 
                  alpha=0.3, color=distance_based_color, label='Distance-based HZ')
        
        # Plot temperature-based habitable zone if available
        if data['temp_based_inner_edge'] is not None and show_both_methods:
            ax.axvspan(data['temp_based_inner_edge'], data['temp_based_outer_edge'], 
                      alpha=0.3, color=temp_based_color, label='Temperature-based HZ')
        
        # Set title and labels
        ax.set_title(f"{data['star_type']} (T={int(data['star_temp'])}K, M={data['star_mass']:.2f}M☉, L={data['star_luminosity']:.2f}L☉)")
        ax.set_ylabel('Temperature (K)')
        ax.set_ylim(0, 600)
        
        # Add custom labels for habitable zone boundaries
        ax.annotate(f"Inner: {data['distance_based_inner_edge']:.2f} AU", 
                   xy=(data['distance_based_inner_edge'], 150), 
                   xytext=(data['distance_based_inner_edge'], 80))
        
        ax.annotate(f"Outer: {data['distance_based_outer_edge']:.2f} AU", 
                   xy=(data['distance_based_outer_edge'], 150), 
                   xytext=(data['distance_based_outer_edge'], 120))
        
        # Add legend (only for the first subplot to avoid repetition)
        if i == 0:
            ax.legend(loc='upper right')
    
    # Set common x-label
    axs[-1].set_xlabel('Distance from Star (AU)')
    axs[-1].set_xlim(0, 5)
    
    # Adjust layout
    plt.tight_layout(rect=[0, 0, 1, 0.97])
    
    # Save the plot
    plt.savefig(os.path.join(viz_dir, 'habitable_zones_by_star_type.png'), dpi=300)
    
    # Show the plot
    plt.show()

# Load exoplanet data if available
def load_exoplanet_data():
    try:
        # Try to load the processed data from previous analysis
        return pd.read_csv(os.path.join(data_dir, 'exoplanets_temporal_analysis.csv'))
    except FileNotFoundError:
        print("Exoplanet data file not found. Visualization will not include actual planets.")
        return None

# Plot exoplanets on the habitable zone chart
def create_advanced_visualization_with_planets():
    """Create an advanced visualization that includes actual exoplanets"""
    # Load exoplanet data
    planets_df = load_exoplanet_data()
    
    if planets_df is None:
        # If no planet data, fall back to basic visualization
        create_habitability_visualization()
        return
    
    # Generate star data
    viz_data = generate_visualization_data()
    
    # Create figure
    fig, axs = plt.subplots(len(viz_data), 1, figsize=(14, 18), sharex=True)
    fig.suptitle('Habitable Zones by Star Type with Exoplanet Candidates', fontsize=16)
    
    # Color scheme
    distance_based_color = 'lightblue'
    temp_based_color = 'lightgreen'
    
    # Plot each star type
    for i, data in enumerate(viz_data):
        ax = axs[i]
        
        # Plot temperature curve
        ax.plot(data['distances'], data['temperatures'], 'r-', linewidth=2, label='Temperature')
        
        # Add horizontal lines for habitable temperature range
        ax.axhline(y=273, color='b', linestyle='--', alpha=0.5, label='Freezing (273K)')
        ax.axhline(y=323, color='r', linestyle='--', alpha=0.5, label='Too Hot (323K)')
        
        # Plot distance-based habitable zone
        ax.axvspan(data['distance_based_inner_edge'], data['distance_based_outer_edge'], 
                  alpha=0.3, color=distance_based_color, label='Distance-based HZ')
        
        # Plot temperature-based habitable zone if available
        if data['temp_based_inner_edge'] is not None:
            ax.axvspan(data['temp_based_inner_edge'], data['temp_based_outer_edge'], 
                      alpha=0.3, color=temp_based_color, label='Temperature-based HZ')
        
        # Filter planets for this star type
        star_min_temp, star_max_temp = star_types[data['star_type']]
        star_planets = planets_df[(planets_df['koi_steff'] >= star_min_temp) & 
                                (planets_df['koi_steff'] <= star_max_temp)]
        
        # Plot planets if any exist for this star type
        if not star_planets.empty:
            # Get semi-major axis and temperature
            semi_major_axis = star_planets['semi_major_axis_au']
            planet_temps = star_planets['koi_teq']
            
            # Get habitability scores for color mapping
            habitability_scores = star_planets['habitability_score']
            
            # Create scatter plot
            scatter = ax.scatter(semi_major_axis, planet_temps, 
                               c=habitability_scores, cmap='viridis', 
                               s=50, alpha=0.7, edgecolors='black',
                               vmin=0, vmax=100)
            
            # Add colorbar (only once)
            if i == 0:
                cbar = fig.colorbar(scatter, ax=ax, pad=0.01)
                cbar.set_label('Habitability Score (%)')
            
            # Annotate high-scoring planets
            for idx, row in star_planets[star_planets['habitability_score'] > 90].iterrows():
                ax.annotate(row['kepoi_name'], 
                           xy=(row['semi_major_axis_au'], row['koi_teq']),
                           xytext=(10, 10), textcoords='offset points',
                           fontsize=8, 
                           bbox=dict(facecolor='white', alpha=0.7, boxstyle='round,pad=0.2'),
                           ha='center', va='center', zorder=200)
        
        # Set title and labels
        ax.set_title(f"{data['star_type']} (T={int(data['star_temp'])}K, M={data['star_mass']:.2f}M☉, L={data['star_luminosity']:.2f}L☉)")
        ax.set_ylabel('Temperature (K)')
        ax.set_ylim(0, 600)
        
        # Add zone boundary labels
        ax.annotate(f"Inner: {data['distance_based_inner_edge']:.2f} AU", 
                   xy=(data['distance_based_inner_edge'], 150), 
                   xytext=(data['distance_based_inner_edge'], 80))
        
        ax.annotate(f"Outer: {data['distance_based_outer_edge']:.2f} AU", 
                   xy=(data['distance_based_outer_edge'], 150), 
                   xytext=(data['distance_based_outer_edge'], 120))
        
        # Add legend (only for the first subplot to avoid repetition)
        if i == 0:
            handles, labels = ax.get_legend_handles_labels()
            ax.legend(handles, labels, loc='upper right')
    
    # Set common x-label
    axs[-1].set_xlabel('Distance from Star (AU)')
    axs[-1].set_xlim(0, 5)
    
    # Add text box explaining visualization
    fig.text(0.02, 0.02, 
            "This visualization shows habitable zones for different star types using both distance-based and temperature-based methods.\n"
            "The temperature-based approach (green shading) shows where planets would have Earth-like temperatures (273K-323K).\n"
            "The traditional distance-based habitable zone (blue shading) is calculated using stellar luminosity.\n"
            "Actual exoplanets are plotted with color indicating habitability score.", 
            fontsize=10, bbox=dict(facecolor='white', alpha=0.7))
    
    # Adjust layout
    plt.tight_layout(rect=[0, 0.04, 1, 0.97])
    
    # Save the plot
    plt.savefig(os.path.join(viz_dir, 'habitable_zones_with_exoplanets.png'), dpi=300)
    
    # Show the plot
    plt.show()

# Create a concentric circles visualization with all habitable zones around a single center
def create_concentric_circles_visualization():
    """Create a radial diagram showing habitable zones for different star types as concentric circles"""
    # Generate data
    viz_data = generate_visualization_data()
    
    # Set up plot
    fig, ax = plt.subplots(figsize=(12, 12), facecolor='black')
    fig.suptitle('Radial Habitable Zone Diagram for Different Star Types', fontsize=16, color='white')
    ax.set_facecolor('black')  # Set subplot background to black
    
    # Sort data by star temperature (coolest to hottest)
    viz_data.sort(key=lambda x: x['star_temp'])
    
    # Color scheme for different star types
    star_colors = {
        "M-dwarf": "#FF5555",  # Red
        "K-dwarf": "#FF9955",  # Orange
        "G-dwarf": "#FFFF55",  # Yellow (Sun-like)
        "F-dwarf": "#FFFFFF",  # White
        "A-dwarf": "#AAAAFF"   # Blue-White
    }
    
    # Add Earth's orbit for reference (1 AU)
    earth_orbit = plt.Circle((0, 0), 1.0, fill=False, color='blue', 
                           linestyle='-', alpha=0.8, linewidth=2.0, zorder=100)
    ax.add_patch(earth_orbit)
    
    # Add Mars' orbit for reference (1.52 AU)
    mars_orbit = plt.Circle((0, 0), 1.52, fill=False, color='red', 
                          linestyle='-', alpha=0.6, linewidth=1.5, zorder=100)
    ax.add_patch(mars_orbit)
    
    # Plot each star's habitable zones as concentric rings
    for i, data in enumerate(viz_data):
        # Calculate star-specific alpha to make overlapping regions visible
        alpha_base = 0.15
        
        # Draw distance-based habitable zone
        inner_edge = data['distance_based_inner_edge']
        outer_edge = data['distance_based_outer_edge']
        
        # Create rings with translucent fill
        hz_ring = plt.Circle((0, 0), outer_edge, 
                           color=star_colors[data['star_type']], 
                           alpha=alpha_base, zorder=i*2)
        ax.add_patch(hz_ring)
        
        # Add circular mask to create ring effect
        mask = plt.Circle((0, 0), inner_edge, 
                         color='black', zorder=i*2+1)
        ax.add_patch(mask)
        
        # Add outline circles
        inner_circle = plt.Circle((0, 0), inner_edge, 
                                fill=False, 
                                edgecolor=star_colors[data['star_type']], 
                                linestyle='--', 
                                linewidth=2, 
                                alpha=0.7, 
                                zorder=50)
        outer_circle = plt.Circle((0, 0), outer_edge, 
                                fill=False, 
                                edgecolor=star_colors[data['star_type']], 
                                linestyle='--', 
                                linewidth=2, 
                                alpha=0.7, 
                                zorder=50)
        ax.add_patch(inner_circle)
        ax.add_patch(outer_circle)
        
        # We no longer need the mid-point calculations since labels are removed
    
    # Add distance markers (AU indicators)
    # Create concentric dotted circles at 1 AU intervals
    for r in range(1, 6):  # Fixed range to show AU from 1 to 5
        circle = plt.Circle((0, 0), r, fill=False, linestyle=':', color='white', alpha=0.2, linewidth=0.5)
        ax.add_patch(circle)
        # Add distance label at bottom
        ax.text(0, -r, f"{r} AU", color="white", fontsize=8, ha='center', va='top', 
                bbox=dict(facecolor='black', alpha=0.7, boxstyle='round,pad=0.2'))
    
    # Add Earth and Mars labels
    # Removing these annotations as requested
    # ax.annotate("Earth (1.0 AU)", xy=(1.0, 0), xytext=(1.1, 0.3), 
    #            color='blue', fontsize=10, 
    #            bbox=dict(facecolor='black', alpha=0.7, boxstyle='round,pad=0.3'))
    
    # ax.annotate("Mars (1.52 AU)", xy=(1.52, 0), xytext=(1.6, -0.3), 
    #            color='red', fontsize=10, 
    #            bbox=dict(facecolor='black', alpha=0.7, boxstyle='round,pad=0.3'))
    
    # Set equal aspect ratio so circles look like circles
    ax.set_aspect('equal')
    
    # Set limits to show all habitable zones with some margin
    max_radius = max([data['distance_based_outer_edge'] for data in viz_data])
    ax.set_xlim(-max_radius*1.2, max_radius*1.2)
    ax.set_ylim(-max_radius*1.2, max_radius*1.2)
    
    # Add grid
    ax.grid(True, linestyle='--', alpha=0.3, color='white')
    
    # Add axis labels
    ax.set_xlabel('Distance (AU)', color='white')
    ax.set_ylabel('Distance (AU)', color='white')
    
    # Add legend for star types
    legend_elements = []
    for star_type, color in star_colors.items():
        legend_elements.append(plt.Line2D([0], [0], color=color, lw=2, label=star_type))
    
    legend_elements.append(plt.Line2D([0], [0], color='blue', lw=2, label="Earth's orbit (1 AU)"))
    legend_elements.append(plt.Line2D([0], [0], color='red', lw=1.5, label="Mars' orbit (1.52 AU)"))
    
    legend = ax.legend(handles=legend_elements, loc='upper right', facecolor='black', framealpha=0.7)
    
    # Set the legend text color to match the line colors
    for i, text in enumerate(legend.get_texts()):
        if i < len(star_colors):  # Star types
            text.set_color(list(star_colors.values())[i])
        elif i == len(star_colors):  # Earth
            text.set_color('blue')
        else:  # Mars
            text.set_color('red')
    
    # Add text explanation
    fig.text(0.02, 0.02, 
            "\nThis radial diagram shows the traditional distance-based habitable zones for different star types.\n"
            "Each colored ring represents the habitable zone for a specific type of star.\n"
            "Note how the habitable zone shifts outward for hotter, more luminous stars.\n"
            "Earth's orbit (1 AU) is shown in blue for reference.",
            fontsize=10, bbox=dict(facecolor='black', alpha=0.7), color='white')
    
    # Adjust layout
    plt.tight_layout(rect=[0, 0.04, 1, 0.96])
    
    # Save the plot
    plt.savefig(os.path.join(viz_dir, 'concentric_habitable_zones.png'), dpi=300)
    
    # Show the plot
    plt.show()

# Create a semicircular visualization showing habitable zones
def create_semicircular_visualization():
    """Create a semicircular visualization of habitable zones for different star types"""
    # Generate visualization data
    viz_data = generate_visualization_data()
    
    # Set up plot
    fig, ax = plt.subplots(figsize=(12, 6), facecolor='black')
    fig.suptitle('Semicircular Habitable Zone Diagram', fontsize=16, color='white')
    ax.set_facecolor('black')  # Set subplot background to black
    
    # Sort star types by temperature (coolest to hottest)
    viz_data.sort(key=lambda x: x['star_temp'])
    
    # Color scheme for different star types (matching the other visualizations)
    star_colors = {
        "M-dwarf": "#FF5555",  # Red
        "K-dwarf": "#FF9955",  # Orange
        "G-dwarf": "#FFFF55",  # Yellow (Sun-like)
        "F-dwarf": "#FFFFFF",  # White
        "A-dwarf": "#AAAAFF"   # Blue-White
    }
    
    # Draw semicircular habitable zones for each star type
    for i, data in enumerate(viz_data):
        star_type = data['star_type']
        
        # Use the same color scheme as in other visualizations
        color = star_colors[star_type]
        
        # Calculate inner and outer edges of habitable zone
        inner_edge = data['distance_based_inner_edge']
        outer_edge = data['distance_based_outer_edge']
        
        # Create theta values for semicircle (0 to 180 degrees)
        theta = np.linspace(0, np.pi, 100)
        
        # Create the outer and inner edges of the habitable zone
        x_outer = outer_edge * np.cos(theta)
        y_outer = outer_edge * np.sin(theta)
        
        x_inner = inner_edge * np.cos(theta)
        y_inner = inner_edge * np.sin(theta)
        
        # Plot the outer arc
        ax.plot(x_outer, y_outer, color=color, linewidth=2)
        
        # Plot the inner arc
        ax.plot(x_inner, y_inner, color=color, linewidth=2)
        
        # Fill between the arcs for the habitable zone
        ax.fill_between(
            np.concatenate([x_inner, x_outer[::-1]]),
            np.concatenate([y_inner, y_outer[::-1]]),
            color=color, alpha=0.3
        )
        
        # Add star type label at the outer edge of the habitable zone
        angle = np.pi / 2  # 90 degrees (top of the semicircle)
        label_x = outer_edge * np.cos(angle) * 1.05
        label_y = outer_edge * np.sin(angle) * 1.05
        
        ax.text(label_x, label_y, star_type, 
                color=color, fontsize=10, 
                ha='center', va='bottom',
                bbox=dict(facecolor='black', alpha=0.7, boxstyle='round,pad=0.2'))
    
    # Add Earth's orbit line for reference
    earth_theta = np.linspace(0, np.pi, 100)
    earth_x = 1.0 * np.cos(earth_theta)
    earth_y = 1.0 * np.sin(earth_theta)
    ax.plot(earth_x, earth_y, 'b--', linewidth=1.5, alpha=0.7)
    
    # Add Mars's orbit line for reference
    mars_theta = np.linspace(0, np.pi, 100)
    mars_x = 1.52 * np.cos(mars_theta)
    mars_y = 1.52 * np.sin(mars_theta)
    ax.plot(mars_x, mars_y, 'r--', linewidth=1.5, alpha=0.7)
    
    # Add distance markers (AU indicators)
    # Create concentric semicircles at 1 AU intervals
    for r in range(1, 6):  # Fixed range to show AU from 1 to 5
        # Draw semicircle
        theta = np.linspace(0, np.pi, 100)
        x = r * np.cos(theta)
        y = r * np.sin(theta)
        ax.plot(x, y, ':', color='white', alpha=0.2, linewidth=0.5)
        
        # Add distance label at the bottom
        ax.text(r, 0.1, f"{r} AU", color="white", fontsize=8, ha='center', va='bottom', 
                bbox=dict(facecolor='black', alpha=0.7, boxstyle='round,pad=0.2'))
    
    # Add horizontal axis
    ax.axhline(y=0, color='white', linestyle='-', alpha=0.3)
    
    # Set equal aspect ratio so semicircles look like semicircles
    ax.set_aspect('equal')
    
    # Set limits
    max_radius = max([data['distance_based_outer_edge'] for data in viz_data])
    ax.set_xlim(-0.5, max_radius*1.2)
    ax.set_ylim(-0.2, max_radius*1.2)
    
    # Add grid
    ax.grid(True, linestyle='--', alpha=0.3, color='white')
    
    # Add axis labels
    ax.set_xlabel('Distance (AU)', color='white')
    ax.set_ylabel('Distance (AU)', color='white')
    
    # Add legend for star types
    legend_elements = []
    for star_type, color in star_colors.items():
        legend_elements.append(plt.Line2D([0], [0], color=color, lw=2, label=star_type))
    
    legend_elements.append(plt.Line2D([0], [0], color='blue', linestyle='--', lw=1.5, label="Earth's orbit (1 AU)"))
    legend_elements.append(plt.Line2D([0], [0], color='red', linestyle='--', lw=1.5, label="Mars' orbit (1.52 AU)"))
    
    legend = ax.legend(handles=legend_elements, loc='upper right', facecolor='black', framealpha=0.7)
    
    # Set the legend text color to match the line colors
    for i, text in enumerate(legend.get_texts()):
        if i < len(star_colors):  # Star types
            text.set_color(list(star_colors.values())[i])
        elif i == len(star_colors):  # Earth
            text.set_color('blue')
        else:  # Mars
            text.set_color('red')
    
    # Add text explanation
    fig.text(0.02, 0.02, 
            "This semicircular diagram shows the habitable zones for different star types.\n"
            "Each colored band represents the habitable zone for a specific type of star.\n"
            "The semicircular format allows for easier comparison between star types.",
            fontsize=10, bbox=dict(facecolor='black', alpha=0.7), color='white')
    
    # Adjust layout
    plt.tight_layout(rect=[0, 0.04, 1, 0.96])
    
    # Save the figure
    plt.savefig(os.path.join(viz_dir, 'semicircular_habitable_zones.png'), 
                dpi=300, bbox_inches='tight', facecolor='black')
    
    # Show the plot
    plt.show()

# Create an interactive visualization comparing distance-based vs temperature-based habitable zones
def create_interactive_hz_comparison():
    """Create an interactive Plotly visualization comparing distance-based and temperature-based habitable zones"""
    # Generate data
    viz_data = generate_visualization_data()
    
    # Create a new figure
    fig = go.Figure()
    
    # Add a black background
    fig.update_layout(
        paper_bgcolor="black",
        plot_bgcolor="black",
        title={
            'text': "Interactive Comparison: Distance vs Temperature-Based Habitable Zones",
            'y':0.95,
            'x':0.5,
            'xanchor': 'center',
            'yanchor': 'top',
            'font': {'color': 'white', 'size': 24}
        },
        width=1000,
        height=1000,
        xaxis=dict(
            title="Distance (AU)",
            gridcolor='rgba(255, 255, 255, 0.1)',
            zerolinecolor='rgba(255, 255, 255, 0.2)',
            tickfont=dict(color='white'),
            titlefont=dict(color='white', size=16)
        ),
        yaxis=dict(
            title="Distance (AU)",
            gridcolor='rgba(255, 255, 255, 0.1)',
            zerolinecolor='rgba(255, 255, 255, 0.2)',
            scaleanchor="x",
            scaleratio=1,
            tickfont=dict(color='white'),
            titlefont=dict(color='white', size=16)
        )
    )
    
    # Color scheme for different star types
    star_colors = {
        "M-dwarf": "#FF5555",  # Red
        "K-dwarf": "#FF9955",  # Orange
        "G-dwarf": "#FFFF55",  # Yellow (Sun-like)
        "F-dwarf": "#FFFFFF",  # White
        "A-dwarf": "#AAAAFF"   # Blue-White
    }
    
    # Sort data by star temperature (coolest to hottest)
    viz_data.sort(key=lambda x: x['star_temp'])
    
    # Add Earth's orbit for reference
    theta = np.linspace(0, 2*np.pi, 100)
    earth_x = 1.0 * np.cos(theta)
    earth_y = 1.0 * np.sin(theta)
    
    fig.add_trace(go.Scatter(
        x=earth_x,
        y=earth_y,
        mode='lines',
        line=dict(color='blue', width=2),
        name="Earth's orbit (1 AU)",
        hoverinfo='name'
    ))
    
    # Add Mars' orbit for reference
    mars_x = 1.52 * np.cos(theta)
    mars_y = 1.52 * np.sin(theta)
    
    fig.add_trace(go.Scatter(
        x=mars_x,
        y=mars_y,
        mode='lines',
        line=dict(color='red', width=2),
        name="Mars' orbit (1.52 AU)",
        hoverinfo='name'
    ))
    
    # Add traces for each star type
    for data in viz_data:
        star_type = data['star_type']
        color = star_colors[star_type]
        
        # Distance-based habitable zone (outer ring)
        theta = np.linspace(0, 2*np.pi, 100)
        
        # Create the outer ring
        outer_x = data['distance_based_outer_edge'] * np.cos(theta)
        outer_y = data['distance_based_outer_edge'] * np.sin(theta)
        
        fig.add_trace(go.Scatter(
            x=outer_x,
            y=outer_y,
            mode='lines',
            line=dict(color=color, width=2, dash='dot'),
            name=f"{star_type} Distance HZ (Outer: {data['distance_based_outer_edge']:.2f} AU)",
            hoverinfo='name',
            hoverlabel=dict(bgcolor=color)
        ))
        
        # Create the inner ring
        inner_x = data['distance_based_inner_edge'] * np.cos(theta)
        inner_y = data['distance_based_inner_edge'] * np.sin(theta)
        
        fig.add_trace(go.Scatter(
            x=inner_x,
            y=inner_y,
            mode='lines',
            line=dict(color=color, width=2, dash='dot'),
            name=f"{star_type} Distance HZ (Inner: {data['distance_based_inner_edge']:.2f} AU)",
            fill='tonext',
            fillcolor=f'rgba({int(color[1:3], 16)}, {int(color[3:5], 16)}, {int(color[5:7], 16)}, 0.2)',
            hoverinfo='name',
            hoverlabel=dict(bgcolor=color)
        ))
        
        # Temperature-based habitable zone if available
        if data['temp_based_inner_edge'] is not None:
            # Create the outer ring
            temp_outer_x = data['temp_based_outer_edge'] * np.cos(theta)
            temp_outer_y = data['temp_based_outer_edge'] * np.sin(theta)
            
            fig.add_trace(go.Scatter(
                x=temp_outer_x,
                y=temp_outer_y,
                mode='lines',
                line=dict(color=color, width=2),
                name=f"{star_type} Temp HZ (Outer: {data['temp_based_outer_edge']:.2f} AU)",
                hoverinfo='name',
                hoverlabel=dict(bgcolor=color)
            ))
            
            # Create the inner ring
            temp_inner_x = data['temp_based_inner_edge'] * np.cos(theta)
            temp_inner_y = data['temp_based_inner_edge'] * np.sin(theta)
            
            fig.add_trace(go.Scatter(
                x=temp_inner_x,
                y=temp_inner_y,
                mode='lines',
                line=dict(color=color, width=2),
                name=f"{star_type} Temp HZ (Inner: {data['temp_based_inner_edge']:.2f} AU)",
                fill='tonext',
                fillcolor=f'rgba({int(color[1:3], 16)}, {int(color[3:5], 16)}, {int(color[5:7], 16)}, 0.3)',
                hoverinfo='name',
                hoverlabel=dict(bgcolor=color)
            ))
    
    # Add explanations as annotations
    fig.add_annotation(
        x=0.02,
        y=0.02,
        xref="paper",
        yref="paper",
        text="<b>Temperature-based</b> habitable zones: solid lines<br><b>Distance-based</b> habitable zones: dotted lines<br>Hover over elements for details",
        showarrow=False,
        font=dict(color="white", size=14),
        align="left",
        bgcolor="rgba(0,0,0,0.7)",
        bordercolor="white",
        borderwidth=1,
        borderpad=4
    )
    
    # Add buttons for toggling visibility
    updatemenus = [
        dict(
            type="buttons",
            direction="left",
            buttons=[
                dict(
                    args=[{'visible': [True] * len(fig.data)}],
                    label="Show All",
                    method="update"
                ),
                dict(
                    args=[{'visible': [True, True] + [trace.name and 'Distance' in trace.name for trace in fig.data[2:]]}],
                    label="Distance-based Only",
                    method="update"
                ),
                dict(
                    args=[{'visible': [True, True] + [trace.name and 'Temp' in trace.name for trace in fig.data[2:]]}],
                    label="Temperature-based Only",
                    method="update"
                ),
            ],
            pad={"r": 10, "t": 10},
            showactive=True,
            x=0.5,
            xanchor="center",
            y=1.1,
            yanchor="top",
            bgcolor="rgba(0,0,0,0.7)",
            font=dict(color="white")
        ),
    ]
    
    fig.update_layout(updatemenus=updatemenus)
    
    # Save the interactive HTML plot
    pio.write_html(fig, file=os.path.join(viz_dir, 'interactive_habitable_zones.html'), auto_open=False)
    
    print("Interactive visualization saved to 'Visualisations/interactive_habitable_zones.html'")
    
    return fig

# Modified main function to include all visualizations
if __name__ == "__main__":
    print("Creating habitable zone visualizations...")
    # Create basic visualization (temperature vs. distance for each star type)
    create_habitability_visualization()
    
    # The new semicircular visualization replaces the old concentric circles visualization
    create_semicircular_visualization()
    
    # We're only showing the two main visualizations to avoid duplication
    # create_advanced_visualization_with_planets()
    
    # Create interactive visualization
    create_interactive_hz_comparison()
    
    print("Visualizations completed and saved to 'Visualisations' folder.")