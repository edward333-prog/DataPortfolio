# Exoplanet Habitability Analysis - "Are We Alone In The Universe?" 

"As we refine our understanding of habitability beyond simple distance-based models, we're discovering that the cosmic real estate suitable for life may be far more extensive than previously thought. This project explores that expanded frontier."
                    
## Project Overview
This project analyzes exoplanet data from NASA's Kepler mission to develop and refine methodologies for assessing exoplanet habitability. Our analysis evolves beyond traditional distance-based models to incorporate temperature-based assessment and star-specific evolutionary considerations. With two potential candidates returned in the results at a 0.02% success rate. This is exciting as I was expecting 0%. 

## Key Innovations

### Temperature-Based Habitability Assessment
We've shifted from conventional distance-based habitable zone calculations to a more nuanced temperature-based model that:
- Defines a habitable temperature range of 273K-323K (0°C-50°C) for complex life
- Resolves contradictions between high Earth-similarity scores and traditional habitable zone boundaries
- Reveals that planets can maintain Earth-like temperatures despite being outside distance-based habitable zones

### Star-Type Specific Evolution Models
Our models account for different stellar evolution patterns:
- M-dwarf: ~5% luminosity increase per billion years
- K-dwarf: ~8% luminosity increase per billion years
- G-dwarf: ~10% luminosity increase per billion years (scaled by mass^0.4)

### Tidal Locking Analysis
We consider how the terminator zone (boundary between day/night sides) of tidally locked planets can provide habitable conditions even when other regions cannot.

## Data & Methodology
This analysis uses data from NASA's Exoplanet Archive, focusing on key parameters:
- Planetary radius (koi_prad)
- Equilibrium temperature (koi_teq)
- Insolation flux (koi_insol)
- Orbital period (koi_period)
- Stellar characteristics (koi_steff, koi_slogg, koi_srad)

[View Complete Methodology Document](./METHODOLOGY.md)

## Key Findings
1. **Temperature Trumps Distance**: Planets can maintain Earth-like temperatures despite being outside traditional habitable zones
2. **Star-Specific Habitable Zones**: Different star types require tailored habitability calculations
3. **Terminator Zone Potential**: The boundary regions of tidally locked planets offer unique habitability potential

## Visualizations

These next two visualisations are part of our Earth Similarity Index (ESI) 

These radar charts visualize our innovative temperature-based approach to exoplanet habitability assessment. Unlike traditional distance-based calculations, our methodology evaluates planets based on actual physical properties that contribute to Earth-like conditions.

Key Categories in the Radar Charts:

Radius: Similarity to Earth's radius (1.0 = exact match), affecting surface gravity and potential for maintaining a suitable atmosphere.
Temperature: Correspondence to Earth's habitable temperature range (273K-323K or 0°C-50°C), which we prioritize as critical for supporting complex life.
Insolation: Similarity to the solar energy Earth receives from the Sun, accounting for different star types and their evolution over time.
Period: Orbital period similarity to Earth, indicating suitable day-night cycles which affect climate stability.
Star Type: Suitability of the host star for supporting habitable conditions. We've developed specific models for:
M-dwarf stars (~5% luminosity increase per billion years)
K-dwarf stars (~8% luminosity increase per billion years)
G-dwarf stars (~10% luminosity increase per billion years)
Escape Velocity: Similarity to Earth's escape velocity, indicating a planet's ability to retain an atmosphere necessary for life.

These visualizations demonstrate how planets like K05499.01 (98.15% habitability) and K03497.01 (96.69% habitability) maintain Earth-like temperatures despite potentially lying outside traditional habitable zones, supporting our key finding that temperature-based assessment provides a more accurate picture of habitability than distance-based calculations alone.

- K05499.01 Radar Chart
<img src="https://raw.githubusercontent.com/YOUR-USERNAME/ExoPlanet/main/Visualisations/K03497.01_radar.png" alt="K05499.01 Radar Chart">

- K03497.01 Radar Chart 
<picture>
  <source srcset="/ExoPlanet/Visualisations/K03497.01_radar.png" media="(prefers-color-scheme: dark)">
  <img src="/ExoPlanet/Visualisations/K03497.01_radar.png" alt="K03497.01 Radar Chart">
</picture>

## Technologies Used
- Python (pandas, numpy)
- Data processing pipeline for NASA Exoplanet Archive datasets
- Custom habitability modeling algorithms

## Repository Structure
- `/Data`: Contains raw and processed exoplanet datasets
- `/Scripts`: Python scripts for data cleaning and analysis
- [METHODOLOGY.md](./METHODOLOGY.md): Comprehensive methodology documentation
- [Other documents]

## Future Directions
- Atmospheric composition and greenhouse effect modeling
- Planetary magnetic field analysis
- Geological activity considerations
- Non-Earth-like biochemistry exploration

## Download Full Reports
- [Comprehensive Analysis Report](./Reports/analysis.pdf) [Placeholder]
- [Detailed Methodology Evolution](./METHODOLOGY.md)
- [Data Processing Documentation](./Documentation/data_processing.pdf) [Placeholder]
