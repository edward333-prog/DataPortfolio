# Exoplanet Habitability Analysis

## Project Overview
This project analyzes exoplanet data from NASA's Kepler mission to develop and refine methodologies for assessing exoplanet habitability. Our analysis evolves beyond traditional distance-based models to incorporate temperature-based assessment and star-specific evolutionary considerations.

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
![Habitable Zone Visualization](./Visualizations/goldilocks_zone.png)

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
