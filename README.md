# Sediment Mixing Analysis - Mataquito River

This repository contains the analytical framework and computational notebooks for sediment mixing analysis in the Mataquito River system, supporting research on sediment provenance and transport dynamics in the Chilean Coastal Cordillera.

## Project Overview

This project investigates sediment mixing processes through multi-proxy analysis combining detrital zircon geochronology, cosmogenic nuclide analysis, and geomorphological assessment. The analysis integrates fertility calculations, spatial watershed analysis, and mixing model scenarios to quantify sediment provenance and transport mechanisms along the Mataquito River system.

### Research Objectives
- Quantify sediment provenance and mixing processes in the Mataquito River system
- Determine fertility ratios for quartz and zircon minerals from different source regions
- Assess erosion rates and sediment transport dynamics using ¹⁰Be cosmogenic nuclides
- Evaluate mixing models and their uncertainties through statistical analysis
- Characterize geologic controls on sediment production and transport

### Study Area
The Mataquito River system in central Chile (35°S) drains from the Andes Mountains to the Pacific Ocean, crossing multiple geological terranes including the Main Cordillera, Central Depression, and Coastal Cordillera. Sample locations range from high-elevation tributaries (>600m) to the river mouth near the coast.

## Repository Structure

### Core Analysis Notebooks

#### Geomorphological Analysis
- **`MataquitoAll30mDEManalysis.ipynb`** - Digital elevation model analysis of the Mataquito River system using 30m resolution DEM, including drainage network extraction, watershed delineation, and topographic metrics calculation for 11 sample locations (CT-1 through CT-11)

#### Cosmogenic Nuclide Analysis
- **`ProductionRates.ipynb`** - Calculates ¹⁰Be production rates for cosmogenic nuclide analysis and erosion rate quantification using sample-specific shielding factors and topographic scaling

#### Fertility Analysis
- **`Mataquito/FertilityAnalysis/Mataquito_quartz_zircon_fertilities.ipynb`** - Calculates fertility indices for quartz and zircon populations using Monte Carlo simulations (100,000 iterations) to assess source rock contributions and mineral abundances across different mixing scenarios

#### Geological Analysis
- **`Mataquito/Geology/SubwatershedGeologicUnits.ipynb`** - Quantifies areal distributions and proportions of geological units within subwatersheds of the Mataquito River system using 1:1,000,000 scale geological mapping

#### Detrital Zircon Mixing Analysis
- **`Mataquito/Zircon/detritalPy_mix/detritalMixer_CT_2025.ipynb`** - Calculates mixing coefficients for different sediment source scenarios and evaluates mixing model performance using detritalPy-mix version 0.1 (beta) with bootstrapping analysis (5,000 iterations)

- **`Mataquito/Zircon/detritalPy/detritalPy_Mataquito_All_2025.ipynb`** - Comprehensive detrital zircon analysis using detritalPy framework v1.4.3, including age distributions, provenance analysis, multi-dimensional scaling, and statistical visualization for all Mataquito samples

#### Downstream Mixing Analysis
- **`Mataquito/Mixing_below_confluence/Below_confluence.ipynb`** - Analyzes downstream evolution of sediment mixing below tributary confluences using both detrital zircon and ¹⁰Be data to understand transport length scales

### Data Files

#### Sample Data
- **`Mataquito/MataquitoSampleData.xlsx`** - Master sample database containing coordinates, watershed areas, cosmogenic nuclide concentrations, erosion rates, and analytical uncertainties for 11 samples
- **`Mataquito/mataquito_area_km2.xlsx`** - Watershed area calculations

#### Geological Data
- **`Mataquito/Geology/subwatershed_geologic_units_data.xlsx`** - Geological unit distributions within each subwatershed
- **`Mataquito/Geology/combined_subwatershed_summary.xlsx`** - Compiled geological statistics for all subwatersheds

#### Detrital Zircon Data
- **`Mataquito/Zircon/Spreadsheets/`** - U-Pb geochronological data for samples CT-1 through CT-11, including individual analytical sessions and compiled datasets
  - `CT-1-11_ZrUPb_datasheets.xlsx` - Master compilation of all detrital zircon U-Pb analyses
  - Individual sample files from multiple analytical sessions

#### Mixing Model Results
- **`Mataquito/Zircon/detritalPy_mix/Spreadsheets/`** - Bootstrap mixing coefficient results (5,000 iterations each) for various parent-child mixing scenarios
- **`Mataquito/Zircon/detritalPy_mix/Figures/`** - Visualization outputs including age distribution plots and mixing result summaries

### Analytical Outputs

#### Figures
- **`Mataquito/FertilityAnalysis/Figures/`** - Violin plots showing fertility ratio distributions for different mixing scenarios
- **`Mataquito/Mixing_below_confluence/Figures/`** - Downstream mixing evolution plots for ¹⁰Be concentrations and erosion rates
- **`Mataquito/Zircon/detritalPy/Figures/`** - Detrital zircon age distribution plots, multi-dimensional scaling results, and provenance visualizations
- **`Mataquito/MataquitoManuscriptFig1.png`** - Sample location map for reference

## Methodology

### Sampling Strategy
Eleven sediment samples (CT-1 through CT-11) were collected from strategic locations throughout the Mataquito River system to capture different tributary contributions and downstream mixing processes. Sample locations were chosen to:
- Represent major tributaries from different geological terranes
- Capture mixing processes at tributary confluences
- Document downstream transport and homogenization

### Analytical Methods

#### Detrital Zircon Geochronology
- U-Pb analyses conducted via laser ablation inductively coupled plasma mass spectrometry (LA-ICP-MS)
- Data processing using detritalPy v1.4.3 for age distribution analysis, provenance assessment, and statistical comparisons
- Mixing model analysis using detritalPy-mix v0.1 with bootstrap uncertainty quantification (5,000 iterations)

#### Cosmogenic Nuclide Analysis
- ¹⁰Be concentrations measured in quartz separates
- Production rate scaling using topographic shielding factors and virtual elevation calculations
- Erosion rate determination incorporating analytical and production rate uncertainties

#### Fertility Calculations
- Monte Carlo simulations (100,000 iterations) incorporating erosion rate uncertainties
- Flux-ordered mixing scenarios to ensure physically realistic combinations
- Separate fertility calculations for quartz and zircon using different analytical frameworks

#### Geospatial Analysis
- Watershed delineation using 30m digital elevation models
- Geological unit area calculations using 1:1,000,000 scale mapping
- Distance calculations for downstream transport analysis

## Software Dependencies

### Required Python Packages
- **detritalPy v1.4.3** - Detrital zircon data analysis and visualization
- **detritalPy-mix v0.1** - Mixing model analysis with uncertainty quantification
- **numpy, pandas** - Data manipulation and numerical analysis
- **matplotlib, seaborn** - Visualization and statistical plotting
- **geopandas, rasterio** - Geospatial data processing
- **geopy** - Distance calculations
- **scipy** - Statistical analysis

### External Software
- **TopoAnalysis** - DEM processing and watershed analysis (custom module)
- **Microsoft Excel** - Data compilation and storage

## Data Processing Workflow

1. **Sample Collection and Preparation**
   - Field sampling at 11 strategic locations
   - Laboratory processing for zircon separation and quartz purification

2. **Analytical Data Generation**
   - U-Pb geochronology of detrital zircons
   - ¹⁰Be cosmogenic nuclide analysis

3. **Geospatial Data Processing**
   - Watershed delineation from 30m DEM
   - Geological unit area calculations
   - Sample location coordinate verification

4. **Statistical Analysis**
   - Bootstrap mixing model analysis (5,000 iterations)
   - Monte Carlo fertility calculations (100,000 iterations)
   - Uncertainty propagation through all calculations

5. **Results Visualization**
   - Age distribution plots and provenance diagrams
   - Mixing coefficient violin plots
   - Downstream transport evolution plots

## Key Results

### Mixing Scenarios
The analysis evaluates multiple mixing scenarios including:
- 2-component mixtures (e.g., CT-5 + CT-6 → various downstream samples)
- 3-component mixtures (e.g., CT-5 + CT-6 + CT-3 → CT-9)
- Downstream evolution of mixing proportions

### Fertility Ratios
Calculated fertility ratios show:
- Significant differences between quartz and zircon fertility
- Systematic variations related to source rock geology
- Uncertainty ranges derived from Monte Carlo analysis

### Transport Dynamics
Downstream analysis reveals:
- Length scales of sediment homogenization
- Erosion rate variations along transport pathways
- Systematic changes in mixing proportions with distance

## Citation and Usage

This repository supports ongoing research on sediment mixing processes in the Chilean Coastal Cordillera. If you use any data or methods from this repository, please cite the associated publications.

### Key Publications
*[Publications list to be updated upon manuscript submission]*

### Software Citations
- Sharman, G.R., and Johnstone, S.A., 2017, Sediment unmixing using detrital geochronology: Earth and Planetary Science Letters, v. 477, p. 183-194.
- Malkowski, M.A., Sharman, G.R., Johnstone, S.A., and Kimbrough, D.L., 2019, Dilution and propagation of provenance signals in siliciclastic sediments: American Journal of Science, v. 319, p. 773-803.

## Contact Information

For questions regarding data access, methodology, or collaboration opportunities, please contact the repository maintainers.
