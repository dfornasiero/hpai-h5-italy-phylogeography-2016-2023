# Continuous phylogeography reveals shifting environmental drivers of Highly Pathogenic Avian Influenza H5 spread in Italy, 2016-2023
Diletta Fornasiero¹˒⁶\*, Paolo Mulatti¹, Alice Fusaro¹, Isabella Monne¹, Fabiana Gambaro², Simon Dellicour²˒³˒⁴, Francesca Scolamacchia¹, Luca Martelli¹, Giulia Faustini⁵, Mariette F. Ducatez⁶, Claire Guinat⁶

¹ Istituto Zooprofilattico Sperimentale delle Venezie, Viale dell’Università 10, 35020, Legnaro, Italy\
² Spatial Epidemiology Lab (SpELL), Université Libre de Bruxelles, Brussels, Belgium\
³ Department of Microbiology, Immunology and Transplantation, Rega Institute, Laboratory for Clinical and Epidemiological Virology, KU Leuven, Leuven, Belgium\
⁴ Interuniversity Institute of Bioinformatics in Brussels, Université Libre de Bruxelles, Vrije Universiteit Brussel, Brussels, Belgium\
⁵ Istituto Zooprofilattico Sperimentale della Lombardia ed Emilia-Romagna, Via A. Bianchi 9, 25124, Brescia, Italy\
⁶ Univ Toulouse, ENVT, INRAE, IHAP, Toulouse, France

\*Corresponding author: **Diletta Fornasiero** - `dfornasiero@izsvenezie.it`
------------------------------------------------------------------------

## Abstract
<small> The ongoing global spread of highly pathogenic avian influenza (HPAI) continues causing major impacts on poultry, wildlife, and public health. Among the most affected countries in Europe, Italy experienced three major epidemic waves in 2016–2023, caused by 2.3.4.4b H5N8 and H5N1 subtypes. Yet, the underlying spatio-temporal dispersal dynamics and drivers of these epidemics remain to be elucidated. Here, we combined viral genome sequencing with continuous phylogeographic inference to reconstruct the evolutionary and spatio-temporal history of these epidemics. By combining genetic data with landscape, poultry-related variables, and wild bird abundances, we explored how the environment was associated with viral diffusion in space and time.
We identified distinct diffusion patterns and environmental associations across the three epidemics. The first two waves (2016–2017 and 2021–2022) were dominated by local poultry-to-poultry transmission in high-density farming areas, with agricultural areas being significantly associated with a relatively higher viral lineage dispersal velocity in 2021–2022. In contrast, the most recent wave (2022–2023) showed evidence of mid- to long-distance dispersal events positively associated with wetlands, waterbodies, and gull abundances, pointing to a larger role of wild birds in HPAI spatial dynamics. The estimated weighted diffusion coefficients and spatial wavefront distances revealed a shift from localised poultry-driven to long-distance spread, likely related to wild birds.
Our findings highlight a progressive shift from poultry-driven epidemics to more complex environmental transmission dynamics involving wild birds and natural habitats. These changes underscore the need to adapt surveillance and control strategies to an evolving and ecologically diverse viral landscape. </small>

------------------------------------------------------------------------

# Spatiotemporal Dynamics of HPAI in Italy (2016-2023)

This repository contains the data and analytical pipeline for reconstructing the phylogeographic dynamics and environmental drivers of HPAI epidemics in Italy. The workflow integrates discrete and continuous BEAST analyses with phylogeography analyses using the [seraphim](https://github.com/sdellicour/seraphim) R package.

---

## 📂 Project Structure
The project is organised into five main functional modules:

* **`data/`**: Cleaned environmental rasters per epidemic wave, and shapefiles for Italy and Europe.
* **`data_preparation/`**: Scripts for sequence cleaning and Maximum Likelihood (ML) tree pre-processing.
* **`DTA/`**: BEAST XML configurations for Discrete Trait Analysis.
* **`continuous_reconstruction/`**: BEAST XML configurations for Continuous Phylogeography.
* **`phylogeography/`**: A set of scripts for tree extraction, dispersal statistics, and phylogeography analyses.

---

## 🚀 Analytical Pipeline
To replicate the analysis, follow the scripts in this order:

### 1. Data Pre-processing
Clean and align sequences, followed by ML tree estimation to identify monophyletic clades.
* `HPAI_sequence_preparation.R`
* `ML_data_preparation.R`

### 2. Bayesian Inference (BEAST)
Run the XML files located in the `DTA/` and `continuous_reconstruction/` folders using **BEAST v1.10.4**.
> [!IMPORTANT]
> **Data Privacy Note:** To comply with privacy regulations, precise geographic coordinates in the continuous phylogeography XML files have been removed. Consequently, these specific XMLs are provided for structural reference of the model parameters and cannot be re-run.

### 3. Posterior Tree Extraction
Extract 1,000 random trees from the posterior distribution and filter for the study area.
* `beast_output_data_extraction_1.R`: For single monophyletic group analysis (i.e., 2016-2017 epidemic).
* `beast_output_data_extraction_2.R`: For multiple monophyletic group scenarios (i.e., 2021-2022 and 2022-2023 epidemics).

### 4. Dispersal Statistics
Estimate diffusion coefficients and spatial expansion over time.
* `phylogeographic_dispersal_statistic.R`: Core dispersal metrics (weighted diffusion coefficient and spatial wavefront distance).
* `spatial_wavefront_distance_per_monophyletic_group.R`: To compute multiple spatial wavefront distances, one for each monophyletic group.

### 5. Landscape Genetics (Environmental Correlation)
Test the impact of environmental factors (Conductance/Resistance) on viral dispersal.
* `impact_dispersal_location.R`: Correlation with specific locations.
* `impact_dispersal_velocity.R`: Using the Least-Cost and Circuitscape path models.
* `impact_dispersal_velocity_processing.R`: Result synthesis and final table generation for dispersal velocity.

---

## 📊 Visualizations
## 📊 Visualizations
The final results are integrated and rendered using the `multipanel_plots.R` script. This script synthesises the outputs from all previous steps to generate the primary publication figures (Figures 1-3), consisting of:
* **Panel A**: **Epidemic Curves**
* **Panel B**: **MCC Continuous Trees**
* **Panel C**: **Spatial Wavefront Distance**
* **Panel D**: **Continuous Spatial Reconstruction**


---

## 🛠 Requirements
* **Software**: R, BEAST v1.10.x, FigTree
* **R Packages**: `seraphim`, `raster`, `sf`, `ggplot2`, `lubridate`, `dplyr`

---

## 💡 Usage Note
When adapting these scripts for different epidemic waves (e.g., 2016-2017, 2021-2022, 2022-2023), ensure you update the `mostRecentSamplingDatum` and the working directory paths in the configuration section of each script.
