# glfBLUP-Procrustes

<img width="907" height="473" alt="graphical_abstract" src="https://github.com/user-attachments/assets/13b7e0b0-867d-41d7-9be0-33dc01bf461f" />

Scripts for generating results in Kunst et al. (xxxx) A latent factor approach to hyperspectral time series data for multivariate genomic prediction of grain yield in wheat.

# Folder Structure
Within the _scripts_ folder you will find the R scripts used for generating the results and figures found in the paper. The _workflow_ folder has numbered scripts to run the sequence of steps for generating factor scores and performing genomic predictions. The data generated through the scripts in the _workflow_ folder can be visualised using the scripts in the _visualisation_ folder. 

Within the _functions_ folder you will find auxillary functions used within the scripts. 

I have added an R project file which should make it possible to run the scripts and function when the data used is available. For getting the data, I refer to [Melsen et al. (2025)](https://doi.org/10.1002/bimj.70081) *Improving Genomic Prediction Using High-Dimensional Secondary Phenotypes: The Genetic Latent Factor Approach* and [Krause et al. (2019)](https://doi.org/10.1534/g3.118.200856) *Hyperspectral Reflectance-Derived Relationship Matrices for Genomic Prediction of Grain Yield in Wheat*.

For the environment names the following abbreviations are used:
Optimal flat: B5IR
Moderate drought: B2IR
Severe drought: DRIP

