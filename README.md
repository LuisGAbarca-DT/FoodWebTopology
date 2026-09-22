#  
#  
#               Food Web Topology Analysis
#   
#                     V. 1.2
**12 SEPTEMBER 2026**

**Authors:**

- Luis Gerardo Abarca     gabarca@uv.mx , luisgaa@gmail.com
- Israel Huesca Domínguez ihuesca@uv.mx
   
**Institution:**

- Instituto de Investigaciones Biológicas 
    - Universidad Veracruzana
    - Veracruz, México

**Repository:** https://github.com/LuisGAbarca-DT/FoodWebTopology

## Overview

### Computes various structural and topological indices for a food web, including:
    - Number of species by category (basal, intermediate, top)
    - Connectivity and centrality values
    - Trophic levels
    - Modularity (using Leiden algorithm with adjustable resolution)
    - And many other topological metrics
    - Multiple Null Models: Erdős–Rényi, Cascade, Niche (two versiona), Randomm links (Cannonball) 
            - Every random food web is analyzed through vaious indices
    - Validation Checks: Automated data integrity verification
    - Reproducible: Complete analysis logging and parameter tracking
    - All results are exported as comma-delimited files.



## Data Format Requirements

- The adjacency matrix should be in csv format with:
- **Columns** = predators
- **Rows** = prey  
- **Values**: 0 (no interaction) or 1 (predation)
- **Requirements**:
  - A CSV file with a square adjacency matrix (same number of rows and columns)
  - Square matrix with the same species order in rows and columns
  - At least one basal species (in-degree = 0)
  - At least one top species (out-degree = 0)
  - No isolated nodes or disconnected groups
  - Nodes names should be alpha-numeric
        - If the row and column names are numeric, the program will add the suffix "SPS_" to the number in order to comply with the alphanumeric characteristic

  
### **Validation:**
- The script will stop with a clear error message if your file does not meet 
these specifications. Please reformat your data accordingly. A template file is 
provided in the repository (`template_foodweb.csv`).
  
### **Example: `example_foodweb.csv`**

- ,Species_A,Species_B,Species_C,Species_D
- Species_A,0,1,0,0
- Species_B,0,0,1,1
- Species_C,0,0,0,0
- Species_D,0,0,1,0

## Quick Start


### 1. Setup
- Download all fwt files to a single folder
- Place your food web csv file in the same folder (Optional)
- Set R's working directory to this folder
- The results will be saved in a subdirectory named after your csv file (without the .csv extension) in this directory



### 2. Run Analysis

    source("fwt.R")```
    
- This will install all dependencies (if not already installed) and load the libraries

### 3. Follow Prompts

- Select the csv file with the food web data
- Choose the randomization model
- Enter the number of random food webs to be generated
- Enter if you want to compute number of chains and trophic levels for each
    random web (the computation is time consuming and uses big amounts of RAM for large webs)
- Enter the resolution for the Leiden Modularity algorithm (higher = more modules)
- Enter if you want a figure of the food web and each of the modules. 
    
### 4. Results

- Considering the data file name, the resuts will be saved according with this
name as a prefix of the type of analysis presented. For example, if the file 
name of the data is **your_web.csv** then the program will generate a directory (**RESULTS**). 
Within it, two more directories: **ORIGINAL** and **RANDOM** The first one with the results for the 
original data and the second one with the results for the randomized webs. A log file
will be placed at the ORIGINAL directory with information related to the analysis.



### Support

For issues or questions, please open an issue on GitHub or contact the authors.

