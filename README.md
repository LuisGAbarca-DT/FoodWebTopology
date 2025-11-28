#  
#  
#               Food Web Topology Analysis
#   
**27 NOVEMBER 2025**

**Authors:**

- Luis Gerardo Abarca     gabarca@uv.mx   luisgaa@gmail.com
- Israel Huesca Domínguez ihuesca@uv.mx
   
**Institution:**

- Instituto de Investigaciones Biológicas 
- Universidad Veracruzana

**Repository:** https://github.com/LuisGAbarca-DT/FoodWebTopology

## Overview

Computes various structural and topological indices for a food web, Including
- Number of species by category (basal, intermediate, top)
- Connectivity and centrality values
- Trophic levels
- Modularity (using Leiden algorithm
- And other topological metrics

Generates random food webs through different null models. All results are 
exported as comma-delimited files.


## Data Format Requirements

- The adjacency matrix should be in csv format with:
- **Columns** = predators
- **Rows** = prey  
- **Values**: 0 (no interaction) or 1 (predation)
- **Requirements**:
  - Square matrix (same species order in rows and columns)
  - At least one basal species (in-degree = 0)
  - At least one top species (out-degree = 0)
  - No isolated nodes or disconnected groups

## Quick Start


### 1. Setup
- Download all FWTopo files to a single folder
- Place your food web csv file in the same folder
- Set R's working directory to this folder


### 2. Run Analysis

    source("run_FWTopo.R")```
    
- This will install all dependencies and load the libraries

### 3. Follow Prompts

- Select the csv file with the food web data
- Choose the randomization model
- Enter the number of random food webs to be generated
- Enter the resolution for the Leiden Modularity algorithm (higher = more modules)
    

Considering the data file name, the resuts will be saved according with this
name as a prefix of the type of analysis presented. For example, if the file 
name of the data is **your_web.csv** then the program will genrate a directory 
with that name and inside it two more directories, one with the results for the 
original data and one more with the results for the randomized webs.

### Features
- Multiple Null Models: Cascade, Niche, Erdős–Rényi
- Comprehensive Metrics: 15+ structural and topological indices
- Modularity Analysis: Leiden algorithm with adjustable resolution
- Validation Checks: Automated data integrity verification
- Reproducible: Complete analysis logging and parameter tracking

### Support

For issues or questions, please open an issue on GitHub or contact the authors.

