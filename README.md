[README_FWTopo.md](https://github.com/user-attachments/files/22982299/README_FWTopo.md)
#  
#  
#               Food Web Topology Analysis
#   
#                   17 OCTOBER 2025
#   
#   Luis Gerardo Abarca     gabarca@uv.mx   luisgaa@gmail.com
#   Israel Huesca Domínguez ihuesca@uv.mx
#   
#   Instituto de Investigaciones Biológicas 
#        Universidad Veracruzana
#   

Computes various structural and topological indices for a food web, like number 
of species by category (basal, intermediate, top), connectivity, centrality 
values, trophic level, modularity, and others.

Calculates modularity using the Leiden algorithm and generate random food webs 
through different models  There is the possibility to choose the modularity 
algorithm. There is the possibility to choose the algorithm to produce a desire 
number of food webs randomly by choosing a series of different algorithms. All
results are saved as a comma-limited file. For the original as well as for all 
new randomly produced food webs structural properties are computed and saved 
as comma-delimited file. For each new food web the number of modules and the
modularity is computed and saved. 

The adyacency matrix should be in csv format (separated by commas) with columns 
as predators and rows as prey. The format is 0 with no prey eaten by a 
predator and 1 if a predator eats the prey.

1. Setup
    Download all FWTopo files to a single folder
    Place your food web CSV file in the same folder
    Set Working Directory as the folder containing FWTopo files

2. Run Analysis
    In R or RStudio console:
    source("run_FWTopo.R")
    This will install all dependencies and load the libraries

3. Follow Prompts
    Select the csv file with the food web data
    Choose the randomization model
    Enter the number of random food webs to be generated
    Enter the resolution for the Leiden Modularity algorithm (if greater than 1 
    more modules)
    
Data Characteristics
    Rows are the prey
    Columns are the predators
    Same species order in rows and columns
    Data file must contain:
        One or more basal species (in-degree = 0)
        One or more top species (out-degreeo = 0)
        No isolated nodes or groups
        
Considering the data file name, the resuts will be saved according with this
name as a prefix of the type of analysis presented. For example, if the file 
name of the data is your_we.csv then the program will genrate a directory with 
that name and inside it two more directories, one with the results for the 
original data and one more with the results for the randomized webs.Following 
is an example of the structure of the directories:

your_web_RESULTS/                  # Main results directory
your_web_ANALYSIS_LOG.txt          # Parameters & validation report
├── ORIGINAL/                      # Original food web analysis
│   ├── your_web_ORIGINAL_STRC.csv
│   ├── your_web_ORIGINAL_MODULES.csv
│   └── ...
└── RANDOM/                        # Random web ensemble analysis
    ├── your_web_CASCADE_stats.csv
    ├── your_web_NICHE_stats.csv
    └── ...
    
    

