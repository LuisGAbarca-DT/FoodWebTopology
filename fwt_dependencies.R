#  
#  
#               Food Web Topology Analysis
#                   Dependencies Script
#                      V 1.2

#                   12 SEPTEMBER 2026
#   
#   Luis Gerardo Abarca     gabarca@uv.mx   luisgaa@gmail.com
#   Israel Huesca Domínguez ihuesca@uv.mx
#   
#               IIB Universidad Veracruzana
#   

#   ____________________________________________________________________________
#   Inicio                                                                  ####

#   dependencies for Food Web Topology Analysis
#CLEAR SCREEN (CONSOLE)
cat("\014")

##  ............................................................................
##  List of packages used                                                   ####

# List of ALL CRAN packages (including ig.degree.betweenness)
packages <- c(
    "igraph", 
    "stringr", 
    "cheddar", 
    "leidenAlg", 
    "plyr",
    "dplyr", 
    "ggplot2", 
    "ggpubr", 
    "sna", 
    "rnetcarto",
    "ATNr",
    "logger",
    "purrr",
    "parallel",
    "progressr",
    "ggnetwork",
    "ggthemes",
    "RColorBrewer"
)

##  ............................................................................
##  Install packages                                                        ####

# Install missing packages with error handling
install_if_missing <- function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        install.packages(pkg, dependencies = TRUE)
    }
}
library(purrr)
purrr::walk(packages, install_if_missing)  

##  ............................................................................
##  Load libraries                                                          ####

#READ ALL LIBRARIES TO BE USED
library(igraph)
library(stringr)
library(cheddar)
library(leidenAlg)
library(plyr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(sna)
library(rnetcarto)
library(ATNr)
library(logger)
library(parallel)
library(progressr)
library(ggnetwork)
library(ggthemes)
library(RColorBrewer)

#CLEAR SCREEN (CONSOLE)
cat("\014")
cat("------ ALL PACKAGES INSTALLED AND READY TO BE USED -----\n")
