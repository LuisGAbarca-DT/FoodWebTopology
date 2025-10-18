#  
#  
#               Food Web Topology Analysis
#                   Master Script
#   
#                   31 AUGUST 2025
#   
#   Luis Gerardo Abarca     gabarca@uv.mx   luisgaa@gmail.com
#   Israel Huesca Domínguez ihuesca@uv.mx
#   
#               IIB Universidad Veracruzana
#   

#   dependencies for Food Web Topology Analysis
#CLEAR SCREEN (CONSOLE)
cat("\014")
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
    "nortest", 
    "sna", 
    "rnetcarto",
    "ATNr",
    "logger",
    "purrr"
)

# Install missing packages with error handling
install_if_missing <- function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        install.packages(pkg, dependencies = TRUE)
    }
}
library(purrr)
purrr::walk(packages, install_if_missing)  # Installs all gracefully

#READ ALL LIBRARIES TO BE USED
library(igraph)
library(stringr)
library(cheddar)
library(leidenAlg)
library(plyr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(nortest)
library(sna)
library(rnetcarto)
library(ATNr)
library(logger)
#CLEAR SCREEN (CONSOLE)
cat("\014")
cat("------ ALL PACKAGES INSTALLED AND READY TO BE USED -----\n")
