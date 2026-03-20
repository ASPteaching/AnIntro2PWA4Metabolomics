if (!require("devtools"))
  install.packages('devtools', repos="http://cran.rstudio.com/")
if (!require("RCurl"))
  install.packages('RCurl', repos="http://cran.rstudio.com/")
if (!require("pacman"))
  install.packages('pacman', repos="http://cran.rstudio.com/")
library(devtools)

library(RCurl)

library(pacman)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

pacman::p_load(GGally)

pacman::p_load(DT)

pacman::p_load(RCurl)

pacman::p_load(RJSONIO)

pacman::p_load(ape)

pacman::p_load(devEMF)

pacman::p_load(dynamicTreeCut)

pacman::p_load(extrafont)

pacman::p_load(ggplot2)

pacman::p_load(ggpubr)

pacman::p_load(ggrepel)

pacman::p_load(grid)

pacman::p_load(htmlwidgets)

pacman::p_load(igraph)

pacman::p_load(magrittr)

pacman::p_load(network)

pacman::p_load(officer)

pacman::p_load(openxlsx)

pacman::p_load(phytools)

pacman::p_load(plotly)

pacman::p_load(plotrix)

pacman::p_load(rcdk)

pacman::p_load(readxl)

pacman::p_load(rvg)

pacman::p_load(sna)

pacman::p_load(visNetwork)

source("https://raw.githubusercontent.com/barupal/ChemRICH/refs/heads/master/chemrich_minimum_analysis.R")


