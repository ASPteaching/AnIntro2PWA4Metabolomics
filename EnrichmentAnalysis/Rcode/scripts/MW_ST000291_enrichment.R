## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval = FALSE-------------------------------------------------------------
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# 
# BiocManager::install("fobitools")

## ----warning = FALSE, message = FALSE, comment = FALSE------------------------
library(fobitools)

## ----warning = FALSE, message = FALSE, comment = FALSE------------------------
# CRAN
library(tidyverse)
library(rvest)
library(ggrepel)
library(kableExtra)

# Bioconductor
library(POMA)
library(metabolomicsWorkbenchR)
library(SummarizedExperiment)

## ----warning = FALSE, message = FALSE, comment = FALSE------------------------
data_negative_mode <- do_query(
  context = "study",
  input_item = "analysis_id",
  input_value = "AN000465",
  output_item = "SummarizedExperiment")

data_positive_mode <- do_query(
  context = "study",
  input_item = "analysis_id",
  input_value = "AN000464",
  output_item = "SummarizedExperiment")

## ----metabolitenames, message = FALSE, warning = FALSE, comment = FALSE, echo = FALSE, fig.align = "center", out.width = "100%", fig.cap = 'Metabolite identifiers of the ST000291 Metabolomics Workbench study.'----
# knitr::include_graphics("figure/ST000291_names.png")

## ----warning = FALSE, message = FALSE, comment = FALSE------------------------
metaboliteNamesURL <- "https://www.metabolomicsworkbench.org/data/show_metabolites_by_study.php?STUDY_ID=ST000291&SEARCH_TYPE=KNOWN&STUDY_TYPE=MS&RESULT_TYPE=1"

metaboliteNames <- metaboliteNamesURL %>% 
  read_html() %>% 
  html_nodes(".datatable")

metaboliteNames_negative <- metaboliteNames %>%
  .[[1]] %>%
  html_table() %>%
  dplyr::select(`Metabolite Name`, PubChemCompound_ID, `Kegg Id`)

metaboliteNames_positive <- metaboliteNames %>%
  .[[2]] %>%
  html_table() %>%
  dplyr::select(`Metabolite Name`, PubChemCompound_ID, `Kegg Id`)

metaboliteNames <- bind_rows(metaboliteNames_negative, metaboliteNames_positive) %>%
  dplyr::rename(names = 1, PubChem = 2, KEGG = 3) %>%
  mutate(KEGG = ifelse(KEGG == "-", "UNKNOWN", KEGG),
         PubChem = ifelse(PubChem == "-", "UNKNOWN", PubChem)) %>%
  filter(!duplicated(PubChem))

## ----warning = FALSE, message = FALSE, comment = FALSE------------------------
## negative mode features 
features_negative <- assay(data_negative_mode) %>%
  dplyr::slice(-n())
rownames(features_negative) <- rowData(data_negative_mode)$metabolite[1:(length(rowData(data_negative_mode)$metabolite)-1)]

## positive mode features
features_positive <- assay(data_positive_mode) %>%
  dplyr::slice(-n())
rownames(features_positive) <- rowData(data_positive_mode)$metabolite[1:(length(rowData(data_positive_mode)$metabolite)-1)]

## combine positive and negative mode and set PubChem IDs as feature names
features <- bind_rows(features_negative, features_positive) %>%
  tibble::rownames_to_column("names") %>%
  right_join(metaboliteNames, by = "names") %>%
  select(-names, -KEGG) %>%
  tibble::column_to_rownames("PubChem")

## metadata
pdata <- colData(data_negative_mode) %>% # or "data_positive_mode". They are equal
  as.data.frame() %>%
  tibble::rownames_to_column("ID") %>%
  mutate(Treatment = case_when(Treatment == "Baseline urine" ~ "Baseline",
                               Treatment == "Urine after drinking apple juice" ~ "Apple",
                               Treatment == "Urine after drinking cranberry juice" ~ "Cranberry"))

## ----warning = FALSE, message = FALSE, comment = FALSE------------------------
pdata$Treatment <- as.character(pdata$Treatment)
pdata$Treatment <- factor(pdata$Treatment,
                          levels = c("Baseline", "Apple", "Cranberry"))

data_sumexp <- PomaCreateObject(metadata = pdata, features = t(features))

## ----warning = FALSE, message = FALSE, comment = FALSE------------------------

data_preprocessed <- data_sumexp %>%
  PomaImpute() %>%             # sense group_by
  PomaNorm()

# data_preprocessed <- data_sumexp %>%
#   PomaImpute(group_by = "Treatment") %>%
#   PomaNorm()

## ----warning = FALSE, message = FALSE, comment = FALSE------------------------

#######################
### Apple vs Baseline
#######################

limma_resAB <- data_preprocessed %>%
  PomaLimma(contrast = "Apple-Baseline", adjust = "fdr") %>%
  dplyr::rename(PubChemCID = feature) %>% 
  dplyr::mutate(PubChemCID = gsub("X", "", PubChemCID)) %>%
  filter(PubChemCID != "UNKNOWN")

# show the first 10 features
limma_resAB %>%
  dplyr::slice(1L:10L) %>%
  kbl(row.names = FALSE, booktabs = TRUE) %>%
  kable_styling(latex_options = c("striped"))

#######################
### Cranberry vs Baseline
#######################

limma_resCB <- data_preprocessed %>%
  PomaLimma(contrast = "Cranberry-Baseline", adjust = "fdr") %>%
  dplyr::rename(PubChemCID = feature) %>% 
  dplyr::mutate(PubChemCID = gsub("X", "", PubChemCID)) %>%
  filter(PubChemCID != "UNKNOWN")

# show the first 10 features
limma_resCB %>%
  dplyr::slice(1L:10L) %>%
  kbl(row.names = FALSE, booktabs = TRUE) %>%
  kable_styling(latex_options = c("striped"))

#######################
### Apple vs Cranberry
#######################

limma_resAC <- data_preprocessed %>%
  PomaLimma(contrast = "Apple-Cranberry", adjust = "fdr") %>%
  dplyr::rename(PubChemCID = feature) %>% 
  dplyr::mutate(PubChemCID = gsub("X", "", PubChemCID))%>%
  filter(PubChemCID != "UNKNOWN")

# show the first 10 features
limma_resAC %>%
  dplyr::slice(1L:10L) %>%
  kbl(row.names = FALSE, booktabs = TRUE) %>%
  kable_styling(latex_options = c("striped"))

limma_res = list(ApplvsBase = limma_resAB, 
                 CranvsBase = limma_resCB,
                 ApplvsCran = limma_resAC)

save(limma_res, file="MW_ST000291-limma_res.Rda")


#########################################################
### Fins aquí tenim el resultat de l'analisi amb limma
### La vinyeta segueix amb conversió a ids de FOBI 
### per fer servir fobitools
#########################################################


## ----warning = FALSE, message = FALSE, comment = FALSE------------------------
limma_FOBI_names <- limma_res %>%
  dplyr::pull("PubChemCID") %>%
  fobitools::id_convert()

# show the ID conversion results
limma_FOBI_names %>%
  head() %>%
  kbl(row.names = FALSE, booktabs = TRUE) %>%
  kable_styling(latex_options = c("striped"))

## ----warning = FALSE, message = FALSE, comment = FALSE------------------------
limma_FOBI_names <- limma_FOBI_names %>%
  dplyr::right_join(limma_res, by = "PubChemCID") %>%
  dplyr::arrange(-dplyr::desc(pvalue))

## ----warning = FALSE, message = FALSE, comment = FALSE------------------------
metaboliteList <- limma_FOBI_names$FOBI[limma_FOBI_names$pvalue < 0.01]
metaboliteUniverse <- limma_FOBI_names$FOBI

fobitools::ora(metaboliteList = metaboliteList,
               metaboliteUniverse = metaboliteUniverse,
               pvalCutoff = 0.5) %>%
  kbl(row.names = FALSE, booktabs = TRUE) %>%
  kable_styling(latex_options = c("striped"))

## ----warning = FALSE, message = FALSE, comment = FALSE------------------------
limma_FOBI_msea <- limma_FOBI_names %>%
  dplyr::select(FOBI, pvalue) %>%
  dplyr::filter(!is.na(FOBI)) %>%
  dplyr::arrange(-dplyr::desc(abs(pvalue)))

FOBI_msea <- as.vector(limma_FOBI_msea$pvalue)
names(FOBI_msea) <- limma_FOBI_msea$FOBI

msea_res <- fobitools::msea(FOBI_msea, pvalCutoff = 0.06)

msea_res %>%
  kbl(row.names = FALSE, booktabs = TRUE) %>%
  kable_styling(latex_options = c("striped"))

## ----warning = FALSE, message = FALSE, comment = FALSE------------------------
fobi_graph(terms = c("FOODON:03301123", "FOODON:03301702"), 
           get = "anc", 
           labels = TRUE, 
           labelsize = 6)

## ----warning = FALSE, message = FALSE, comment = FALSE, fig.width = 12, fig.height = 8----
ggplot(msea_res, aes(x = -log10(pval), y = NES, color = NES, size = classSize, label = className)) +
  xlab("-log10(P-value)") +
  ylab("NES (Normalized Enrichment Score)") +
  geom_point() +
  ggrepel::geom_label_repel(color = "black", size = 7) +
  theme_bw() +
  theme(legend.position = "top",
        text = element_text(size = 22)) +
  scale_color_viridis_c() +
  scale_size(guide = "none")

## ----warning = FALSE, message = FALSE, comment = FALSE, fig.width = 12, fig.height = 10----
FOBI_terms <- msea_res %>% 
  unnest(cols = leadingEdge)

fobitools::fobi %>%
  filter(FOBI %in% FOBI_terms$leadingEdge) %>%
  pull(id_code) %>%
  fobi_graph(get = "anc", 
             labels = TRUE, 
             legend = TRUE, 
             labelsize = 6, 
             legendSize = 20)

## -----------------------------------------------------------------------------
sessionInfo()
