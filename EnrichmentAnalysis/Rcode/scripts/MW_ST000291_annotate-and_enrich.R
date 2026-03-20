
#########################################################
##       PIPELINE FINAL — Metabolòmica + ORA/MSEA
##   KEGG (PubChem) + SMPDB (HMDB) + SMPDB→PubChem
#########################################################

library(tidyverse)
library(KEGGREST)
library(localEnrichment)
library(metaboliteIDmapping)

#########################################################
## 0. INPUTS NECESSARIS
#########################################################

load(file="MW_ST000291-limma_res.Rda")
limma_resAB <- limma_res$ApplvsBase
limma_resAC <- limma_res$ApplvsCran
limma_resCB <- limma_res$CranvsBase

load("SMPDB_pathways.Rda")

#########################################################
## 1. PREPARAR LIMMA (PubChem com ID únic)
#########################################################

fix_limma <- function(df) {
  df %>%
    filter(!is.na(PubChemCID), PubChemCID != "UNKNOWN") %>%
    mutate(PubChemCID = as.character(PubChemCID)) %>%
    distinct(PubChemCID, .keep_all = TRUE)
}

limma_resAB <- fix_limma(limma_resAB)
limma_resCB <- fix_limma(limma_resCB)
limma_resAC <- fix_limma(limma_resAC)

#########################################################
## 2. PREPARAR SMPDB EN FORMAT HMDB (id natiu)
#########################################################

SMPDB_path2hmdb <- SMPDBset_df %>%
  separate_rows(Metabolites, sep = ",") %>%
  mutate(HMDB = str_trim(Metabolites)) %>%
  select(pathway = set_name, HMDB)

#########################################################
## 3. SMPDB → PubChem (opcional: per integrar amb KEGG PubChem)
#########################################################

hmdb_to_pubchem <- metabolitesMapping %>%
  filter(!is.na(HMDB), !is.na(CID)) %>%
  transmute(HMDB, PubChemCID = as.character(CID)) %>%
  distinct()

SMPDB_path2cid <- SMPDB_path2hmdb %>%
  left_join(hmdb_to_pubchem, by = "HMDB") %>%
  filter(!is.na(PubChemCID)) %>%
  distinct(pathway, PubChemCID)

#########################################################
## 4. KEGG — Pathways humans → PubChem CID
#########################################################

hsa_ids <- KEGGREST::keggList("pathway", "hsa") %>%
  names() %>%
  sub("^path:", "", .)

extract_pathway_compounds <- function(pid) {
  k <- KEGGREST::keggGet(pid)[[1]]
  if (is.null(k$COMPOUND)) return(NULL)
  tibble(pathway = pid, kegg_compound = names(k$COMPOUND))
}

kegg_compounds <- map_dfr(hsa_ids, extract_pathway_compounds)

kegg_conv <- KEGGREST::keggConv("pubchem", "compound")

# Compte que el seguent pas tarda molt
kegg_cpd2cid <- tibble(
  kegg_compound = sub("^cpd:", "", names(kegg_conv)),
  PubChemCID = sub("^pubchem:", "", unname(kegg_conv))
)

KEGG_path2cid <- kegg_compounds %>%
  inner_join(kegg_cpd2cid, by = "kegg_compound") %>%
  select(pathway, PubChemCID) %>%
  distinct()

#########################################################
## 5. CONSTRUIR ENRICHMENT SETS PER localEnrichmentR
#########################################################
# Converteix un mapping "long" (pathway, ID) -> format wide per EnrichmentSet
make_ES_wide <- function(df, id_col, metadata, sep = ";") {
  df_wide <- df %>%
    group_by(pathway) %>%
    summarise(
      feature_ids = paste(unique(.data[[id_col]]), collapse = sep),
      .groups = "drop"
    ) %>%
    transmute(
      set_id = pathway,
      set_name = pathway,
      feature_ids = feature_ids
    )
  
  ES <- EnrichmentSet(
    df_wide,
    metadata = metadata,
    sep = sep
  )
  
  validateEnrichmentSet(ES)
  ES
}

## KEGG → PubChem
ES_KEGG_CID <- make_ES_wide(
  df = KEGG_path2cid,
  id_col = "PubChemCID",
  metadata = list(
    mapping_name = "KEGG_PubChem",
    feature_id_type = "PubChemCID",
    feature_species = "human",
    set_source = "KEGG",
    version = "2025.01",
    description = "Human KEGG pathways mapped to PubChem CID"
  ),
  sep = ";"
)

## SMPDB → HMDB (natiu)
# ES_SMPDB_HMDB <- make_ES_wide(
#   df = SMPDB_path2hmdb,
#   id_col = "HMDB",
#   metadata = list(
#     mapping_name = "SMPDB_HMDB",
#     feature_id_type = "HMDB",
#     feature_species = "human",
#     set_source = "SMPDB",
#     version = "1.0",
#     description = "SMPDB pathways in HMDB"
#   ),
#   sep = ";"
# )

## SMPDB → PubChem
ES_SMPDB_CID <- make_ES_wide(
  df = SMPDB_path2cid,
  id_col = "PubChemCID",
  metadata = list(
    mapping_name = "SMPDB_PubChem",
    feature_id_type = "PubChemCID",
    feature_species = "human",
    set_source = "SMPDB",
    version = "1.0",
    description = "SMPDB pathways mapped to PubChem CID"
  ),
  sep = ";"
)

# validateEnrichmentSet(ES_SMPDB_HMDB)
validateEnrichmentSet(ES_SMPDB_CID)
validateEnrichmentSet(ES_KEGG_CID)


#########################################################
## 6. PREPARAR SIG / BACKGROUND PER LIMMA
#########################################################

prep_sig_bg <- function(df, p_cut = 0.05) {
  sig <- df %>%
    filter(adj_pvalue < p_cut) %>%
    pull(PubChemCID) %>%
    unique()
  
  bg  <- df %>%
    pull(PubChemCID) %>%
    unique()
  
  list(sig = sig, bg = bg)
}

AB <- prep_sig_bg(limma_resAB, p_cut = 0.25)
CB <- prep_sig_bg(limma_resCB, p_cut = 0.25)
AC <- prep_sig_bg(limma_resAC, p_cut = 0.25)


filter_ES_by_background <- function(ES, background) {
  
  mp <- as(ES, "list")  # mapping: set_name -> vector de IDs
  
  keep_names <- names(mp)[vapply(mp, function(v)
    any(v %in% background), logical(1))]
  
  message("Sets abans del filtratge: ", length(mp))
  message("Sets després del filtratge: ", length(keep_names))
  
  df <- ES@data %>%
    dplyr::filter(set_name %in% keep_names)
  
  # SEP NO ÉS UN SLOT → cal tornar-lo a passar manualment
  EnrichmentSet(df, ES@metadata, sep = ";")
}


### Background is the same for all three lists!!!

ES_SMPDB_CID_filt <- filter_ES_by_background(ES_SMPDB_CID, AB$bg)
# ES_SMPDB_HMDB_filt <- filter_ES_by_background(ES_SMPDB_HMDB, AB$bg)
ES_KEGG_CID_filt <- filter_ES_by_background(ES_KEGG_CID, AB$bg)


#########################################################
## 7. ORA AMB localEnrichmentR
#########################################################

run_ORA <- function(sig, bg, ES) {
  ora_test(
    eset      = ES,
    selected  = sig,
    background = bg
  )
}


# Apple vs Baseline

ORA_AB_SMPDB_CID  <- run_ORA(AB$sig, AB$bg, ES_SMPDB_CID_filt)
ORA_AB_KEGG       <- run_ORA(AB$sig, AB$bg, ES_KEGG_CID_filt)

# Cranberry vs Baseline

ORA_CB_SMPDB_CID  <- run_ORA(CB$sig, CB$bg, ES_SMPDB_CID_filt)
ORA_CB_KEGG       <- run_ORA(CB$sig, CB$bg, ES_KEGG_CID_filt)

# Apple vs Cranberry

ORA_AC_SMPDB_CID  <- run_ORA(AC$sig, AC$bg, ES_SMPDB_CID_filt)
ORA_AC_KEGG       <- run_ORA(AC$sig, AC$bg, ES_KEGG_CID_filt)

#########################################################
## 8. MSEA (fgsea-style) amb localEnrichmentR
#########################################################

prep_ranks <- function(df) {
  r <- sign(df$logFC) * -log10(df$P.Value)
  names(r) <- df$PubChemCID
  r
}

rAB <- prep_ranks(limma_resAB)
rCB <- prep_ranks(limma_resCB)
rAC <- prep_ranks(limma_resAC)

run_RSEA <- function(ranks, ES) {
  localEnrichment(
    sig    = ranks,
    set    = ES,
    method = "RSEA"
  )
}


RSEA_AB_KEGG <- run_RSEA(rAB, ES_KEGG_CID)
RSEA_CB_KEGG <- run_RSEA(rCB, ES_KEGG_CID)
RSEA_AC_KEGG <- run_RSEA(rAC, ES_KEGG_CID)

RSEA_AB_SMPDB <- run_RSEA(rAB, ES_SMPDB_CID)
RSEA_CB_SMPDB <- run_RSEA(rCB, ES_SMPDB_CID)
RSEA_AC_SMPDB <- run_RSEA(rAC, ES_SMPDB_CID)

#########################################################
## 9. GUARDAR TOTS ELS RESULTATS
#########################################################

save(
  ES_SMPDB_CID, ES_KEGG_CID,
  ORA_AB_SMPDB_CID, ORA_AB_KEGG,
  ORA_CB_SMPDB_CID, ORA_CB_KEGG,
  ORA_AC_SMPDB_CID, ORA_AC_KEGG,
  RSEA_AB_KEGG, RSEA_AB_SMPDB,
  RSEA_CB_KEGG, RSEA_CB_SMPDB,
  RSEA_AC_KEGG, RSEA_AC_SMPDB, 
  file = "testResults/Pipeline_Final_ST000291_Enrichment.Rda"
)
