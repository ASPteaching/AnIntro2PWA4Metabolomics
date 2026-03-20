#' Prepare a metabolite list for pathway or enrichment analysis
#'
#' @description
#' This function filters and extracts a metabolite list from a single dimension
#' (i.e., one element of a list such as `metab_list[[i]]`) based on user-defined
#' thresholds for adjusted p-values and correlation values, or alternatively
#' by selecting the top percentage of metabolites according to the *p*-value.
#' It is designed to generate lists for ORA, MSEA, or ChemRich analyses,
#' depending on the selected fields.
#'
#' @param metab_table A data frame containing metabolite information for one dimension.
#'        It must include at least the columns `p.value` and `p.adj`,
#'        and optionally others such as `correlation`, `HMDB_ID`, `Variable`,
#'        `Common_name`, or `chemical_families`.
#' @param fdr_cutoff Numeric. Maximum adjusted p-value allowed for selection.
#'        Default is 1 (no filtering by significance).
#' @param cor_cutoff Numeric. Minimum absolute correlation required.
#'        Default is 0 (no filtering by correlation). Ignored if correlation not present.
#' @param cutoffPercentile Numeric between 0 and 100.
#'        If not NULL, this overrides `fdr_cutoff` and `cor_cutoff`
#'        and selects the top X percent of metabolites ranked by increasing `p.value`.
#' @param fields Character vector indicating which columns to return.
#'        Default is `c("HMDB_ID")`.
#' @param save_file Logical. If `TRUE`, saves the resulting table as a TSV file.
#'        Default is `FALSE`.
#' @param prefix Character. Prefix for the output file name (e.g., "Dim1_ORA").
#'        If `NULL`, the unique value from `metab_table$Dimension` will be used.
#' @param outdir Character. Output directory where the file will be saved.
#'        Default is `"results/metabolite_lists"`.
#' @param verbose Logical. If `TRUE`, displays informative messages about the selection.
#'        Default is `TRUE`.
#'
#' @return A tibble containing the filtered metabolites with the selected fields.
#'
#' @examples
#' \dontrun{
#' # Percentile-based selection (top 20% by p-value, with message)
#' top20 <- prepare_metabolite_list(
#'   metab_table = metab_dim1,
#'   cutoffPercentile = 20,
#'   fields = c("HMDB_ID", "p.value")
#' )
#' }
#'
#' @export
prepare_metabolite_list <- function(metab_table,
                                    fdr_cutoff = 1,
                                    cor_cutoff = 0,
                                    cutoffPercentile = NULL,
                                    fields = c("HMDB_ID"),
                                    save_file = FALSE,
                                    prefix = NULL,
                                    outdir = "results/metabolite_lists",
                                    verbose = TRUE) {
  # -------------------------------------------------------------
  # Step 1: Verify that requested fields exist in the table
  # -------------------------------------------------------------
  missing_fields <- setdiff(fields, names(metab_table))
  if (length(missing_fields) > 0) {
    stop(paste("Missing columns in the input table:",
               paste(missing_fields, collapse = ", ")))
  }
  
  # Ensure required columns exist
  required_cols <- c("p.value", "p.adj")
  missing_req <- setdiff(required_cols, names(metab_table))
  if (length(missing_req) > 0) {
    stop(paste("Missing required columns:", paste(missing_req, collapse = ", ")))
  }
  
  # -------------------------------------------------------------
  # Step 2: Filter or rank metabolites
  # -------------------------------------------------------------
  n_total <- nrow(metab_table)
  
  if (!is.null(cutoffPercentile)) {
    # ---- Percentile-based selection ----
    if (cutoffPercentile <= 0 || cutoffPercentile > 100)
      stop("cutoffPercentile must be between 0 and 100.")
    
    # Expected number of hits using the standard thresholds
    if ("correlation" %in% names(metab_table)) {
      n_expected <- metab_table %>%
        dplyr::filter(p.adj <= fdr_cutoff,
                      abs(correlation) >= cor_cutoff) %>%
        nrow()
    } else {
      n_expected <- metab_table %>%
        dplyr::filter(p.adj <= fdr_cutoff) %>%
        nrow()
    }
    
    # Actual percentile-based selection
    n_select <- max(1, floor(n_total * (cutoffPercentile / 100)))
    filtered <- metab_table %>%
      dplyr::arrange(p.value) %>%
      dplyr::slice_head(n = n_select)
    
    if (verbose) {
      perc_expected <- 100 * n_expected / n_total
      message(sprintf(
        "Selected top %.1f%% (%d of %d metabolites) by p-value.\n",
        cutoffPercentile, n_select, n_total
      ))
      message(sprintf(
        "Using standard criteria (FDR ≤ %.3f, |cor| ≥ %.3f) would select %d metabolites (%.1f%% of total).",
        fdr_cutoff, cor_cutoff, n_expected, perc_expected
      ))
    }
    
  } else {
    # ---- Standard FDR and correlation-based filtering ----
    if ("correlation" %in% names(metab_table)) {
      filtered <- metab_table %>%
        dplyr::filter(p.adj <= fdr_cutoff,
                      abs(correlation) >= cor_cutoff) %>%
        dplyr::arrange(p.adj, desc(abs(correlation)))
    } else {
      filtered <- metab_table %>%
        dplyr::filter(p.adj <= fdr_cutoff) %>%
        dplyr::arrange(p.adj)
    }
    
    if (verbose) {
      n_selected <- nrow(filtered)
      message(sprintf(
        "Selected %d metabolites (%.1f%% of %d total) using FDR ≤ %.3f and |cor| ≥ %.3f.",
        n_selected, 100 * n_selected / n_total, n_total, fdr_cutoff, cor_cutoff
      ))
    }
  }
  
  # -------------------------------------------------------------
  # Step 3: Select only the requested fields
  # -------------------------------------------------------------
  result <- filtered %>%
    dplyr::select(all_of(fields))
  
  # -------------------------------------------------------------
  # Step 4: Optionally save the filtered list as a TSV file
  # -------------------------------------------------------------
  if (save_file) {
    if (is.null(prefix)) prefix <- unique(metab_table$Dimension)
    if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
    outname <- file.path(outdir, paste0(prefix, "_filtered.tsv"))
    readr::write_tsv(result, outname, col_names = FALSE)
    if (verbose) message("Saved file: ", outname)
  }
  
  # -------------------------------------------------------------
  # Step 5: Return the filtered tibble
  # -------------------------------------------------------------
  return(result)
}


#' @title Prepare metabolite data for enrichment analysis
#' @description
#' Clean, normalize and annotate metabolite correlation tables obtained for each MFA dimension.
#' The function performs HMDB ID normalization, joins with mapping files to add common names and
#' chemical classes, and reports potential inconsistencies (duplicates, missing data).
#'
#' @param res_df A data frame containing results for one MFA dimension, with columns:
#'   Variable, correlation, p.value, Analyte_Name, HMDB_ID, and p.adj.
#' @param mappings Data frame with columns Query, Match, PubChem, ChEBI (from Food4Brain-name_map.csv).
#' @param megaINFO Data frame with columns r_name, chemical_families, method (from MEGA info file).
#' @param dim_name Character name for the dimension (used in saved filenames and messages).
#' @param save_files Logical; if TRUE, intermediate cleaned and duplicate files will be saved as Excel/Rda.
#'
#' @return A cleaned and annotated data frame ready for enrichment analysis.
#' @examples
#' res_final <- prepare_enrichment_data(resList[[1]], mappings, megaINFO, dim_name="Dim_1")
#'
prepare_enrichment_data <- function(res_df, mappings, megaINFO, dim_name="Dim", save_files=TRUE) {
  
  library(dplyr)
  library(stringr)
  library(openxlsx)
  
  message("Processing ", dim_name, "...")
  
  # 1. Simplify and rename
  res_simpl <- res_df %>%
    select(Variable, correlation, p.value, Analyte_Name, HMDB_ID, p.adj) %>%
    rename(r_name = Variable)
  
  # 2. Normalize HMDB IDs
  res_simpl <- res_simpl %>%
    mutate(HMDB_ID = str_trim(toupper(HMDB_ID)))
  
  valid_hmdb <- function(x) str_detect(x, "^HMDB[0-9]{5,7}$")
  
  removed <- res_simpl %>%
    filter(is.na(HMDB_ID) | HMDB_ID == "")
  
  res_simpl <- res_simpl %>%
    filter(!is.na(HMDB_ID) & HMDB_ID != "") %>%
    mutate(
      HMDB_ID_clean = str_replace_all(HMDB_ID, "[,;]", " "),
      HMDB_ID = str_extract(HMDB_ID_clean, "^\\S+"),
      HMDB_ID_other = str_trim(str_remove(HMDB_ID_clean, "^\\S+"))
    ) %>%
    select(-HMDB_ID_clean)
  
  if (save_files) {
    save(res_simpl, file = paste0("res_simpl_", dim_name, "_clean.Rda"))
    save(removed, file = paste0("res_removed_", dim_name, ".Rda"))
  }
  
  # 3. Check duplicated HMDB_IDs in mappings
  mappings_selected <- mappings %>%
    select(Query, Match, PubChem, ChEBI) %>%
    rename(Common_name = Match)
  
  dups <- mappings_selected %>%
    count(Query, sort = TRUE) %>%
    filter(n > 1)
  
  if (nrow(dups) > 0) {
    message("Found ", nrow(dups), " duplicated HMDB_IDs in mappings.")
    dup_details <- mappings_selected %>%
      semi_join(dups, by = "Query") %>%
      arrange(Query)
    if (save_files)
      write.xlsx(dup_details, paste0("duplicated_HMDB_mappings_", dim_name, ".xlsx"), overwrite = TRUE)
  }
  
  # 4. Merge with mapping
  n_before <- nrow(res_simpl)
  res_merged <- res_simpl %>%
    left_join(mappings_selected, by = c("HMDB_ID" = "Query")) %>%
    select(r_name, Common_name, Analyte_Name, HMDB_ID,
           correlation, p.value, p.adj, PubChem, ChEBI)
  n_after <- nrow(res_merged)
  
  if (n_after > n_before) {
    join_dups <- res_merged %>%
      count(HMDB_ID, sort = TRUE) %>%
      filter(n > 1)
    if (nrow(join_dups) > 0 && save_files) {
      join_dup_details <- res_merged %>%
        semi_join(join_dups, by = "HMDB_ID") %>%
        arrange(HMDB_ID)
      write.xlsx(join_dup_details, paste0("duplicated_after_join_", dim_name, ".xlsx"), overwrite = TRUE)
    }
  }
  
  # 5. Join MEGA chemical info
  mega_selected <- megaINFO %>%
    select(r_name, chemical_families, method)
  
  res_final <- res_merged %>%
    left_join(mega_selected, by = "r_name") %>%
    select(r_name, Common_name, Analyte_Name, HMDB_ID,
           correlation, p.value, p.adj, PubChem, ChEBI,
           chemical_families, method)
  
  # 6. Reporting
  message("Final table for ", dim_name, ": ", nrow(res_final), " metabolites")
  message("With chemical class info: ", sum(!is.na(res_final$chemical_families)))
  
  return(res_final)
}

