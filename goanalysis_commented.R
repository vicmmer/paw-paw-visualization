# Trying to figure out TOPGO analysis based on Suzy's code with some modifications
## tutorial she linked:
# http://avrilomics.blogspot.com/2015/07/using-topgo-to-test-for-go-term.html

# ---- Load packages ----
setwd("C:/Users/vicme/OneDrive/Desktop/PhD/pawpaw_project/Suzy_pipeline/Goanalysis")  # set my working folder so files read/write here
library(topGO)       # GO enrichment functions live here
library(Rgraphviz)   # topGO needs this to draw the little GO network PDFs
library(readr)

# ---- Define the gene universe (mapping from gene/OG -> GO terms) ----
geneID2GO <- readMappings(file = "all_go_unique_clean.tsv")  # read the mapping (each line: gene <tab> GO1,GO2,...)
geneUniverse <- names(geneID2GO)                              # all gene/OG IDs that appear in the mapping = my “universe”

# helper to read a plain one-ID-per-line file
read_interest <- function(path) {
  x <- read.table(path, header = FALSE, stringsAsFactors = FALSE)  # no header, keep as characters
  as.character(x[[1]])                                             # return just the first column as a character vector
}

# ---- 1. Pawpaw figure (BP and MF) ----
#load the pawpaw gene/OG list 
pawpaw_ids <- read_interest("new_pawpaw_list.txt")                 # this is the “selected” set for pawpaw

# build the 0/1 selection vector over the full universe (1 = gene is in pawpaw_ids; 0 = not).
# topGO uses this to know which genes are "selected" when it builds the Fisher 2x2 tables per GO term.
pawpaw_flag <- factor(as.integer(geneUniverse %in% pawpaw_ids))    # 1 if in pawpaw_ids, 0 otherwise
names(pawpaw_flag) <- geneUniverse                                 # names must be the gene IDs

# Pawpaw  - BP (Biological Processes)
myGOdataBP <- new("topGOdata",               # build the topGO object for Biological Process
                  description = "IC unique OGs - Pawpaw BP", # this label shows up on plots
                  ontology    = "BP",                         # BP = Biological Process
                  allGenes    = pawpaw_flag,                  # my 0/1 vector saying which genes are “selected”
                  annot       = annFUN.gene2GO,               # tell topGO we’re passing a gene->GO list
                  gene2GO     = geneID2GO)                    # the actual mapping

resultFisher <- runTest(myGOdataBP, algorithm = "classic", statistic = "fisher")  # classic Fisher exact test

bp_scores <- score(resultFisher)                                    # numeric p-values for all BP terms (some can be NA)
bp_n_all  <- sum(!is.na(bp_scores))                                 # how many BP terms actually have a score

if (bp_n_all > 0) {                                                 # only make a table if there’s something to report
  pawpaw_BP <- GenTable(myGOdataBP,                                 # grab ALL terms (not just top 10)
                        classicFisher = resultFisher,
                        orderBy       = "classicFisher",                                # sort by that column
                        ranksOf       = "classicFisher",
                        topNodes      = bp_n_all)
  
  # add numeric p + FDR and clean up column names
  pawpaw_BP$P_Value <- bp_scores[pawpaw_BP$GO.ID]                   # align numeric p-values by GO.ID
  pawpaw_BP$FDR_BH  <- p.adjust(bp_scores, method = "BH")[pawpaw_BP$GO.ID]  # BH FDR across all BP terms
  
  names(pawpaw_BP)[names(pawpaw_BP) == "GO.ID"]       <- "GO_ID"                 # clearer names
  names(pawpaw_BP)[names(pawpaw_BP) == "Term"]        <- "GO_Term"
  names(pawpaw_BP)[names(pawpaw_BP) == "Annotated"]   <- "N_Annotated_Universe"
  names(pawpaw_BP)[names(pawpaw_BP) == "Significant"] <- "N_Significant_Selected"
  names(pawpaw_BP)[names(pawpaw_BP) == "Expected"]    <- "Expected_Significant"
  
  pawpaw_BP$classicFisher <- NULL                                     # drop the pretty p-value string column
  pawpaw_BP$Label    <- "pawpaw"                                      # tag which set this is
  pawpaw_BP$Category <- "BP"                                          # tag the ontology
  
  pawpaw_BP <- pawpaw_BP[, c("Label","Category","GO_ID","GO_Term",    # reorder to my preferred column order
                             "N_Annotated_Universe","N_Significant_Selected","Expected_Significant",
                             "P_Value","FDR_BH"), drop = FALSE]
  
  write.csv(pawpaw_BP, "pawpaw_BP_topgo.csv", row.names = FALSE)      # write the per-ontology CSV
} else {
  # if no terms were feasible, still write an empty file with headers so the combine step won’t break
  write.csv(data.frame(Label=character(),Category=character(),GO_ID=character(),GO_Term=character(),
                       N_Annotated_Universe=integer(),N_Significant_Selected=integer(),Expected_Significant=numeric(),
                       P_Value=numeric(),FDR_BH=numeric()),
            "pawpaw_BP_topgo.csv", row.names = FALSE)
}

# draw a small GO graph PDF (cap nodes so it doesn’t get messy)
showSigOfNodes(myGOdataBP, bp_scores, firstSigNodes = min(10, max(1, bp_n_all)), useInfo = "all")
printGraph(myGOdataBP, resultFisher, firstSigNodes = min(10, max(1, bp_n_all)),
           fn.prefix = "tGO_pawpaw_BP", useInfo = "all", pdfSW = TRUE)

# Pawpaw  - MF (molecular function)
myGOdataMF <- new("topGOdata",               # same as above but for Molecular Function
                  description = "IC unique OGs - Pawpaw MF",
                  ontology    = "MF",
                  allGenes    = pawpaw_flag,
                  annot       = annFUN.gene2GO,
                  gene2GO     = geneID2GO)


resultFisher <- runTest(myGOdataMF, algorithm = "classic", statistic = "fisher")

mf_scores <- score(resultFisher)                                    # numeric p-values for MF terms
mf_n_all  <- sum(!is.na(mf_scores))

if (mf_n_all > 0) {
  pawpaw_MF <- GenTable(myGOdataMF,
                        classicFisher = resultFisher,
                        orderBy       = "classicFisher",
                        ranksOf       = "classicFisher",
                        topNodes      = mf_n_all)
  
  pawpaw_MF$P_Value <- mf_scores[pawpaw_MF$GO.ID]
  pawpaw_MF$FDR_BH  <- p.adjust(mf_scores, method = "BH")[pawpaw_MF$GO.ID]
  
  names(pawpaw_MF)[names(pawpaw_MF) == "GO.ID"]       <- "GO_ID"
  names(pawpaw_MF)[names(pawpaw_MF) == "Term"]        <- "GO_Term"
  names(pawpaw_MF)[names(pawpaw_MF) == "Annotated"]   <- "N_Annotated_Universe"
  names(pawpaw_MF)[names(pawpaw_MF) == "Significant"] <- "N_Significant_Selected"
  names(pawpaw_MF)[names(pawpaw_MF) == "Expected"]    <- "Expected_Significant"
  
  pawpaw_MF$classicFisher <- NULL
  pawpaw_MF$Label    <- "pawpaw"
  pawpaw_MF$Category <- "MF"
  
  pawpaw_MF <- pawpaw_MF[, c("Label","Category","GO_ID","GO_Term",
                             "N_Annotated_Universe","N_Significant_Selected","Expected_Significant",
                             "P_Value","FDR_BH"), drop = FALSE]
  
  write.csv(pawpaw_MF, "pawpaw_MF_topgo.csv", row.names = FALSE)
} else {
  write.csv(data.frame(Label=character(),Category=character(),GO_ID=character(),GO_Term=character(),
                       N_Annotated_Universe=integer(),N_Significant_Selected=integer(),Expected_Significant=numeric(),
                       P_Value=numeric(),FDR_BH=numeric()),
            "pawpaw_MF_topgo.csv", row.names = FALSE)
}

showSigOfNodes(myGOdataMF, mf_scores, firstSigNodes = min(10, max(1, mf_n_all)), useInfo = "all")
printGraph(myGOdataMF, resultFisher, firstSigNodes = min(10, max(1, mf_n_all)),
           fn.prefix = "tGO_pawpaw_MF", useInfo = "all", pdfSW = TRUE)


#---- 2. Annonaceae figure (BP, CC and MF) ----

# load the annonaceae gene/OG list 
annon_ids <- read_interest("new_annonaceae_list.txt")                # the “selected” set for annonaceae

#make the 0/1 vector over the universe for annonaceae
annon_flag <- factor(as.integer(geneUniverse %in% annon_ids))
names(annon_flag) <- geneUniverse

# Annonaceae and BP 
myGOdataBP <- new("topGOdata",
                  description = "IC unique OGs - Annonaceae BP",
                  ontology    = "BP",
                  allGenes    = annon_flag,
                  annot       = annFUN.gene2GO,
                  gene2GO     = geneID2GO)

resultFisher <- runTest(myGOdataBP, algorithm = "classic", statistic = "fisher")

bp_scores <- score(resultFisher)
bp_n_all  <- sum(!is.na(bp_scores))

if (bp_n_all > 0) {
  annon_BP <- GenTable(myGOdataBP,
                       classicFisher = resultFisher,
                       orderBy       = "classicFisher",
                       ranksOf       = "classicFisher",
                       topNodes      = bp_n_all)
  
  annon_BP$P_Value <- bp_scores[annon_BP$GO.ID]
  annon_BP$FDR_BH  <- p.adjust(bp_scores, method = "BH")[annon_BP$GO.ID]
  
  names(annon_BP)[names(annon_BP) == "GO.ID"]       <- "GO_ID"
  names(annon_BP)[names(annon_BP) == "Term"]        <- "GO_Term"
  names(annon_BP)[names(annon_BP) == "Annotated"]   <- "N_Annotated_Universe"
  names(annon_BP)[names(annon_BP) == "Significant"] <- "N_Significant_Selected"
  names(annon_BP)[names(annon_BP) == "Expected"]    <- "Expected_Significant"
  
  annon_BP$classicFisher <- NULL
  annon_BP$Label    <- "annonaceae"
  annon_BP$Category <- "BP"
  
  annon_BP <- annon_BP[, c("Label","Category","GO_ID","GO_Term",
                           "N_Annotated_Universe","N_Significant_Selected","Expected_Significant",
                           "P_Value","FDR_BH"), drop = FALSE]
  
  write.csv(annon_BP, "annonaceae_BP_topgo.csv", row.names = FALSE)
} else {
  write.csv(data.frame(Label=character(),Category=character(),GO_ID=character(),GO_Term=character(),
                       N_Annotated_Universe=integer(),N_Significant_Selected=integer(),Expected_Significant=numeric(),
                       P_Value=numeric(),FDR_BH=numeric()),
            "annonaceae_BP_topgo.csv", row.names = FALSE)
}

showSigOfNodes(myGOdataBP, bp_scores, firstSigNodes = min(10, max(1, bp_n_all)), useInfo = "all")
printGraph(myGOdataBP, resultFisher, firstSigNodes = min(10, max(1, bp_n_all)),
           fn.prefix = "tGO_annonaceae_BP", useInfo = "all", pdfSW = TRUE)

# Annonaceae and MF 
myGOdataMF <- new("topGOdata",
                  description = "IC unique OGs - Annonaceae MF",
                  ontology    = "MF",
                  allGenes    = annon_flag,
                  annot       = annFUN.gene2GO,
                  gene2GO     = geneID2GO)

resultFisher <- runTest(myGOdataMF, algorithm = "classic", statistic = "fisher")

mf_scores <- score(resultFisher)
mf_n_all  <- sum(!is.na(mf_scores))

if (mf_n_all > 0) {
  annon_MF <- GenTable(myGOdataMF,
                       classicFisher = resultFisher,
                       orderBy       = "classicFisher",
                       ranksOf       = "classicFisher",
                       topNodes      = mf_n_all)
  
  annon_MF$P_Value <- mf_scores[annon_MF$GO.ID]
  annon_MF$FDR_BH  <- p.adjust(mf_scores, method = "BH")[annon_MF$GO.ID]
  
  names(annon_MF)[names(annon_MF) == "GO.ID"]       <- "GO_ID"
  names(annon_MF)[names(annon_MF) == "Term"]        <- "GO_Term"
  names(annon_MF)[names(annon_MF) == "Annotated"]   <- "N_Annotated_Universe"
  names(annon_MF)[names(annon_MF) == "Significant"] <- "N_Significant_Selected"
  names(annon_MF)[names(annon_MF) == "Expected"]    <- "Expected_Significant"
  
  annon_MF$classicFisher <- NULL
  annon_MF$Label    <- "annonaceae"
  annon_MF$Category <- "MF"
  
  annon_MF <- annon_MF[, c("Label","Category","GO_ID","GO_Term",
                           "N_Annotated_Universe","N_Significant_Selected","Expected_Significant",
                           "P_Value","FDR_BH"), drop = FALSE]
  
  write.csv(annon_MF, "annonaceae_MF_topgo.csv", row.names = FALSE)
} else {
  write.csv(data.frame(Label=character(),Category=character(),GO_ID=character(),GO_Term=character(),
                       N_Annotated_Universe=integer(),N_Significant_Selected=integer(),Expected_Significant=numeric(),
                       P_Value=numeric(),FDR_BH=numeric()),
            "annonaceae_MF_topgo.csv", row.names = FALSE)
}

showSigOfNodes(myGOdataMF, mf_scores, firstSigNodes = min(10, max(1, mf_n_all)), useInfo = "all")
printGraph(myGOdataMF, resultFisher, firstSigNodes = min(10, max(1, mf_n_all)),
           fn.prefix = "tGO_annonaceae_MF", useInfo = "all", pdfSW = TRUE)

# Annonaceae and CC 
myGOdataCC <- new("topGOdata",
                  description = "IC unique OGs - Annonaceae CC",
                  ontology    = "CC",
                  allGenes    = annon_flag,
                  annot       = annFUN.gene2GO,
                  gene2GO     = geneID2GO)

resultFisher <- runTest(myGOdataCC, algorithm = "classic", statistic = "fisher")

cc_scores <- score(resultFisher)
cc_n_all  <- sum(!is.na(cc_scores))

if (cc_n_all > 0) {
  annon_CC <- GenTable(myGOdataCC,
                       classicFisher = resultFisher,
                       orderBy       = "classicFisher",
                       ranksOf       = "classicFisher",
                       topNodes      = cc_n_all)
  
  annon_CC$P_Value <- cc_scores[annon_CC$GO.ID]
  annon_CC$FDR_BH  <- p.adjust(cc_scores, method = "BH")[annon_CC$GO.ID]
  
  names(annon_CC)[names(annon_CC) == "GO.ID"]       <- "GO_ID"
  names(annon_CC)[names(annon_CC) == "Term"]        <- "GO_Term"
  names(annon_CC)[names(annon_CC) == "Annotated"]   <- "N_Annotated_Universe"
  names(annon_CC)[names(annon_CC) == "Significant"] <- "N_Significant_Selected"
  names(annon_CC)[names(annon_CC) == "Expected"]    <- "Expected_Significant"
  
  annon_CC$classicFisher <- NULL
  annon_CC$Label    <- "annonaceae"
  annon_CC$Category <- "CC"
  
  annon_CC <- annon_CC[, c("Label","Category","GO_ID","GO_Term",
                           "N_Annotated_Universe","N_Significant_Selected","Expected_Significant",
                           "P_Value","FDR_BH"), drop = FALSE]
  
  write.csv(annon_CC, "annonaceae_CC_topgo.csv", row.names = FALSE)
} else {
  write.csv(data.frame(Label=character(),Category=character(),GO_ID=character(),GO_Term=character(),
                       N_Annotated_Universe=integer(),N_Significant_Selected=integer(),Expected_Significant=numeric(),
                       P_Value=numeric(),FDR_BH=numeric()),
            "annonaceae_CC_topgo.csv", row.names = FALSE)
}

showSigOfNodes(myGOdataCC, cc_scores, firstSigNodes = min(10, max(1, cc_n_all)), useInfo = "all")
printGraph(myGOdataCC, resultFisher, firstSigNodes = min(10, max(1, cc_n_all)),
           fn.prefix = "tGO_annonaceae_CC", useInfo = "all", pdfSW = TRUE)


# ---- 3. Combine per-ontology CSVS into one tsv per label  ----


combine_topgo <- function(label) {
  # find files like "pawpaw_BP_topgo.csv", "pawpaw_MF_topgo.csv", etc.
  files <- list.files(pattern = paste0("^", label, "_(BP|MF|CC)_topgo\\.csv$"),
                      ignore.case = TRUE)
  if (!length(files)) { message("No files for ", label); return(invisible(NULL)) }
  
  # read each CSV and drop empties (0-row)
  dfs <- lapply(files, function(f) {
    d <- read.csv(f, check.names = FALSE)
    if (nrow(d) == 0) return(NULL)
    d
  })
  keep <- !vapply(dfs, is.null, logical(1))
  if (!any(keep)) { message("All files empty for ", label); return(invisible(NULL)) }
  dfs <- dfs[keep]
  
  # align columns (union of names) then stack them
  all_cols <- Reduce(union, lapply(dfs, names))
  dfs <- lapply(dfs, function(d) {
    miss <- setdiff(all_cols, names(d))
    if (length(miss)) for (m in miss) d[[m]] <- NA
    d[all_cols]
  })
  out <- do.call(rbind, dfs)
  
  # put the columns in the order I like
  front <- c("Label","Category","GO_ID","GO_Term",
             "N_Annotated_Universe","N_Significant_Selected","Expected_Significant",
             "P_Value","FDR_BH")
  front <- front[front %in% names(out)]
  out <- out[, c(front, setdiff(names(out), front)), drop = FALSE]
  
  # write the combined TSV
  write.table(out, paste0(label, "_ALL_topgo.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  message("Wrote ", label, "_ALL_topgo.tsv with ", nrow(out), " rows (", length(dfs), " files)")
}


# run the combiner for each label to get the final TSVs
combine_topgo("pawpaw")
combine_topgo("annonaceae")

#  ---- 4. Load the combined TSVs back in and look at them  ----

annonaceae <- read.delim("annonaceae_ALL_topgo.tsv", sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
View(annonaceae)   # reload to see if it looks nice and pretty :P 

pawpaw <- read.delim("pawpaw_ALL_topgo.tsv", sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
View(pawpaw)       # same thing for pawpaw

