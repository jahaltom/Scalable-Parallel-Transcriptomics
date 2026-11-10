############################################################
# COMPLETE DESEQ2 WORKFLOW
# Input:
#   Gene_Level_Raw_Counts.txt
#   metadata.tsv
############################################################


############################################################
# 0. PACKAGES
############################################################

# Install once if needed:
#
# install.packages("BiocManager")

# BiocManager::install(c(
#     "DESeq2",
#     "EnhancedVolcano",
#     "clusterProfiler",
#     "org.Hs.eg.db",
#     "org.Mm.eg.db",
#     "enrichplot",
#     "biomaRt"
# ))

# install.packages(c(
#     "ggplot2",
#     "dplyr",
#     "tibble",
#     "pheatmap",
#     "enrichR",
#     "readxl",
#      "msigdbr"
# ))

library(yaml)

config <- read_yaml("config.yaml")

ORGDB_PACKAGE <- config$OrgDb
GENE_ID_TYPE <- config$gene_id_type
ENRICHR_DATABASES <- unlist(config$enrichr_databases)
COMPARISONS <- config$comparisons

# Optional covariates from config.yaml
COVARIATES <- unlist(config$covariates, use.names = FALSE)

if (is.null(COVARIATES)) {
    COVARIATES <- character(0)
}

COVARIATES <- unique(as.character(COVARIATES))

if ("Treatment" %in% COVARIATES) {
    stop("Do not include Treatment in covariates; it is added automatically.")
}

############################################################
# LOAD ORGANISM DATABASE
############################################################

if(!requireNamespace(
    ORGDB_PACKAGE,
    quietly = TRUE
)){
    stop(
        paste(
            "Required annotation package is not installed:",
            ORGDB_PACKAGE
        )
    )
}

library(
    ORGDB_PACKAGE,
    character.only = TRUE
)

OrgDb <- get(
    ORGDB_PACKAGE
)


library(httr)
library(msigdbr)
library(biomaRt)
library(readxl)
library(DESeq2)
library(ggplot2)
library(dplyr)
library(tibble)
library(pheatmap)
library(EnhancedVolcano)
library(clusterProfiler)
library(enrichplot)
library(enrichR)
library(limma)

############################################################
# 1. CREATE OUTPUT DIRECTORIES
############################################################

dir.create("DESeq2_results", showWarnings = FALSE)
dir.create("DESeq2_results/QC", showWarnings = FALSE)
dir.create("DESeq2_results/PCA", showWarnings = FALSE)
dir.create("DESeq2_results/DGE", showWarnings = FALSE)
dir.create("DESeq2_results/Volcano", showWarnings = FALSE)
dir.create("DESeq2_results/Heatmaps", showWarnings = FALSE)
dir.create(
    "DESeq2_results/Heatmaps/Major_Pathways",
    showWarnings = FALSE,
    recursive = TRUE
)
dir.create("DESeq2_results/GSEA", showWarnings = FALSE)


############################################################
# 2. READ RAW COUNTS
############################################################

counts <- read.delim(
    "Gene_Level_Raw_Counts.txt",
    header = TRUE,
    check.names = FALSE
)

cat("Raw count matrix dimensions:\n")
print(dim(counts))

cat("\nFirst few rows:\n")
print(head(counts))


############################################################
# 3. READ METADATA
############################################################

metadata <- read.delim(
    "metadata.tsv",
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE
)

############################################################
# STANDARDIZE AGE COLUMN
############################################################

if ("age (months)" %in% colnames(metadata)) {
    colnames(metadata)[
        colnames(metadata) == "age (months)"
    ] <- "AgeMonths"
}

if (!"AgeMonths" %in% colnames(metadata)) {
    stop("AgeMonths column is missing from metadata.")
}

metadata$AgeMonths <- as.numeric(metadata$AgeMonths)

if (anyNA(metadata$AgeMonths)) {
    stop("AgeMonths contains missing or non-numeric values.")
}

cat("\nMetadata:\n")
print(metadata)

cat("\nSamples per group:\n")
print(table(metadata$Treatment))





############################################################
# 6. CHECK THAT ALL METADATA sample EXIST
############################################################

missing_samples <- setdiff(
    metadata$sample,
    colnames(counts)
)

if(length(missing_samples) > 0){

    stop(
        paste(
            "These metadata sample are missing from the count matrix:",
            paste(missing_samples, collapse = ", ")
        )
    )
}


############################################################
# 7. KEEP ONLY SAMPLES PRESENT IN METADATA
############################################################

counts <- counts[
    ,
    c(
        "Gene",
        metadata$sample
    )
]


############################################################
# 8. PUT SAMPLES IN METADATA ORDER
############################################################

counts <- counts[
    ,
    c(
        "Gene",
        metadata$sample
    )
]


############################################################
# 9. CONVERT GENE COLUMN TO ROW NAMES
############################################################

rownames(counts) <- counts$Gene

counts$Gene <- NULL

counts <- as.matrix(counts)

storage.mode(counts) <- "numeric"


############################################################
# 10. VERIFY RAW COUNTS
############################################################

if(any(is.na(counts))){
    stop("NA values found in count matrix.")
}

if(any(counts < 0)){
    stop("Negative values found in count matrix.")
}

if(any(counts %% 1 != 0)){
    warning(
        "Non-integer values detected. Raw counts should normally be integers."
    )
}


############################################################
# 11. SET UP METADATA FOR DESEQ2
############################################################

rownames(metadata) <- metadata$sample

metadata <- metadata[
    colnames(counts),
    ,
    drop = FALSE
]


############################################################
# PREPARE DESEQ2 COVARIATES
############################################################


############################################################
# STANDARDIZE MODEL VARIABLES
############################################################

metadata$Treatment <- factor(metadata$Treatment)

# Prevent uppercase/lowercase sex labels becoming
# different factor levels.
metadata$Sex <- factor(
    toupper(trimws(as.character(metadata$Sex)))
)

# Batch is categorical, even if its labels are numeric.
metadata$Batch <- factor(
    trimws(as.character(metadata$Batch))
)

cat("\nSex distribution:\n")
print(table(metadata$Sex, useNA = "ifany"))

cat("\nBatch distribution:\n")
print(table(metadata$Batch, useNA = "ifany"))


# Verify all requested covariates exist
missing_covariates <- setdiff(
    COVARIATES,
    colnames(metadata)
)

if (length(missing_covariates) > 0) {
    stop(
        "Covariates missing from metadata.tsv: ",
        paste(missing_covariates, collapse = ", ")
    )
}

# Verify no missing values in model variables
model_variables <- c("Treatment", COVARIATES)

for (variable in model_variables) {

    if (any(is.na(metadata[[variable]]))) {
        stop(
            "Missing values found in model variable: ",
            variable
        )
    }
}

# Convert categorical covariates to factors.
# Leave numerical covariates numeric.
for (variable in COVARIATES) {

    if (is.character(metadata[[variable]]) ||
        is.logical(metadata[[variable]])) {

        metadata[[variable]] <- factor(metadata[[variable]])
    }

    if (is.factor(metadata[[variable]]) &&
        nlevels(droplevels(metadata[[variable]])) < 2) {

        stop(
            "Categorical covariate has fewer than 2 levels: ",
            variable
        )
    }

    if (is.numeric(metadata[[variable]]) &&
        length(unique(metadata[[variable]])) < 2) {

        stop(
            "Numeric covariate has no variation: ",
            variable
        )
    }
}

############################################################
# BUILD DESEQ2 DESIGN FORMULA
############################################################

design_variables <- c(COVARIATES, "Treatment")

DESEQ2_DESIGN <- stats::reformulate(
    design_variables
)

cat("\nDESeq2 design formula:\n")
print(DESEQ2_DESIGN)



############################################################
# 12. VERIFY SAMPLE ORDER
############################################################

if(!all(colnames(counts) == rownames(metadata))){
    stop("Count matrix and metadata sample order do not match.")
}

cat("\nSample matching successful.\n")


############################################################
# 13. SAVE SAMPLE INFORMATION
############################################################

write.csv(
    metadata,
    "DESeq2_results/Sample_Metadata.csv",
    row.names = FALSE
)


############################################################
# 14. LIBRARY SIZE QC
############################################################

library_sizes <- colSums(counts)

library_df <- data.frame(
    sample = names(library_sizes),
    Reads = as.numeric(library_sizes)
)

library_df <- left_join(
    library_df,
    metadata,
    by = "sample"
)

write.csv(
    library_df,
    "DESeq2_results/QC/Library_Sizes.csv",
    row.names = FALSE
)


p <- ggplot(
    library_df,
    aes(
        x = reorder(sample, Reads),
        y = Reads,
        fill = Treatment
    )
) +
    geom_col() +
    coord_flip() +
    theme_bw() +
    labs(
        title = "Raw Library Size",
        x = "Sample",
        y = "Total Raw Counts"
    ) +
    theme(
        text = element_text(size = 9),
        axis.text.y = element_text(size = 6),
        axis.text.x = element_text(size = 8),
        axis.title = element_text(size = 10),
        plot.title = element_text(size = 12),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9)
    )

ggsave(
    "DESeq2_results/QC/Library_Size.png",
    p,
    width = 8,
    height = 7,
    dpi = 300
)


############################################################
# 15. FILTER LOW-COUNT GENES
############################################################

# Keep a gene if it has at least 10 counts
# in at least min sample group size.
#
# Smallest group here appears to contain 5 samples,
# so 5 is a reasonable starting threshold.

min_group_n <- min(table(metadata$Treatment))

keep <- rowSums(counts >= 10) >= min_group_n

cat("\nFiltering summary:\n")
print(table(keep))

cat(
    "\nGenes before filtering:",
    nrow(counts),
    "\n"
)

counts_filtered <- counts[
    keep,
    ,
    drop = FALSE
]

cat(
    "Genes after filtering:",
    nrow(counts_filtered),
    "\n"
)


write.table(
    counts_filtered,
    "DESeq2_results/Filtered_Raw_Counts.txt",
    sep = "\t",
    quote = FALSE,
    col.names = NA
)


############################################################
# 16. CREATE DESEQ2 OBJECT
############################################################


############################################################
# CREATE DESEQ2 OBJECT WITH CONFIGURABLE COVARIATES
############################################################


############################################################
# CHECK BATCH CONFOUNDING AND CREATE DESEQ2 OBJECT
############################################################

# Treatment-by-batch distribution
batch_table <- table(
    Treatment = metadata$Treatment,
    Batch = metadata$Batch
)

cat("\nTreatment x Batch distribution:\n")
print(batch_table)

write.csv(
    as.data.frame.matrix(batch_table),
    "DESeq2_results/QC/Treatment_Batch_Table.csv"
)

# Construct model matrix before creating DESeq2 object
design_matrix <- model.matrix(
    DESEQ2_DESIGN,
    data = metadata
)

model_rank <- qr(design_matrix)$rank
model_columns <- ncol(design_matrix)

cat("\nDESeq2 model rank:", model_rank, "\n")
cat("Number of coefficients:", model_columns, "\n")

if (model_rank < model_columns) {
    stop(
        "Design matrix is not full rank. ",
        "Batch, Treatment, or other covariates are ",
        "perfectly confounded. Review the batch table."
    )
}

cat("\nDESeq2 model terms:\n")
print(colnames(design_matrix))

cat("\nDesign matrix condition number:\n")
print(kappa(design_matrix))

# Save design matrix for troubleshooting
write.csv(
    design_matrix,
    "DESeq2_results/QC/Design_Matrix.csv"
)

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(
    countData = round(counts_filtered),
    colData = metadata,
    design = DESEQ2_DESIGN
)





############################################################
# 17. RUN DESEQ2
############################################################

dds <- DESeq(dds)


############################################################
# 18. SAVE SIZE FACTORS
############################################################

size_factor_df <- data.frame(
    sample = names(sizeFactors(dds)),
    SizeFactor = sizeFactors(dds)
)

size_factor_df <- left_join(
    size_factor_df,
    metadata,
    by = "sample"
)

write.csv(
    size_factor_df,
    "DESeq2_results/QC/DESeq2_Size_Factors.csv",
    row.names = FALSE
)


############################################################
# 19. SAVE DESEQ2 NORMALIZED COUNTS
############################################################

normalized_counts <- counts(
    dds,
    normalized = TRUE
)

write.table(
    normalized_counts,
    "DESeq2_results/DESeq2_Normalized_Counts.txt",
    sep = "\t",
    quote = FALSE,
    col.names = NA
)


############################################################
# 20. VST TRANSFORMATION
############################################################

vsd <- vst(
    dds,
    blind = TRUE
)

vst_matrix <- assay(vsd)

write.table(
    vst_matrix,
    "DESeq2_results/VST_Counts.txt",
    sep = "\t",
    quote = FALSE,
    col.names = NA
)



############################################################
# BATCH-ADJUSTED VST FOR VISUALIZATION ONLY
############################################################

library(limma)

vst_matrix_batch_adjusted <- NULL

if ("Batch" %in% COVARIATES) {

    # Preserve the modeled biological Treatment effects.
    # Age and Sex are retained as covariates in the model.
    visual_design <- model.matrix(
        ~ AgeMonths + Sex + Treatment,
        data = metadata
    )

    visual_full_design <- model.matrix(
        ~ AgeMonths + Sex + Batch + Treatment,
        data = metadata
    )

    if (qr(visual_full_design)$rank == ncol(visual_full_design)) {

        vst_matrix_batch_adjusted <- limma::removeBatchEffect(
            vst_matrix,
            batch = metadata$Batch,
            design = visual_design
        )

        write.table(
            vst_matrix_batch_adjusted,
            "DESeq2_results/VST_Batch_Adjusted_Visualization.txt",
            sep = "\t",
            quote = FALSE,
            col.names = NA
        )

    } else {
        warning(
            "Skipping batch-adjusted visualization: ",
            "batch and biological variables are confounded."
        )
    }
}

############################################################
# SELECT MATRIX FOR HEATMAP VISUALIZATION
############################################################

heatmap_matrix <- if (!is.null(vst_matrix_batch_adjusted)) {
    vst_matrix_batch_adjusted
} else {
    vst_matrix
}

cat("\nHeatmap matrix selected: ")

if (!is.null(vst_matrix_batch_adjusted)) {
    cat("Batch-adjusted VST\n")
} else {
    cat("Original VST (not batch-adjusted)\n")
}

############################################################
# HALLMARK PATHWAY HEATMAPS
############################################################

dir.create(
    "DESeq2_results/Heatmaps/Hallmark",
    recursive = TRUE,
    showWarnings = FALSE
)

HALLMARK_PATHWAYS <- c(
    "HALLMARK_WNT_BETA_CATENIN_SIGNALING",
    "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
    "HALLMARK_HEDGEHOG_SIGNALING",
    "HALLMARK_TGF_BETA_SIGNALING",
    "HALLMARK_NOTCH_SIGNALING",
    "HALLMARK_HYPOXIA"
)

# Retrieve Hallmark gene sets for the analysis species
species_name <- if (
    ORGDB_PACKAGE == "org.Mm.eg.db"
) "Mus musculus" else "Homo sapiens"

hallmark_db <- msigdbr(
    species = species_name,
    collection = "H"
)

# Map expression gene IDs to symbols
if (GENE_ID_TYPE == "SYMBOL") {

    gene_symbols <- rownames(vst_matrix)

} else {

    gene_map <- AnnotationDbi::mapIds(
        OrgDb,
        keys = rownames(vst_matrix),
        column = "SYMBOL",
        keytype = GENE_ID_TYPE,
        multiVals = "first"
    )

    gene_symbols <- as.character(gene_map)
}

names(gene_symbols) <- rownames(vst_matrix)

annotation_hallmark <- data.frame(
    Treatment = metadata$Treatment,
    row.names = rownames(metadata)
)



############################################################
# 21. PCA
############################################################

pca_data <- plotPCA(
    vsd,
    intgroup = "Treatment",
    returnData = TRUE
)

percent_var <- round(
    100 * attr(
        pca_data,
        "percentVar"
    )
)


# # Add actual biological sample names
# pca_data$Name <- metadata[
#     rownames(pca_data),
#     "Name"
# ]


p <- ggplot(
    pca_data,
    aes(
        x = PC1,
        y = PC2,
        color = Treatment
    )
) +
    geom_point(size = 4) +
    xlab(
        paste0(
            "PC1: ",
            percent_var[1],
            "% variance"
        )
    ) +
    ylab(
        paste0(
            "PC2: ",
            percent_var[2],
            "% variance"
        )
    ) +
    theme_bw() +
    labs(
        title = "PCA of VST Gene Expression"
    )

ggsave(
    "DESeq2_results/PCA/PCA_Group.png",
    p,
    width = 8,
    height = 6,
    dpi = 300
)


write.csv(
    pca_data,
    "DESeq2_results/PCA/PCA_coordinates.csv",
    row.names = FALSE
)


############################################################
# 21B. BATCH-ADJUSTED PCA
############################################################

if (!is.null(vst_matrix_batch_adjusted)) {

    cat("\nGenerating batch-adjusted PCA...\n")

    # Select the 500 most variable genes
    gene_variances <- apply(
        vst_matrix_batch_adjusted,
        1,
        var
    )

    top_pca_genes <- names(sort(
        gene_variances,
        decreasing = TRUE
    ))[seq_len(min(500, length(gene_variances)))]

    # PCA using batch-adjusted VST values
    pca_adjusted <- prcomp(
        t(vst_matrix_batch_adjusted[
            top_pca_genes,
            ,
            drop = FALSE
        ]),
        center = TRUE,
        scale. = FALSE
    )

    percent_var_adjusted <- round(
        100 * pca_adjusted$sdev^2 /
        sum(pca_adjusted$sdev^2),
        1
    )

    # Ensure metadata matches PCA sample order
    pca_samples <- rownames(pca_adjusted$x)

    pca_adjusted_df <- data.frame(
        sample = pca_samples,
        PC1 = pca_adjusted$x[, 1],
        PC2 = pca_adjusted$x[, 2],
        Treatment = metadata[pca_samples, "Treatment"],
        Batch = metadata[pca_samples, "Batch"],
        row.names = NULL
    )

    write.csv(
        pca_adjusted_df,
        "DESeq2_results/PCA/PCA_Batch_Adjusted_Coordinates.csv",
        row.names = FALSE
    )

    p_adjusted <- ggplot(
        pca_adjusted_df,
        aes(
            x = PC1,
            y = PC2,
            color = Treatment
        )
    ) +
        geom_point(size = 4) +
        xlab(
            paste0(
                "PC1: ",
                percent_var_adjusted[1],
                "% variance"
            )
        ) +
        ylab(
            paste0(
                "PC2: ",
                percent_var_adjusted[2],
                "% variance"
            )
        ) +
        theme_bw() +
        labs(
            title = "PCA of Batch-Adjusted VST Expression"
        )

    ggsave(
        "DESeq2_results/PCA/PCA_Batch_Adjusted.png",
        p_adjusted,
        width = 8,
        height = 6,
        dpi = 300
    )
}

############################################################
# 22. SAMPLE CORRELATION HEATMAP
############################################################

sample_cor <- cor(
    assay(vsd),
    method = "pearson"
)


annotation <- data.frame(
    Group = metadata$Treatment
)

rownames(annotation) <- metadata$sample


png(
    "DESeq2_results/QC/Sample_Correlation_Heatmap.png",
    width = 2200,
    height = 2000,
    res = 250
)

pheatmap(
    sample_cor,
    annotation_col = annotation,
    annotation_row = annotation,
    main = "Sample Correlation",
    fontsize = 7,
    fontsize_row = 5,
    fontsize_col = 5,
    fontsize_main = 10,
    fontsize_legend = 7,
    fontsize_annotation = 7
)

dev.off()


############################################################
# 23. SAMPLE DISTANCE HEATMAP
############################################################

sample_dist <- dist(
    t(assay(vsd))
)

sample_dist_matrix <- as.matrix(
    sample_dist
)


png(
    "DESeq2_results/QC/Sample_Distance_Heatmap.png",
    width = 2200,
    height = 2000,
    res = 250
)

pheatmap(
    sample_dist_matrix,
    annotation_col = annotation,
    annotation_row = annotation,
    main = "Sample-to-Sample Distance",
    fontsize = 7,
    fontsize_row = 5,
    fontsize_col = 5,
    fontsize_main = 10,
    fontsize_legend = 7,
    fontsize_annotation = 7
)

dev.off()


############################################################
# 24. LOAD PAIRWISE COMPARISONS FROM CONFIG
############################################################

comparisons <- lapply(
    COMPARISONS,
    function(x){
        c(
            x$group1,
            x$group2
        )
    }
)

cat("\nPairwise comparisons:\n")
print(comparisons)









############################################################
# LOAD MAJOR PATHWAY GENE LIST
############################################################

pathway_genes <- readxl::read_excel(
    "mtInfGenes.xlsx",
    sheet = "mtInfGenes"
)


pathway_genes <- pathway_genes %>%
    dplyr::filter(
        !is.na(HumanGene),
        HumanGene != "",
        !is.na(`Major Pathways`),
        `Major Pathways` != ""
    ) %>%
    dplyr::mutate(
        HumanGene = trimws(HumanGene),
        `Major Pathways` = trimws(`Major Pathways`)
    )


############################################################
# DETERMINE ANALYSIS SPECIES
############################################################

if(ORGDB_PACKAGE == "org.Hs.eg.db"){

    ANALYSIS_SPECIES <- "human"

} else if(ORGDB_PACKAGE == "org.Mm.eg.db"){

    ANALYSIS_SPECIES <- "mouse"

} else {

    stop(
        paste(
            "Major pathway heatmaps currently support",
            "org.Hs.eg.db or org.Mm.eg.db.",
            "Received:",
            ORGDB_PACKAGE
        )
    )
}

cat(
    "\nAnalysis species:",
    ANALYSIS_SPECIES,
    "\n"
)


############################################################
# CREATE PATHWAY GENE MAP
############################################################

if(ANALYSIS_SPECIES == "human"){

    ########################################################
    # HUMAN: USE HUMAN SYMBOLS DIRECTLY
    ########################################################

    pathway_gene_map <- pathway_genes %>%
        dplyr::transmute(
            HumanGene = HumanGene,
            AnalysisGene = HumanGene,
            `Major Pathways` = `Major Pathways`
        )


} else if (ANALYSIS_SPECIES == "mouse") {

    ########################################################
    # MOUSE: MAP HUMAN GENES TO MOUSE ORTHOLOGS
    # Direct Ensembl BioMart API (bypasses biomaRt bug)
    ########################################################

    cat("\nMapping human pathway genes to mouse orthologs...\n")

    biomart_url <- "https://asia.ensembl.org/biomart/martservice"

    output_dir <- "DESeq2_results/Heatmaps/Major_Pathways"

    dir.create(
        output_dir,
        recursive = TRUE,
        showWarnings = FALSE
    )

    ########################################################
    # FUNCTION: QUERY ONE HUMAN GENE
    ########################################################

    get_mouse_ortholog <- function(human_gene, max_attempts = 3) {

        human_gene <- trimws(as.character(human_gene))

        # Escape XML-sensitive characters
        human_gene_xml <- gsub("&", "&amp;", human_gene, fixed = TRUE)
        human_gene_xml <- gsub('"', "&quot;", human_gene_xml, fixed = TRUE)
        human_gene_xml <- gsub("<", "&lt;", human_gene_xml, fixed = TRUE)
        human_gene_xml <- gsub(">", "&gt;", human_gene_xml, fixed = TRUE)

        xml_query <- paste0(
            '<?xml version="1.0" encoding="UTF-8"?>',
            '<!DOCTYPE Query>',
            '<Query virtualSchemaName="default" ',
            'formatter="TSV" header="1" uniqueRows="1" ',
            'datasetConfigVersion="0.6">',
            '<Dataset name="hsapiens_gene_ensembl" ',
            'interface="default">',
            '<Filter name="hgnc_symbol" value="',
            human_gene_xml,
            '"/>',
            '<Attribute name="mmusculus_homolog_associated_gene_name"/>',
            '</Dataset></Query>'
        )

        for (attempt in seq_len(max_attempts)) {

            response <- tryCatch(
                httr::GET(
                    biomart_url,
                    query = list(query = xml_query),
                    httr::timeout(60)
                ),
                error = function(e) e
            )

            if (!inherits(response, "error") &&
                httr::status_code(response) == 200) {

                output <- httr::content(
                    response,
                    as = "text",
                    encoding = "UTF-8"
                )

                # HTTP 200 can still contain a BioMart error
                if (!grepl("Query ERROR", output, fixed = TRUE)) {

                    genes <- tryCatch(
                        read.delim(
                            text = output,
                            check.names = FALSE,
                            stringsAsFactors = FALSE
                        ),
                        error = function(e) e
                    )

                    if (!inherits(genes, "error") &&
                        ncol(genes) >= 1) {

                        mouse_genes <- unique(
                            trimws(as.character(genes[[1]]))
                        )

                        mouse_genes <- mouse_genes[
                            !is.na(mouse_genes) &
                            nzchar(mouse_genes)
                        ]

                        return(mouse_genes)
                    }
                }
            }

            cat(
                "Retrying", human_gene,
                "- attempt", attempt, "of", max_attempts, "\n"
            )

            if (attempt < max_attempts) {
                Sys.sleep(2 * attempt)
            }
        }

        stop(
            "BioMart query failed for ", human_gene,
            " after ", max_attempts, " attempts. ",
            "Stopping to avoid silently losing pathway genes."
        )
    }

    ########################################################
    # CACHE SUCCESSFUL LOOKUPS
    ########################################################

    cache_file <- file.path(
        output_dir,
        "Ortholog_Query_Cache.rds"
    )

    if (file.exists(cache_file)) {
        ortholog_cache <- readRDS(cache_file)
    } else {
        ortholog_cache <- list()
    }

    human_genes <- unique(pathway_genes$HumanGene)

    ortholog_list <- vector("list", length(human_genes))

    for (i in seq_along(human_genes)) {

        human_gene <- human_genes[i]

        if (human_gene %in% names(ortholog_cache)) {

            mouse_genes <- ortholog_cache[[human_gene]]

        } else {

            mouse_genes <- get_mouse_ortholog(human_gene)

            ortholog_cache[human_gene] <- list(mouse_genes)

            saveRDS(ortholog_cache, cache_file)
        }

        ortholog_list[[i]] <- data.frame(
            HumanGene = human_gene,
            MouseGene = if (length(mouse_genes) > 0) {
                mouse_genes
            } else {
                NA_character_
            },
            stringsAsFactors = FALSE
        )

        if (i %% 25 == 0 || i == length(human_genes)) {
            cat(
                "Ortholog mapping progress:",
                i, "/", length(human_genes), "\n"
            )
        }
    }

    ortholog_map <- dplyr::bind_rows(ortholog_list)

    ########################################################
    # REMOVE EMPTY MAPPINGS
    ########################################################

    ortholog_map <- ortholog_map %>%
        dplyr::filter(
            !is.na(MouseGene),
            MouseGene != ""
        ) %>%
        dplyr::distinct(HumanGene, MouseGene)

    ########################################################
    # ATTACH PATHWAY INFORMATION
    ########################################################

    pathway_gene_map <- pathway_genes %>%
        dplyr::left_join(
            ortholog_map,
            by = "HumanGene"
        ) %>%
        dplyr::filter(
            !is.na(MouseGene),
            MouseGene != ""
        ) %>%
        dplyr::transmute(
            HumanGene = HumanGene,
            AnalysisGene = MouseGene,
            `Major Pathways` = `Major Pathways`
        )

    ########################################################
    # SAVE HUMAN -> MOUSE ORTHOLOG MAP
    ########################################################

    write.csv(
        pathway_gene_map,
        file.path(output_dir, "Human_to_Mouse_Ortholog_Map.csv"),
        row.names = FALSE
    )

    ########################################################
    # SAVE UNMAPPED HUMAN GENES
    ########################################################

    unmapped_human_genes <- pathway_genes %>%
        dplyr::anti_join(
            ortholog_map,
            by = "HumanGene"
        )

    write.csv(
        unmapped_human_genes,
        file.path(output_dir, "Unmapped_Human_Pathway_Genes.csv"),
        row.names = FALSE
    )

    ########################################################
    # MAPPING SUMMARY
    ########################################################

    cat(
        "\nHuman genes requested:",
        length(human_genes), "\n"
    )

    cat(
        "Human genes with mouse orthologs:",
        length(unique(ortholog_map$HumanGene)), "\n"
    )

    cat(
        "Unique mouse orthologs:",
        length(unique(ortholog_map$MouseGene)), "\n"
    )

    cat(
        "Mouse orthologs present in filtered expression data:",
        sum(
            toupper(unique(ortholog_map$MouseGene)) %in%
            toupper(rownames(vst_matrix))
        ),
        "\n"
    )

    if (nrow(ortholog_map) == 0) {
        stop("No human-to-mouse orthologs were mapped.")
    }
}

############################################################
# 25. FUNCTION FOR DIFFERENTIAL EXPRESSION
############################################################

run_comparison <- function(group1, group2){

    comparison_name <- paste0(
        group1,
        "_vs_",
        group2
    )
    comparison_name <- gsub("[^A-Za-z0-9._-]", "_", comparison_name)
    cat(
        "\nRunning:",
        comparison_name,
        "\n"
    )



    

    ########################################################
    # DESEQ2 RESULTS
    ########################################################

    res <- results(
        dds,
        contrast = c(
            "Treatment",
            group1,
            group2
        ),
        alpha = 0.05
    )


    res_df <- as.data.frame(res)

    res_df$Gene <- rownames(res_df)

    ########################################################
    # ADD GENE SYMBOL
    ########################################################
        
    if(GENE_ID_TYPE == "SYMBOL"){
    
        res_df$Symbol <- res_df$Gene
    
    } else {
    
        symbol_conversion <- bitr(
            res_df$Gene,
            fromType = GENE_ID_TYPE,
            toType = "SYMBOL",
            OrgDb = OrgDb
        ) %>%
            distinct(
                .data[[GENE_ID_TYPE]],
                .keep_all = TRUE
            ) %>%
            rename(
                Symbol = SYMBOL
            )
    
        res_df <- res_df %>%
            left_join(
                symbol_conversion,
                by = setNames(
                    GENE_ID_TYPE,
                    "Gene"
                )
            )
    }
            
    res_df <- res_df %>%
        dplyr::select(
            Gene,
            Symbol,
            baseMean,
            log2FoldChange,
            lfcSE,
            stat,
            pvalue,
            padj
        ) %>%
        dplyr::arrange(padj)


    ########################################################
    # SAVE ALL RESULTS
    ########################################################

    write.csv(
        res_df,
        paste0(
            "DESeq2_results/DGE/",
            comparison_name,
            "_All_Genes.csv"
        ),
        row.names = FALSE
    )


    ########################################################
    # SIGNIFICANT GENES
    ########################################################

    sig <- res_df %>%
        filter(
            !is.na(padj),
            padj < 0.05,
            abs(log2FoldChange) >= 1
        )


    up <- sig %>%
        filter(
            log2FoldChange >= 1
        )


    down <- sig %>%
        filter(
            log2FoldChange <= -1
        )


    write.csv(
        sig,
        paste0(
            "DESeq2_results/DGE/",
            comparison_name,
            "_Significant.csv"
        ),
        row.names = FALSE
    )


    write.csv(
        up,
        paste0(
            "DESeq2_results/DGE/",
            comparison_name,
            "_Upregulated.csv"
        ),
        row.names = FALSE
    )


    write.csv(
        down,
        paste0(
            "DESeq2_results/DGE/",
            comparison_name,
            "_Downregulated.csv"
        ),
        row.names = FALSE
    )


    cat(
        "Significant:",
        nrow(sig),
        "\n"
    )

    cat(
        "Upregulated:",
        nrow(up),
        "\n"
    )

    cat(
        "Downregulated:",
        nrow(down),
        "\n"
    )


    ########################################################
    # MA PLOT
    ########################################################

    png(
        paste0(
            "DESeq2_results/DGE/",
            comparison_name,
            "_MAplot.png"
        ),
        width = 1800,
        height = 1500,
        res = 250
    )

    plotMA(
        res,
        ylim = c(-5, 5),
        main = comparison_name
    )

    dev.off()


    ########################################################
    # VOLCANO PLOT
    ########################################################

    png(
        paste0(
            "DESeq2_results/Volcano/",
            comparison_name,
            "_Volcano.png"
        ),
        width = 2200,
        height = 1800,
        res = 250
    )

    print(
        EnhancedVolcano(
            res_df,
            lab = ifelse(is.na(res_df$Symbol) | res_df$Symbol == "",res_df$Gene,res_df$Symbol),
            x = "log2FoldChange",
            y = "padj",
            pCutoff = 0.05,
            FCcutoff = 1,
            title = comparison_name,
            subtitle = paste(
                group1,
                "vs",
                group2
            )
        )
    )

    dev.off()


    ########################################################
    # TOP 50 HEATMAP
    ########################################################
    
    top_genes <- res_df %>%
        filter(
            !is.na(padj),
            padj < 0.05,
            abs(log2FoldChange) >= 1
        ) %>%
        arrange(padj) %>%
        head(50) %>%
        pull(Gene)
    
    top_genes <- top_genes[
        top_genes %in% rownames(vst_matrix)
    ]
    
    
    if(length(top_genes) >= 2){
    

        # Keep group1 samples first, followed by group2.
        # Preserve original metadata order within each group.
      
        comparison_samples <- c(
          rownames(metadata)[metadata$Treatment == group1],
          rownames(metadata)[metadata$Treatment == group2]
        )
      
          
        # Subset expression matrix

        hm <- heatmap_matrix[
            top_genes,
            comparison_samples,
            drop = FALSE
        ]

        
        ########################################################
        # USE GENE SYMBOLS AS HEATMAP ROW LABELS
        ########################################################
        
        heatmap_labels <- res_df$Symbol[
            match(
                rownames(hm),
                res_df$Gene
            )
        ]
        
        missing_labels <- is.na(heatmap_labels) |
                          heatmap_labels == ""
        
        heatmap_labels[missing_labels] <- rownames(hm)[
            missing_labels
        ]
        
        rownames(hm) <- make.unique(heatmap_labels)
    
        # Subset annotation
        annotation_comp <- annotation[
            comparison_samples,
            ,
            drop = FALSE
        ]
    
        png(
            paste0(
                "DESeq2_results/Heatmaps/",
                comparison_name,
                "_Top50.png"
            ),
            width = 2500,
            height = 2400,
            res = 250
        )
    
        pheatmap(
            hm,
            scale = "row",
            annotation_col = annotation_comp,
            cluster_cols = FALSE,
            show_rownames = TRUE,
            fontsize_row = 6,
            main = paste(
                "Top 50 DE Genes:",
                comparison_name
            )
        )
    
        dev.off()
    }





    ########################################################
    # HALLMARK HEATMAPS WITH SIGNIFICANCE MARKERS
    ########################################################
    

    # Keep group1 samples first, followed by group2.
    # Preserve original metadata order within each group.

    comparison_samples <- c(
        rownames(metadata)[metadata$Treatment == group1],
        rownames(metadata)[metadata$Treatment == group2]
    )

    
    hallmark_dir <- file.path(
        "DESeq2_results/Heatmaps/Hallmark",
        comparison_name
    )
    
    dir.create(
        hallmark_dir,
        recursive = TRUE,
        showWarnings = FALSE
    )
    
    for (pathway in HALLMARK_PATHWAYS) {
    
        pathway_symbols <- unique(
            hallmark_db$gene_symbol[
                hallmark_db$gs_name == pathway
            ]
        )
    
        matched_genes <- names(gene_symbols)[
            !is.na(gene_symbols) &
            gene_symbols %in% pathway_symbols
        ]
    
        if (length(matched_genes) < 2) next
    

        hm <- heatmap_matrix[
            matched_genes,
            comparison_samples,
            drop = FALSE
        ]
 
    
        # Remove zero-variance genes
        keep <- apply(hm, 1, var) > 0
        keep[is.na(keep)] <- FALSE
        hm <- hm[keep, , drop = FALSE]
    
        if (nrow(hm) < 2) next
    
        # Match DESeq2 results by original gene ID
        gene_stats <- res_df[
            match(rownames(hm), res_df$Gene),
            ,
            drop = FALSE
        ]
    
        sig_marker <- rep("", nrow(hm))
    
        sig_marker[
            !is.na(gene_stats$padj) &
            gene_stats$padj < 0.05
        ] <- " *"
    
        sig_marker[
            !is.na(gene_stats$padj) &
            gene_stats$padj < 0.05 &
            !is.na(gene_stats$log2FoldChange) &
            abs(gene_stats$log2FoldChange) >= 1
        ] <- " **"
    
        labels <- gene_symbols[rownames(hm)]
    
        missing <- is.na(labels) | labels == ""
        labels[missing] <- rownames(hm)[missing]
    
        rownames(hm) <- make.unique(
            paste0(labels, sig_marker)
        )
    
        annotation_comp <- data.frame(
            Treatment = metadata[comparison_samples, "Treatment"],
            row.names = comparison_samples
        )
    
        png(
            file.path(hallmark_dir, paste0(pathway, ".png")),
            width = 2500,
            height = max(1800, nrow(hm) * 35),
            res = 250
        )
    
        pheatmap(
            hm,
            scale = "row",
            annotation_col = annotation_comp,
            cluster_cols = FALSE,
            show_rownames = TRUE,
            fontsize_row = 6,
            main = paste0(
                pathway, " - ", comparison_name,
                "\n* padj < 0.05; ** padj < 0.05 & |log2FC| >= 1"
            )
        )
    
        dev.off()
    }




    ########################################################
    # MAJOR PATHWAY HEATMAPS
    ########################################################
    

    comparison_samples <- c(
        rownames(metadata)[metadata$Treatment == group1],
        rownames(metadata)[metadata$Treatment == group2]
    )

    
    
    ########################################################
    # CREATE EXPRESSION GENE -> SYMBOL MAP
    ########################################################
    
    expression_gene_map <- res_df %>%
        dplyr::select(
            Gene,
            Symbol
        ) %>%
        dplyr::filter(
            !is.na(Symbol),
            Symbol != ""
        ) %>%
        dplyr::mutate(
            Symbol_match = toupper(Symbol)
        )
    
    
    ########################################################
    # LOOP THROUGH MAJOR PATHWAYS
    ########################################################
    
    major_pathways <- unique(
        pathway_gene_map$`Major Pathways`
    )
    
    
    for(pathway in major_pathways){
    
        cat(
            "\nCreating major pathway heatmap:",
            pathway,
            "\n"
        )
    
    
        ####################################################
        # PATHWAY GENES
        ####################################################
    
        pathway_symbols <- pathway_gene_map %>%
            dplyr::filter(
                `Major Pathways` == pathway
            ) %>%
            dplyr::pull(
                AnalysisGene
            ) %>%
            unique()
    
    
        ####################################################
        # CASE-INSENSITIVE SYMBOL MATCH
        ####################################################
    
        pathway_symbols_match <- toupper(
            pathway_symbols
        )
    
        pathway_matches <- expression_gene_map %>%
            dplyr::filter(
                Symbol_match %in% pathway_symbols_match,
                Gene %in% rownames(vst_matrix)
            ) %>%
            dplyr::distinct(
                Gene,
                .keep_all = TRUE
            )
    
    
        ####################################################
        # ADD DE STATISTICS
        ####################################################
    
        pathway_matches <- pathway_matches %>%
            dplyr::left_join(
                res_df %>%
                    dplyr::select(
                        Gene,
                        log2FoldChange,
                        pvalue,
                        padj
                    ),
                by = "Gene"
            )
    
    
        ####################################################
        # REQUIRE AT LEAST 2 EXPRESSED GENES
        ####################################################
    
        if(nrow(pathway_matches) < 2){
    
            cat(
                "Skipping",
                pathway,
                "- fewer than 2 genes found in expression data.\n"
            )
    
            next
        }
    
    
        ####################################################
        # SUBSET VST MATRIX
        ####################################################
    

        hm_pathway <- heatmap_matrix[
            pathway_matches$Gene,
            comparison_samples,
            drop = FALSE
        ]

    
    
        ####################################################
        # SYMBOL ROW LABELS
        ####################################################
    
        pathway_labels <- pathway_matches$Symbol[
            match(
                rownames(hm_pathway),
                pathway_matches$Gene
            )
        ]
    
        missing_labels <- is.na(pathway_labels) |
                          pathway_labels == ""
    
        pathway_labels[missing_labels] <-
            rownames(hm_pathway)[
                missing_labels
            ]
    
    
        ####################################################
        # GET DE STATS IN HEATMAP ROW ORDER
        ####################################################
    
        heatmap_stats <- pathway_matches[
            match(
                rownames(hm_pathway),
                pathway_matches$Gene
            ),
            ,
            drop = FALSE
        ]
    
    
        ####################################################
        # ADD SIGNIFICANCE MARKERS
        #
        # *  = padj < 0.05
        # ** = padj < 0.05 AND |log2FC| >= 1
        ####################################################
    
        sig_marker <- rep(
            "",
            nrow(heatmap_stats)
        )
    
        sig_marker[
            !is.na(heatmap_stats$padj) &
            heatmap_stats$padj < 0.05
        ] <- " *"
    
        sig_marker[
            !is.na(heatmap_stats$padj) &
            heatmap_stats$padj < 0.05 &
            abs(heatmap_stats$log2FoldChange) >= 1
        ] <- " **"
    
    
        ####################################################
        # ADD MARKERS TO ROW LABELS
        ####################################################
    
        pathway_labels <- paste0(
            pathway_labels,
            sig_marker
        )
    
        rownames(hm_pathway) <- make.unique(
            pathway_labels
        )
    
    
        ####################################################
        # REMOVE ZERO-VARIANCE GENES
        #
        # scale="row" cannot meaningfully Z-score genes
        # with no variance across samples.
        ####################################################
    
        pathway_variance <- apply(
            hm_pathway,
            1,
            var
        )
    
        keep_variable <- !is.na(pathway_variance) &
                         pathway_variance > 0
    
        hm_pathway <- hm_pathway[
            keep_variable,
            ,
            drop = FALSE
        ]
    
    
        if(nrow(hm_pathway) < 2){
    
            cat(
                "Skipping",
                pathway,
                "- fewer than 2 variable genes remain.\n"
            )
    
            next
        }
    
    
        ####################################################
        # SAMPLE ANNOTATION
        ####################################################
    
        annotation_comp <- annotation[
            comparison_samples,
            ,
            drop = FALSE
        ]
    
    
        ####################################################
        # CLEAN PATHWAY NAME
        ####################################################
    
        pathway_file <- gsub(
            "[^A-Za-z0-9._-]",
            "_",
            pathway
        )
    
    
        ####################################################
        # ADD SIGNIFICANCE CATEGORY TO OUTPUT TABLE
        ####################################################
    
        pathway_output_df <- pathway_matches %>%
            dplyr::mutate(
                Significance = dplyr::case_when(
                    !is.na(padj) &
                    padj < 0.05 &
                    abs(log2FoldChange) >= 1 ~
                        "** padj < 0.05 & |log2FC| >= 1",
    
                    !is.na(padj) &
                    padj < 0.05 ~
                        "* padj < 0.05",
    
                    TRUE ~
                        "Not significant"
                )
            ) %>%
            dplyr::select(
                Gene,
                Symbol,
                log2FoldChange,
                pvalue,
                padj,
                Significance
            )
    
    
        ########################################################
        # CREATE COMPARISON-SPECIFIC PATHWAY DIRECTORY
        ########################################################
        
        comparison_pathway_dir <- file.path(
            "DESeq2_results/Heatmaps/Major_Pathways",
            comparison_name
        )
        
        dir.create(
            comparison_pathway_dir,
            showWarnings = FALSE,
            recursive = TRUE
        )       
    
        ####################################################
        # SAVE PATHWAY GENE LIST USED
        ####################################################
    
        write.csv(
            pathway_output_df,
            file.path(
                comparison_pathway_dir,
                paste0(
                    pathway_file,
                    "_Genes_Used.csv"
                )
            ),
            row.names = FALSE
        )
    
    
        ####################################################
        # SAVE HEATMAP
        ####################################################
    
        png(
        file.path(
            comparison_pathway_dir,
            paste0(
                pathway_file,
                "_Heatmap.png"
            )
        ),
            width = 2500,
            height = max(
                1800,
                nrow(hm_pathway) * 35
            ),
            res = 250
        )
    
        pheatmap(
            hm_pathway,
            scale = "row",
            annotation_col = annotation_comp,
            cluster_cols = FALSE,
            show_rownames = TRUE,
            fontsize_row = 6,
            main = paste0(
                pathway,
                " - ",
                comparison_name,
                "\n* padj < 0.05; ** padj < 0.05 & |log2FC| >= 1"
            )
        )
    
        dev.off()
    }




    ########################################################
    # TOP 50 UP / DOWN PATHWAY ENRICHMENT
    ########################################################
    
    # Select top 50 significant UP genes
    top_up_genes <- up %>%
        arrange(padj) %>%
        head(50) %>%
        pull(Gene)
    
    # Select top 50 significant DOWN genes
    top_down_genes <- down %>%
        arrange(padj) %>%
        head(50) %>%
        pull(Gene)
    
    

    
    

    
    ########################################################
    # FUNCTION TO RUN ENRICHR
    ########################################################
    
    run_enrichr <- function(gene_vector, direction){
    
    
        ####################################################
        # CONVERT TO GENE SYMBOLS
        ####################################################
    
        if(GENE_ID_TYPE == "SYMBOL"){
    
            enrichr_genes <- unique(gene_vector)
    
        } else {
    
            gene_conversion <- bitr(
                gene_vector,
                fromType = GENE_ID_TYPE,
                toType = "SYMBOL",
                OrgDb = OrgDb
            )
    
            enrichr_genes <- unique(
                gene_conversion$SYMBOL
            )
        }
    
    
        ####################################################
        # VERIFY ENOUGH GENES
        ####################################################
    
        if(length(enrichr_genes) < 5){
    
            cat(
                "\nSkipping",
                direction,
                "pathway enrichment - fewer than 5 mapped genes.\n"
            )
    
            return(NULL)
        }
    
    
        cat(
            "\nRunning Top",
            length(enrichr_genes),
            direction,
            "pathway enrichment...\n"
        )
    
    
        ####################################################
        # RUN ENRICHR
        ####################################################
    
        pathway_results <- enrichr(
            enrichr_genes,
            ENRICHR_DATABASES
        )
    
    
        ####################################################
        # SAVE EACH DATABASE AUTOMATICALLY
        ####################################################
    
        for(database in ENRICHR_DATABASES){
    
            result_df <- pathway_results[[database]]
    
            database_file <- gsub(
                "[^A-Za-z0-9_-]",
                "_",
                database
            )
    
            output_file <- paste0(
                "DESeq2_results/GSEA/",
                comparison_name,
                "_Top50_",
                direction,
                "_",
                database_file,
                ".csv"
            )
    
            write.csv(
                result_df,
                output_file,
                row.names = FALSE
            )
        }
    }
    
    
    ########################################################
    # RUN UP AND DOWN ENRICHMENT
    ########################################################
    
    run_enrichr(
        top_up_genes,
        "UP"
    )
    
    run_enrichr(
        top_down_genes,
        "DOWN"
    )
    
    
    
    ########################################################
    # GSEA - ALL GENES
    ########################################################
    
    ranked_df <- res_df %>%
        filter(
            !is.na(stat),
            !is.na(Gene)
        )
    
    ########################################################
    # CONVERT INPUT GENE IDs -> ENTREZ
    ########################################################
    
    gene_conversion <- bitr(
        ranked_df$Gene,
        fromType = GENE_ID_TYPE,
        toType = "ENTREZID",
        OrgDb = OrgDb
    )

    ########################################################
    # JOIN ENTREZ IDS TO RANKED GENES
    ########################################################
        
    gsea_df <- ranked_df %>%
        inner_join(
            gene_conversion,
            by = setNames(
                GENE_ID_TYPE,
                "Gene"
            )
        ) %>%
        arrange(
            desc(abs(stat))
        ) %>%
        distinct(
            ENTREZID,
            .keep_all = TRUE
        )
    
    ########################################################
    # CREATE RANKED GENE LIST
    ########################################################
    
    gene_list <- gsea_df$stat
    
    names(gene_list) <- gsea_df$ENTREZID
    
    gene_list <- sort(
        gene_list,
        decreasing = TRUE
    )
    
    
    ########################################################
    # RUN GO GSEA
    ########################################################
    
    if(length(gene_list) > 100){
    
        go_ontologies <- c(
            "BP",
            "CC",
            "MF"
        )
    
    
        for(go_ont in go_ontologies){
    
            cat(
                "\nRunning GSEA GO",
                go_ont,
                "for",
                comparison_name,
                "\n"
            )
    
    
            gsea_go <- gseGO(
                geneList = gene_list,
                OrgDb = OrgDb,
                ont = go_ont,
                keyType = "ENTREZID",
                minGSSize = 10,
                maxGSSize = 500,
                pvalueCutoff = 0.05,
                verbose = FALSE
            )
    
    
            ################################################
            # SAVE GSEA RESULTS
            ################################################
    
            write.csv(
                as.data.frame(gsea_go),
                paste0(
                    "DESeq2_results/GSEA/",
                    comparison_name,
                    "_GSEA_GO_",
                    go_ont,
                    ".csv"
                ),
                row.names = FALSE
            )
    
    
            ################################################
            # GSEA DOTPLOT
            ################################################
    
            if(
                !is.null(gsea_go) &&
                nrow(as.data.frame(gsea_go)) > 0
            ){
    
                png(
                    paste0(
                        "DESeq2_results/GSEA/",
                        comparison_name,
                        "_GSEA_GO_",
                        go_ont,
                        "_Dotplot.png"
                    ),
                    width = 2600,
                    height = 1800,
                    res = 250
                )
    
                print(
                    dotplot(
                        gsea_go,
                        showCategory = 20,
                        font.size = 8
                    ) +
                    ggtitle(
                        paste(
                            "GSEA GO",
                            go_ont,
                            "-",
                            comparison_name
                        )
                    ) +
                    theme(
                        axis.text.y = element_text(size = 7),
                        axis.text.x = element_text(size = 8),
                        plot.title = element_text(size = 10)
                    )
                )
                    
                dev.off()
            }
        }
}}


############################################################
# 26. RUN ALL PAIRWISE GROUP COMPARISONS
############################################################

for(comp in comparisons){

    run_comparison(
        comp[1],
        comp[2]
    )
}


############################################################
# 27. SESSION INFO
############################################################

sink(
    "DESeq2_results/sessionInfo.txt"
)

sessionInfo()

sink()


############################################################
# 28. DONE
############################################################

cat("\n========================================\n")
cat("ANALYSIS COMPLETE\n")
cat("========================================\n")

cat(
    "\nSamples:",
    ncol(counts),
    "\n"
)

cat(
    "Genes before filtering:",
    nrow(counts),
    "\n"
)

cat(
    "Genes after filtering:",
    nrow(counts_filtered),
    "\n"
)

cat(
    "Groups:",
    paste(
        levels(metadata$Treatment),
        collapse = ", "
    ),
    "\n"
)

cat(
    "Comparisons:",
    length(comparisons),
    "\n"
)

cat(
    "\nResults saved in: DESeq2_results/\n"
)
