## Run MSstats for the 30 SPD runs.
## Need to locally load the dataset beforehand because it is enormous. This step will be eliminated in the future due to the PDP MSstats integration.
setwd("~/Desktop")
## also commented out all .RData file saves except for final one, to load into .Rmd file

#### 30 SPD first through MSstats ####
df_raw_30 <- read.table("Report_CardioCRISPR_FullScreen_Plates1-4_DIA-30SPD.tsv", sep = "\t", header = T, na.strings = c("NA", "NULL", "Null", "null", "NaN"))

## drop rows that MSstats throws out
df_raw_30 <- df_raw_30[which(df_raw_30$F.ExcludedFromQuantification=="False" & df_raw_30$F.FrgLossType=="noloss"),]

## Create annotations and run Converter
annot_precursor_30 <- unique(df_raw_30[,c("R.Condition", "R.FileName", "R.Replicate")])
annot_30 <- data.frame(
  Run = annot_precursor_30$R.FileName,
  Condition = annot_precursor_30$R.Condition,
  BioReplicate = annot_precursor_30$R.Replicate
)
## save memory
annot_precursor_30 <- NULL

quant_30 <- MSstats::SpectronauttoMSstatsFormat(df_raw_30,
                                                annotation = annot_30,
                                                filter_with_Qvalue = T,
                                                qvalue_cutoff = 0.01,
                                                removeProtein_with1Feature = T)
## save memory
df_raw_30 <- NULL
char_filename_date <- format(Sys.time(), "%Y_%m_%d_%H_%M_%S")
# save.image(paste0("./data_images/", char_filename_date, "quant_30.RData"))

## Run dataProcess step, "actual math" as Liang says
## If running Mac/Linux, can use multiple cores.
if (.Platform$OS.type=="unix") {
  spectronaut_proposed_30 <- MSstats::dataProcess(quant_30, 
                                                  normalization = 'EQUALIZEMEDIANS',
                                                  summaryMethod = "TMP",
                                                  # cutoffCensored = "minFeature",
                                                  censoredInt = "0",
                                                  ## suggested by Devon Kohler for large datasets
                                                  MBimpute = FALSE,
                                                  ## for DIA datasets, use topN
                                                  featureSubset = "topN",
                                                  ## for DIA datasets, use topN
                                                  n_top_feature = 20,
                                                  ## for MacOS or Linux, can assign multiple cores
                                                  numberOfCores = 4,
                                                  maxQuantileforCensored = 0.999)
}
if (.Platform$OS.type=="windows") {
  spectronaut_proposed_30 <- MSstats::dataProcess(quant_30, 
                                                  normalization = 'EQUALIZEMEDIANS',
                                                  summaryMethod = "TMP",
                                                  # cutoffCensored = "minFeature",
                                                  censoredInt = "0",
                                                  ## suggested by Devon Kohler for large datasets
                                                  MBimpute = FALSE,
                                                  ## for DIA datasets, use topN
                                                  featureSubset = "topN",
                                                  ## for DIA datasets, use topN
                                                  n_top_feature = 20,
                                                  ## for MacOS or Linux, can assign multiple cores
                                                  # numberOfCores = 4,
                                                  maxQuantileforCensored = 0.999)
}
## save memory, compost
quant_30 <- NULL
char_filename_date <- format(Sys.time(), "%Y_%m_%d_%H_%M_%S")
# save.image(paste0("./data_images/", char_filename_date, "spectronaut_proposed_30.RData"))

## assign protein level data a dataframe
df_30 <- spectronaut_proposed_30$ProteinLevelData


## an optional data summarization sanity check
# MSstats::dataProcessPlots(data = spectronaut_proposed_30, type = "QCplot", width = 40, height = 10, which.Protein="allonly")

## Next step for MSstats is groupComparison. Have to make the comparison matrix. Automated.
## make the controls, aka denominator
vec_control_conditions <- c("AAVS1_G1G2", "HPRT_G1G2", "HPRT_G1")
## create a list to collapse at end of building, one matrix for each control condition
## each element of this list will become a matrix, which we will rbind together at end
list_mat_comparison_msstats_30 <- vector("list", length= length(vec_control_conditions))
## make the experimental conditions
vec_experimental_conditions_30 <- levels(spectronaut_proposed_30$ProteinLevelData$GROUP)[!(levels(spectronaut_proposed_30$ProteinLevelData$GROUP) %in% vec_control_conditions)]
## how to automate this matrix, hmm

## First create matrices that are length(vec_experimental_conditions) rows long with rownames == vec_experimental_conditions and colnames == all.
for (i in 1:length(list_mat_comparison_msstats_30)) {
  ## for now, the rownames do not contain the name of the control condition. We will replace later after grepl step. 
  list_mat_comparison_msstats_30_names01 <- list(mat_comparison_rows = paste0(vec_experimental_conditions_30),
                                                 mat_comparison_cols = levels(spectronaut_proposed_30$ProteinLevelData$GROUP))
  list_mat_comparison_msstats_30[[i]] <- matrix(ncol = length(levels(spectronaut_proposed_30$ProteinLevelData$GROUP)),
                                                nrow = length(vec_experimental_conditions_30),
                                                dimnames = list_mat_comparison_msstats_30_names01)
}

## This is the major matrix creation step
for (j in 1:length(list_mat_comparison_msstats_30)) {
  for (i in 1:length(vec_experimental_conditions_30)) {
    a <- vec_experimental_conditions_30[i]
    b <- list_mat_comparison_msstats_30[[j]][i,]
    list_mat_comparison_msstats_30[[j]][i,] <- as.integer(grepl(a, names(b)))
  }
}

## Now that we have that, we can assign the controls a value of -1
for (i in 1:length(vec_control_conditions)) {
  list_mat_comparison_msstats_30[[i]][,grep(paste0("^", vec_control_conditions[i], "$"), colnames(list_mat_comparison_msstats_30[[i]]))] <- -1
}

## last step is to add the control to the rownames
for (i in 1:length(list_mat_comparison_msstats_30)) {
  list_mat_comparison_msstats_30_names02 <- paste0(vec_experimental_conditions_30, "-", vec_control_conditions[i])
  rownames(list_mat_comparison_msstats_30[[i]]) <-  list_mat_comparison_msstats_30_names02
  
}

## best to build each comparison matrix for each CONTROL at a time. Each control row gets a -1 as the denominator. ALso want to fill in the 

mat_comparison_msstats_30 <- do.call(rbind, list_mat_comparison_msstats_30)
## BOOM


## Run the differential analysis using MSstats
DA_msstats_30 <- MSstats::groupComparison(contrast.matrix = mat_comparison_msstats_30, data = spectronaut_proposed_30)

## save memory, compost
spectronaut_proposed_30 <- NULL
char_filename_date <- format(Sys.time(), "%Y_%m_%d_%H_%M_%S")
save.image(paste0("./data_images/", char_filename_date, "DA_msstats_30.RData"))

## Pull out the comparison table
df_comp_30 <- DA_msstats_30$ComparisonResult

## save memory, compost
DA_msstats_30 <- NULL
char_filename_date <- format(Sys.time(), "%Y_%m_%d_%H_%M_%S")
# save.image(paste0("./data_images/", char_filename_date, "df_30_df_comp_30.RData"))





## Run MSstats for 60 SPD runs.

#### 60 SPD first through MSstats ####
df_raw_60 <- read.table("Report_CardioCRISPR_FullScreen_Plates1-4_DIA-60SPD.tsv", sep = "\t", header = T, na.strings = c("NA", "NULL", "Null", "null", "NaN"))

## drop rows that MSstats throws out
df_raw_60 <- df_raw_60[which(df_raw_60$F.ExcludedFromQuantification=="False" & df_raw_60$F.FrgLossType=="noloss"),]

## Create annotations and run Converter
annot_precursor_60 <- unique(df_raw_60[,c("R.Condition", "R.FileName", "R.Replicate")])
annot_60 <- data.frame(
  Run = annot_precursor_60$R.FileName,
  Condition = annot_precursor_60$R.Condition,
  BioReplicate = annot_precursor_60$R.Replicate
)
## save memory
annot_precursor_60 <- NULL

quant_60 <- MSstats::SpectronauttoMSstatsFormat(df_raw_60,
                                                annotation = annot_60,
                                                filter_with_Qvalue = T,
                                                qvalue_cutoff = 0.01,
                                                removeProtein_with1Feature = T)
## save memory
df_raw_60 <- NULL
char_filename_date <- format(Sys.time(), "%Y_%m_%d_%H_%M_%S")
# save.image(paste0("./data_images/", char_filename_date, "quant_60.RData"))

## Run dataProcess step, "actual math" as Liang says
## If running Mac/Linux, can use multiple cores.
if (.Platform$OS.type=="unix") {
  spectronaut_proposed_60 <- MSstats::dataProcess(quant_60, 
                                                  normalization = 'EQUALIZEMEDIANS',
                                                  summaryMethod = "TMP",
                                                  # cutoffCensored = "minFeature",
                                                  censoredInt = "0",
                                                  ## suggested by Devon Kohler for large datasets
                                                  MBimpute = FALSE,
                                                  ## for DIA datasets, use topN
                                                  featureSubset = "topN",
                                                  ## for DIA datasets, use topN
                                                  n_top_feature = 20,
                                                  ## for MacOS or Linux, can assign multiple cores
                                                  numberOfCores = 4,
                                                  maxQuantileforCensored = 0.999)
}
if (.Platform$OS.type=="windows") {
  spectronaut_proposed_60 <- MSstats::dataProcess(quant_60, 
                                                  normalization = 'EQUALIZEMEDIANS',
                                                  summaryMethod = "TMP",
                                                  # cutoffCensored = "minFeature",
                                                  censoredInt = "0",
                                                  ## suggested by Devon Kohler for large datasets
                                                  MBimpute = FALSE,
                                                  ## for DIA datasets, use topN
                                                  featureSubset = "topN",
                                                  ## for DIA datasets, use topN
                                                  n_top_feature = 20,
                                                  ## for MacOS or Linux, can assign multiple cores
                                                  # numberOfCores = 4,
                                                  maxQuantileforCensored = 0.999)
}
## save memory, compost
quant_60 <- NULL
char_filename_date <- format(Sys.time(), "%Y_%m_%d_%H_%M_%S")
# save.image(paste0("./data_images/", char_filename_date, "spectronaut_proposed_60.RData"))

## assign protein level data a dataframe
df_60 <- spectronaut_proposed_60$ProteinLevelData


## an optional data summarization sanity check
# MSstats::dataProcessPlots(data = spectronaut_proposed_60, type = "QCplot", width = 40, height = 10, which.Protein="allonly")

## Next step for MSstats is groupComparison. Have to make the comparison matrix. Automated.
## make the controls, aka denominator
vec_control_conditions <- c("AAVS1_G1G2", "HPRT_G1G2", "HPRT_G1")
## create a list to collapse at end of building, one matrix for each control condition
## each element of this list will become a matrix, which we will rbind together at end
list_mat_comparison_msstats_60 <- vector("list", length= length(vec_control_conditions))
## make the experimental conditions
vec_experimental_conditions_60 <- levels(spectronaut_proposed_60$ProteinLevelData$GROUP)[!(levels(spectronaut_proposed_60$ProteinLevelData$GROUP) %in% vec_control_conditions)]
## how to automate this matrix, hmm

## First create matrices that are length(vec_experimental_conditions) rows long with rownames == vec_experimental_conditions and colnames == all.
for (i in 1:length(list_mat_comparison_msstats_60)) {
  ## for now, the rownames do not contain the name of the control condition. We will replace later after grepl step. 
  list_mat_comparison_msstats_60_names01 <- list(mat_comparison_rows = paste0(vec_experimental_conditions_60),
                                                 mat_comparison_cols = levels(spectronaut_proposed_60$ProteinLevelData$GROUP))
  list_mat_comparison_msstats_60[[i]] <- matrix(ncol = length(levels(spectronaut_proposed_60$ProteinLevelData$GROUP)),
                                                nrow = length(vec_experimental_conditions_60),
                                                dimnames = list_mat_comparison_msstats_60_names01)
}

## This is the major matrix creation step
for (j in 1:length(list_mat_comparison_msstats_60)) {
  for (i in 1:length(vec_experimental_conditions_60)) {
    a <- vec_experimental_conditions_60[i]
    b <- list_mat_comparison_msstats_60[[j]][i,]
    list_mat_comparison_msstats_60[[j]][i,] <- as.integer(grepl(a, names(b)))
  }
}

## Now that we have that, we can assign the controls a value of -1
for (i in 1:length(vec_control_conditions)) {
  list_mat_comparison_msstats_60[[i]][,grep(paste0("^", vec_control_conditions[i], "$"), colnames(list_mat_comparison_msstats_60[[i]]))] <- -1
}

## last step is to add the control to the rownames
for (i in 1:length(list_mat_comparison_msstats_60)) {
  list_mat_comparison_msstats_60_names02 <- paste0(vec_experimental_conditions_60, "-", vec_control_conditions[i])
  rownames(list_mat_comparison_msstats_60[[i]]) <-  list_mat_comparison_msstats_60_names02
  
}

## best to build each comparison matrix for each CONTROL at a time. Each control row gets a -1 as the denominator. ALso want to fill in the 

mat_comparison_msstats_60 <- do.call(rbind, list_mat_comparison_msstats_60)
## BOOM


## Run the differential analysis using MSstats
DA_msstats_60 <- MSstats::groupComparison(contrast.matrix = mat_comparison_msstats_60, data = spectronaut_proposed_60)

## save memory, compost
spectronaut_proposed_60 <- NULL
char_filename_date <- format(Sys.time(), "%Y_%m_%d_%H_%M_%S")
# save.image(paste0("./data_images/", char_filename_date, "DA_msstats_60.RData"))

## Pull out the comparison table
df_comp_60 <- DA_msstats_60$ComparisonResult

## save memory, compost
DA_msstats_60 <- NULL
char_filename_date <- format(Sys.time(), "%Y_%m_%d_%H_%M_%S")
save.image(paste0("./", char_filename_date, "df_comp_30_df_comp_60.RData"))
