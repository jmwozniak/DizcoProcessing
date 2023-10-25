#install.packages("scales")
#install.packages("vioplot")
library(scales)
library(ggplot2)
library(reshape2)
library(VennDiagram)
library(pROC)
library(mzR)
library(protViz)
library(OrgMassSpecR)
library(Epi)
library(mzID)
library(stringr)

extract_labeled <- function(uData, modIDs){
  
  if(is.null(nrow(uData)))
  {
    labeled_psms <- rep(1, length(uData))
    
    for(i in 1:length(uData)){
      labeled_psms[i] <- list(uData[[i]][which(grepl(paste(modIDs,collapse="|"),uData[[i]]$Modifications)),])
      rownames(labeled_psms[[i]]) <- seq(from = 1, to = nrow(labeled_psms[[i]]))
    }
  } else {
    labeled_psms <- uData[which(grepl(paste(modIDs,collapse="|"),uData$Modifications)),]
    rownames(labeled_psms) <- seq(from = 1, to = nrow(labeled_psms))
  }
  
  return(labeled_psms)
  
}
extract_proteins <- function(uData, uni = FALSE){
  
  if(is.null(nrow(uData)))
  {
    prot_ids <- rep("", length(uData))
    simple_prots <- NULL
    
    for(i in 1:length(uData)){
      if(uni == TRUE){
        simple_prots <- gsub(";.*","",uData[[i]]$Master.Protein.Accessions)
        prot_ids[i] <- list(unique(simple_prots))
      }
      else{
        simple_prots <- gsub(";.*","",uData[[i]]$Master.Protein.Accessions)
        prot_ids[i] <- list(simple_prots)
      }
    }
  } else {
    if(uni == TRUE){
      simple_prots <- gsub(";.*","",uData$Master.Protein.Accessions)
      prot_ids <- list(unique(simple_prots))
    }
    else{
      simple_prots <- gsub(";.*","",uData$Master.Protein.Accessions)
      prot_ids <- list(simple_prots)
    }
  }
  
  return(prot_ids)
  
}

processPSMs_isotope <- function(uData, modIDs, ptmRS=TRUE){
  #Remove unnecessary columns
  if(ptmRS == TRUE){
    uData <- uData[,c(2,4,5:6,8:13,16:20,25:29,33)]
  } else {
    uData <- uData[,c(2,4,5:6,8:13,16:20,25:29)]
  }
  
  
  #Shift "Protein.Accessions" to "Master.Protein.Accessions" if missing
  uData$Master.Protein.Accessions[which(uData$Master.Protein.Accessions=="")] <- uData$Protein.Accessions[which(uData$Master.Protein.Accessions=="")]
  
  #Keep only first "Master.Protein.Accessions"
  uData$Master.Protein.Accessions <- ifelse (regexpr(";", uData$Master.Protein.Accessions) != -1, 
                                             substr(uData$Master.Protein.Accessions, 1, regexpr(";", uData$Master.Protein.Accessions)-1), 
                                             uData$Master.Protein.Accessions)
  
  #Process peptide sequence and extract flanking AAs 
  split_seqs <- matrix(unlist(strsplit(uData[,"Annotated.Sequence"], "\\.")), ncol = 3, byrow = TRUE)
  
  uData[,"Upper_Seq"] <- toupper(split_seqs[,2])
  
  ######## Process "Modifications" column to remove Oxidation/Carbamidomethylation ########
  probe_sites <- gsub(";","",uData$Modifications)
  probe_sites <- gsub(" ","", probe_sites)
  probe_sites <- gsub("M.{1,2}[(]Oxidation[)]","", probe_sites)
  probe_sites <- gsub("C.{1,2}[(]Carbamidomethyl[)]","",probe_sites)
  probe_sites <- gsub("[(].*","",probe_sites)
  
  ######## Find out what residue is labeled (AA and location) ######## 
  uData[,"label_site"] <- probe_sites
  
  ######## Find out what amino acids are labeled ######## 
  uData[,"label_AA"] <- substr(probe_sites, 1, 1)
  
  ######## Find out what peptide locations are labeled ######## 
  uData[,"label_loc"] <- substr(probe_sites, 2, nchar(probe_sites))
  
  #Create unique ID with scan number and file number
  uData[,"scanID"] <- paste(uData[,"First.Scan"], uData[,"File.ID"], sep = "_")
  
  #Make new ID for scanID and peptide sequence
  uData[,"scanID_pep"] <- paste(uData[,"scanID"], toupper(uData[,"Annotated.Sequence"]), sep = "_")
  
  #Make unique ID with scanID_pep and Mods (will get around duplicates that just have different mod locations)
  uData[,"uniqueID"] <- paste(uData$scanID_pep, uData$Modifications, sep = "_")
  
  #Calculate peptide length
  uData[,"pepLength"] <- nchar(uData[,"Upper_Seq"])
  
  #Make new data frames with labeled peptides      
  labeled_psms <- extract_labeled(uData, modIDs)
  labeled_psms[,"labeled"] <- 1
  
  #Mark labeled peptdies in main uData data frame
  uData[,"labeled"] <- ifelse((uData[,"uniqueID"] %in% labeled_psms[,"uniqueID"]),1,0)
  
  #Make new data frames with unlabeled peptides     
  unlabeled_psms <- uData[uData$labeled==0,] 
  unlabeled_psms[,"paired"] <- 0
  
  light_mod <- modIDs[1]
  heavy_mod <- modIDs[2]
  
  #Make new data frames with only heavy or light psms
  light_psms <- labeled_psms[grep(light_mod, labeled_psms[,"Modifications"]),]
  heavy_psms <- labeled_psms[grep(heavy_mod, labeled_psms [,"Modifications"]),]
  
  # light_only <- unique(light_psms$Upper_Seq[ which(!(light_psms$Upper_Seq %in% heavy_psms$Upper_Seq)) ])
  # heavy_only <- unique(light_psms$Upper_Seq[ which(!(heavy_psms$Upper_Seq %in% light_psms$Upper_Seq)) ])
  # shared <- unique(light_psms$Upper_Seq[which(light_psms$Upper_Seq %in% heavy_psms$Upper_Seq)])
  
  # venn_lists <- list(unique(light_psms$Upper_Seq), unique(heavy_psms$Upper_Seq))
  # names(venn_lists) <- c("Light","Heavy")
  # venn.diagram(venn_lists, filename="H-Lpairing_VD")
  
  
  #Mark which label type is on a given peptide
  uData[,"label_type"] <- ifelse((uData[,"scanID"] %in% light_psms[,"scanID"]),"L",
                                 ifelse((uData[,"scanID"] %in% heavy_psms[,"scanID"]), "H", "NA"))
  
  
  #Determine if a peptide is paired based on its presence in both the light and heavy psms data frames
  labeled_psms[,"paired"] <- ifelse((labeled_psms[,"Upper_Seq"] %in% light_psms[,"Upper_Seq"])&(labeled_psms[,"Upper_Seq"] %in% heavy_psms[,"Upper_Seq"]),1,0)
  
  #Extract data frame "paired_psms" based on whether a peptide is paired or not
  paired_psms <- labeled_psms[which(labeled_psms$paired==1),]
  
  #Determine if a peptide is paired based on the "paired_psms" data frame
  uData[,"paired"] <- ifelse((uData[,"Annotated.Sequence"] %in% paired_psms[,"Annotated.Sequence"])&(uData[,"labeled"]==1),1,0)
  
  #Count # of PSMs for each scan ID
  uData[,"numPSMs"] <- table(uData$scanID)[uData$scanID]
  
  #Z score numPSMs for model building
  test_numPSMs_mean <- mean(uData$numPSMs)
  test_numPSMs_sd <- sd(uData$numPSMs)
  
  uData[,"numPSMs_scaled"] <- (uData$numPSMs - test_numPSMs_mean)/test_numPSMs_sd
  
  #Count % agree of PSMs for each scan ID
  uData[,"agreePSMs"] <- table(uData$scanID_pep)[uData$scanID_pep] / table(uData$scanID)[uData$scanID]
  
  #Reformat deltascore and deltacn into 1 column
  uData[,"scoreDiff"] <- ifelse(is.na(uData[,"Delta.Score"]),uData[,"Delta.Cn"],uData[,"Delta.Score"])
  
  #Z score scoreDiff for model building
  test_scoreDiff_mean <- mean(uData$scoreDiff)
  test_scoreDiff_sd <- sd(uData$scoreDiff)
  
  uData$scoreDiff_scaled <- (uData$scoreDiff - test_scoreDiff_mean)/test_scoreDiff_sd
  
  #Create a matrix for average RTs for labeled or unlabeled peptides of any length 
  RT_mat <- matrix(data = NA, nrow = length(unique(uData$pepLength)), ncol= 3)
  RT_mat[,1] <- unique(uData$pepLength)
  for(i in 1:length(unique(uData$pepLength))){
    RT_mat[i,2] <- mean(unlabeled_psms$RT.in.min[(unlabeled_psms$pepLength==unique(uData$pepLength)[i])])
    RT_mat[i,3] <- mean(paired_psms$RT.in.min[(paired_psms$pepLength==unique(uData$pepLength)[i])])
  }
  RT_mat <- RT_mat[order(RT_mat[,1]),]
  
  #Rename RT matrix columns accordingly
  colnames(RT_mat)<-c("pepLength", "RTavg_fromUPL", "RTavg_fromLPL")
  
  #Add RT information to main uData data frame
  uData <- merge(uData, RT_mat, by = 'pepLength')
  
  #Calculate RT Differnece from unlabeled peptides of same length
  uData[,"RT_Diff_fromUPL"] <- uData$RTavg_fromUPL - uData$RT.in.min
  
  #Z score RT Difference for model building
  test_RT_mean <- mean(uData$RT_Diff_fromUPL[!is.na(uData$RT_Diff_fromUPL)])
  test_RT_sd <- sd(uData$RT_Diff_fromUPL[!is.na(uData$RT_Diff_fromUPL)])
  
  uData[,"RT_Diff_fromUPL_scaled"] <- (uData$RT_Diff_fromUPL- test_RT_mean)/test_RT_sd
  
  #Make new ID for Master protein and peptide sequence
  uData[,"ProtID_Pep"] <- paste(uData[,"Master.Protein.Accessions"], toupper(uData[,"Upper_Seq"]), sep = "_")
  
  #Count total PSMs per protein
  uData <- merge(uData, table(uData$Master.Protein.Accessions), by.x = "Master.Protein.Accessions", by.y = "Var1")
  colnames(uData)[ncol(uData)] <- "Prot_totalPSMs"
  
  #Count labeled PSMs per protein
  uData <- merge(uData, table(labeled_psms$Master.Protein.Accessions), by.x = "Master.Protein.Accessions", by.y = "Var1", all = TRUE)
  colnames(uData)[ncol(uData)] <- "Prot_labeledPSMs"
  
  #Count paired PSMs per protein
  uData <- merge(uData, table(paired_psms$Master.Protein.Accessions), by.x = "Master.Protein.Accessions", by.y = "Var1", all = TRUE)
  colnames(uData)[ncol(uData)] <- "Prot_pairedPSMs"
  
  #Count unique peptides per protein
  uData <- merge(uData, table(unique(uData[,c("Master.Protein.Accessions","Upper_Seq")])$Master.Protein.Accessions), by.x = "Master.Protein.Accessions", by.y = "Var1", all = TRUE)
  colnames(uData)[ncol(uData)] <- "Prot_uniquePeptides"
  
  #Count labeled unique peptides per protein
  uData <- merge(uData, table(unique(labeled_psms[,c("Master.Protein.Accessions","Upper_Seq")])$Master.Protein.Accessions), by.x = "Master.Protein.Accessions", by.y = "Var1", all = TRUE)
  colnames(uData)[ncol(uData)] <- "Prot_labeledPeptides"
  
  #Count labeled paired peptides per protein
  uData <- merge(uData, table(unique(paired_psms[,c("Master.Protein.Accessions","Upper_Seq")])$Master.Protein.Accessions), by.x = "Master.Protein.Accessions", by.y = "Var1", all = TRUE)
  colnames(uData)[ncol(uData)] <- "Prot_pairedPeptides"
  
  #Split unlabeled psms in half (by every other row) to make unique negative sets for models
  trainNegative <- unlabeled_psms[seq(1, nrow(unlabeled_psms), 2), ]
  testNegative <- unlabeled_psms[seq(2, nrow(unlabeled_psms), 2), ]
  
  trainSet <- rbind(uData[which(uData$paired==1), ], uData[which(uData$uniqueID %in% trainNegative$uniqueID),])
  
  testSet <- rbind(uData[intersect(which((uData$paired==0)),which((uData$labeled==1))), ], uData[which(uData$uniqueID %in% testNegative$uniqueID),])
  
  glm.fit <- glm(labeled ~ numPSMs_scaled + scoreDiff_scaled + agreePSMs + RT_Diff_fromUPL_scaled,
                 data = trainSet,
                 family = binomial)
  
  glm.probs_train <- predict(glm.fit,
                             newdata = trainSet,
                             type = "response")
  
  glm.probs_test <- predict(glm.fit,
                            newdata = testSet,
                            type = "response")
  
  glm.pred_train = ifelse(glm.probs_train > 0.5, "Labeled", "Unlabeled")
  
  train_out <- cbind(glm.probs_train, glm.pred_train, trainSet)
  colnames(train_out)[1:2] <- c("glm.probs", "glm.pred")
  train_out["model_set"] <- "train"
  #write.table(train_out, file="train_preds.txt", sep="\t")
  
  glm.pred_test = ifelse(glm.probs_test > 0.5, "Labeled", "Unlabeled")
  
  test_out <- cbind(glm.probs_test, glm.pred_test, testSet)
  colnames(test_out)[1:2] <- c("glm.probs", "glm.pred")
  test_out["model_set"] <- "test"
  #write.table(test_out, file="test_preds.txt",sep="\t")
  
  merge_out <- rbind(train_out, test_out)
  
  merge_out[,"label_confidence"] <- ifelse(merge_out[,"glm.probs"] >= 0.85,"High",
                                           ifelse(merge_out[,"glm.probs"] >= 0.50, "Medium", "Low"))
  
  return(merge_out)
}
processPeps_isotope <- function(processedPSMs){
  #Create unique scan IDs from spectrum file and scan #
  processedPSMs$unique_scanID <- paste(processedPSMs$First.Scan,processedPSMs$Spectrum.File,sep="_")
  
  #Extract label probe information from spectrum file (and only keep info if peptide is labeled)
  processedPSMs$probeID <- substring(processedPSMs$Spectrum.File, regexpr("SOL_",processedPSMs$Spectrum.File)+4, regexpr("_R",processedPSMs$Spectrum.File)-1)
  processedPSMs$probeID[which(processedPSMs$labeled==0)] <- ""
  
  #Create unique ID with unique peptide sequence and label status
  processedPSMs$pep_label <- paste(processedPSMs$Upper_Seq, processedPSMs$probeID, sep = "_")
  
  
  #Filter Data for Unique Scans -> keep peptide matched to protein with most PSMs and keep PSM with highest XCorr
  processedPSMs_sorted <- processedPSMs[with(processedPSMs, order(unique_scanID, Prot_totalPSMs, XCorr)),]
  #write.table(processedPSMs_sorted, file="processedPSMs_sorted_Check.txt", sep = "\t", quote = FALSE, row.names = FALSE)
  
  uniqueScans <- unique(processedPSMs$unique_scanID)
  
  #write.table(uniqueScans, file="uniqueScan_Check.txt", sep = "\t", quote = FALSE, row.names = FALSE)
  
  uniqueScan_DF <- data.frame(matrix(nrow=length(uniqueScans), ncol = ncol(processedPSMs_sorted)))
  
  for(i in 1:length(uniqueScans)){
    uniqueScan_DF[i,] <- processedPSMs_sorted[grep(paste(c("\\b", uniqueScans[i], "\\b"), collapse=""), processedPSMs_sorted$unique_scanID)[1],]
    print(i)
  }
  
  colnames(uniqueScan_DF) <- colnames(processedPSMs_sorted) 
  write.table(uniqueScan_DF, file = paste(udata_path, "uniqueScanDF_LFQmerged_", Sys.Date(), ".txt", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE)
  
  #Read in uniqueScan_DF
  #uniqueScan_DF <- read.table("uniqueScanDF_LFQmerged_2021-11-17.txt", header=TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  #Extract label probe information from spectrum file (and only keep info if peptide is labeled)
  uniqueScan_DF$probeID <- substring(uniqueScan_DF$Spectrum.File, regexpr("SOL_",uniqueScan_DF$Spectrum.File)+4, regexpr("_R",uniqueScan_DF$Spectrum.File)-1)
  uniqueScan_DF$probeID[which(uniqueScan_DF$labeled==0)] <- ""
  
  #Create unique ID with unique peptide sequence and label status
  uniqueScan_DF$pep_label <- paste(uniqueScan_DF$Upper_Seq, uniqueScan_DF$probeID, sep = "_")
  
  #Create new DF using unique peptides/FFF label status
  pepDF <- data.frame(unique(uniqueScan_DF$pep_label))
  colnames(pepDF)[1] <- "Unique_ID"
  
  #Extract protein IDs from each unique peptide/FFF label status
  # for MS2 pepDF <- unique(merge(pepDF, uniqueScan_DF[,c(3,46)], by.x = "Unique_ID", by.y = "pep_label", all.x = TRUE))
  pepDF <- unique(merge(pepDF, uniqueScan_DF[,c("Master.Protein.Accessions","pep_label")], by.x = "Unique_ID", by.y = "pep_label", all.x = TRUE))
  
  colnames(pepDF)[ncol(pepDF)] <- "Uniprot_ID"
  
  #Count the number of unique scans for each peptide
  pepDF <- merge(pepDF, table(uniqueScan_DF$pep_label), by.x = "Unique_ID", by.y = "Var1", all.x = TRUE)
  colnames(pepDF)[ncol(pepDF)] <- "#_Spectra"
  
  #Count the number of PSMs for each peptide
  pepDF <- merge(pepDF, table(processedPSMs$pep_label), by.x = "Unique_ID", by.y = "Var1", all.x = TRUE)
  colnames(pepDF)[ncol(pepDF)] <- "#_PSMs"
  
  ####### Average Chimeric Spectra Features per peptide #######
  
  # avgXCorr <- aggregate(processedPSMs$XCorr, by = list(Unique_ID = processedPSMs$pep_label), FUN = mean)
  # avgDscore <- aggregate(processedPSMs$scoreDiff, by = list(Unique_ID = processedPSMs$pep_label), FUN = mean)
  # avgNumPSMs <- aggregate(processedPSMs$numPSMs, by = list(Unique_ID = processedPSMs$pep_label), FUN = mean)
  # avgRTshift <- aggregate(processedPSMs$RT_Diff_fromUPL, by = list(Unique_ID = processedPSMs$pep_label), FUN = mean)
  
  featureCols <- c("glm.probs", "XCorr", "scoreDiff", "numPSMs", "RT_Diff_fromUPL")
  
  averagedFeatures <- aggregate(processedPSMs[,featureCols], by = list(Unique_ID = processedPSMs$pep_label), FUN = mean)
  pepDF <- merge(pepDF, averagedFeatures, by = "Unique_ID")
  
  max_glm <- aggregate(processedPSMs[,"glm.probs"], by = list(Unique_ID = processedPSMs$pep_label), FUN = max)
  
  pepDF <- merge(pepDF, max_glm, by = "Unique_ID")
  
  colnames(pepDF)[ncol(pepDF)] <- "max.glm.probs"
  
  ####### Aggregate Label Localization Info per peptide #######
  for(i in 1:nrow(pepDF)){
    all_sites <- processedPSMs$label_site[which(processedPSMs$pep_label == pepDF$Unique_ID[i])]
    
    sites <- unique(processedPSMs$label_site[which(processedPSMs$pep_label == pepDF$Unique_ID[i])])
    locs <- unique(processedPSMs$label_loc[which(processedPSMs$pep_label == pepDF$Unique_ID[i])])
    aas <- substr(sites, start = 1, stop = 1)
    
    siteDF <- data.frame(sites)
    siteDF <- cbind(siteDF, aas, locs)
    siteDF <- merge(siteDF, table(all_sites), by.x = "sites", by.y = "all_sites")
    
    siteDF <- siteDF[order(siteDF$locs),]
    
    pepDF$sites[i] <- paste(siteDF$sites, collapse = ";")
    pepDF$aas[i] <- paste(siteDF$aas, collapse = ";")
    pepDF$locs[i] <- paste(siteDF$locs, collapse = ";")
    pepDF$site_psms[i] <- paste(siteDF$Freq, collapse = ";")
    
    print(i)
  }
  
  pepDF$PEPTIDE.SEQUENCE <- substring(pepDF$Unique_ID, 1, regexpr("_",pepDF$Unique_ID)-1)
  
  return(pepDF)
}

processPSMs_TMT <- function(uData, modIDs, train_set_forPred, ptmRS = TRUE, plex = "16plex", quant_level = "MS3"){
  
    #############################################################################################################################
    #########################               TMT Text processing and predictive modeling                     #########################
    #############################################################################################################################
    
    #Only retain necessary columns
    uData_cols <- c("Confidence", "PSM.Ambiguity", "Annotated.Sequence", "Modifications", "Master.Protein.Accessions", "Protein.Accessions", "Number.of.Missed.Cleavages", "Charge", 
                    "Delta.Score", "Delta.Cn", "mz.in.Da", "MHplus.in.Da", "Theo.MHplus.in.Da", "Delta.M.in.ppm", "Delta.mz.in.Da", "Isolation.Interference.in.Percent", 
                    "Average.Reporter.SN", "Ion.Inject.Time.in.ms", "RT.in.min", "First.Scan", "Spectrum.File", "File.ID", "Quan.Info", "XCorr")  
    
    if(plex == "10plex"){
      uData_cols <- c(uData_cols, "Abundance.126", "Abundance.127N", "Abundance.127C", "Abundance.128N", "Abundance.128C", 
                      "Abundance.129N", "Abundance.129C", "Abundance.130N", "Abundance.130C", "Abundance.131") 
      
      #Process "Modifications" column to remove Oxidation/Carbamidomethylation
      probe_sites <- gsub(";","",uData$Modifications)
      probe_sites <- gsub(" ","", probe_sites)
      probe_sites <- gsub("M.{1,2}[(]Oxidation[)]","", probe_sites)
      probe_sites <- gsub("C.{1,2}[(]Carbamidomethyl[)]","",probe_sites)
      probe_sites <- gsub("N-Term[(]TMT6plex[)]","", probe_sites)
      probe_sites <- gsub("K.{1,2}[(]TMT6plex[)]","",probe_sites)
      probe_sites <- gsub("[(].*","",probe_sites)
    }
    
    if(plex == "16plex") {
      uData_cols <- c(uData_cols, "Abundance.126", "Abundance.127N", "Abundance.127C", "Abundance.128N", "Abundance.128C", 
                      "Abundance.129N", "Abundance.129C", "Abundance.130N", "Abundance.130C", "Abundance.131N", 
                      "Abundance.131C", "Abundance.132N", "Abundance.132C", "Abundance.133N", "Abundance.133C",
                      "Abundance.134N") 
      
      #Process "Modifications" column to remove Oxidation/Carbamidomethylation
      probe_sites <- gsub(";","",uData$Modifications)
      probe_sites <- gsub(" ","", probe_sites)
      probe_sites <- gsub("M.{1,2}[(]Oxidation[)]","", probe_sites)
      probe_sites <- gsub("C.{1,2}[(]Carbamidomethyl[)]","",probe_sites)
      probe_sites <- gsub("N-Term[(]TMTpro[)]","", probe_sites)
      probe_sites <- gsub("K.{1,2}[(]TMTpro[)]","",probe_sites)
      probe_sites <- gsub("[(].*","",probe_sites)
    }
    
    if(quant_level == "MS3"){
      uData_cols <- c(uData_cols, "SPS.Mass.Matches.in.Percent")
    }
    
    if(ptmRS == TRUE){ 
      uData_cols <- c(uData_cols, "ptmRS.Best.Site.Probabilities")
    }
    
    uData <- uData[,uData_cols]
    
    #Shift "Protein.Accessions" to "Master.Protein.Accessions" if missing
    uData$Master.Protein.Accessions[which(uData$Master.Protein.Accessions=="")] <- uData$Protein.Accessions[which(uData$Master.Protein.Accessions=="")]
    
    #Keep only first "Master.Protein.Accessions"
    uData$Master.Protein.Accessions <- ifelse (regexpr(";", uData$Master.Protein.Accessions) != -1, 
                                               substr(uData$Master.Protein.Accessions, 1, regexpr(";", uData$Master.Protein.Accessions)-1), 
                                               uData$Master.Protein.Accessions)
    
    #Process peptide sequence and extract flanking AAs 
    split_seqs <- matrix(unlist(strsplit(uData[,"Annotated.Sequence"], "\\.")), ncol = 3, byrow = TRUE)
    
    uData[,"Upper_Seq"] <- toupper(split_seqs[,2])
    
    #Find out what residue is labeled (AA and location)
    uData[,"label_site"] <- probe_sites
    
    #Find out what amino acids are labeled
    uData[,"label_AA"] <- substr(probe_sites, 1, 1)
    
    #Find out what peptide locations are labeled
    uData[,"label_loc"] <- substr(probe_sites, 2, nchar(probe_sites))
    
    #Find out which label is on each peptide
    for(i in 1:length(modIDs)){
      uData[,modIDs[i]] <- as.integer(grepl(modIDs[i], uData$Modifications))
    }
    
    #Create sum ID with unique peptide sequence and label status
    cols_to_paste <- c(grep(paste(c("Upper_Seq", modIDs), collapse = "|"), colnames(uData)))
    
    uData[,"pep_label"] <- do.call(paste, c(uData[,cols_to_paste], sep = "_"))
    
    #Create unique ID with scan number and file number
    uData[,"scanID"] <- paste(uData[,"First.Scan"], uData[,"File.ID"], sep = "_")
    
    #Make new ID for scanID and peptide sequence
    uData[,"scanID_pep"] <- paste(uData[,"scanID"], toupper(uData[,"Annotated.Sequence"]), sep = "_")
    
    #Make unique ID with scanID_pep and Mods (will get around duplicates that just have different mod locations)
    uData[,"uniqueID"] <- paste(uData$scanID_pep, uData$Modifications, sep = "_")
    
    #Calculate peptide length
    uData[,"pepLength"] <- nchar(uData[,"Upper_Seq"])
    
    #Count # TMT labels
    uData[,"numTMT"] <- str_count(uData$Modifications, "TMT")
    
    #Combine pepLength and TMT labels
    uData[,"pepLength_numTMT"] <- paste(uData$pepLength, uData$numTMT, sep = "_")
    
    uData <- uData[!grepl("CP266", uData$Modifications),]
    #Make new data frames with labeled peptides      
    labeled_psms <- extract_labeled(uData, modIDs)
    labeled_psms[,"labeled"] <- 1
    
    #Mark labeled peptides in main uData data frame
    uData[,"labeled"] <- ifelse((uData[,"uniqueID"] %in% labeled_psms[,"uniqueID"]),1,0)
    
    #Make new data frames with unlabeled peptides     
    unlabeled_psms <- uData[uData$labeled==0,]
    
    #Count # of PSMs for each scan ID
    uData[,"numPSMs"] <- table(uData$scanID)[uData$scanID]
    
    test_numPSMs_mean <- mean(uData$numPSMs)
    test_numPSMs_sd <- sd(uData$numPSMs)
    
    uData[,"numPSMs_scaled"] <- (uData$numPSMs - test_numPSMs_mean)/test_numPSMs_sd
    
    #Count % agree of PSMs for each scan ID
    uData[,"agreePSMs"] <- table(uData$scanID_pep)[uData$scanID_pep] / table(uData$scanID)[uData$scanID]
    
    #Reformat deltascore and deltacn into 1 column
    uData[,"scoreDiff"] <- ifelse(is.na(uData[,"Delta.Score"]),uData[,"Delta.Cn"],uData[,"Delta.Score"])
    
    test_scoreDiff_mean <- mean(uData$scoreDiff)
    test_scoreDiff_sd <- sd(uData$scoreDiff)
    
    uData$scoreDiff_scaled <- (uData$scoreDiff - test_scoreDiff_mean)/test_scoreDiff_sd
    
    #Create a matrix for average RTs for labeled or unlabeled peptides of any length with any number of TMT labels (1-4)
    RT_mat <- matrix(data = NA, nrow = length(unique(uData$pepLength_numTMT)), ncol= 3)
    RT_mat[,1] <- unique(uData$pepLength_numTMT)
    for(i in 1:length(unique(uData$pepLength_numTMT))){
      RT_mat[i,2] <- mean(unlabeled_psms$RT.in.min[(unlabeled_psms$pepLength_numTMT==unique(uData$pepLength_numTMT)[i])])
      RT_mat[i,3] <- mean(labeled_psms$RT.in.min[(labeled_psms$pepLength_numTMT==unique(uData$pepLength_numTMT)[i])])
    }
    RT_mat <- RT_mat[order(RT_mat[,1]),]
    
    #Rename RT matrix columns accordingly
    colnames(RT_mat)<-c("pepLength_numTMT", "RTavg_fromUPL", "RTavg_fromLPL")
    
    #Add RT information to main uData data frame
    uData <- merge(uData, RT_mat, by = 'pepLength_numTMT')
    
    #Calculate RT Differnece from unlabeled peptides of same length
    uData[,"RT_Diff_fromUPL"] <- as.numeric(uData$RTavg_fromUPL) - as.numeric(uData$RT.in.min)
    
    test_RT_mean <- mean(uData$RT_Diff_fromUPL[!is.na(uData$RT_Diff_fromUPL)])
    test_RT_sd <- sd(uData$RT_Diff_fromUPL[!is.na(uData$RT_Diff_fromUPL)])
    
    uData[,"RT_Diff_fromUPL_scaled"] <- (uData$RT_Diff_fromUPL- test_RT_mean)/test_RT_sd
    
    #Make new ID for Master protein and peptide sequence
    uData[,"ProtID_Pep"] <- paste(uData[,"Master.Protein.Accessions"], toupper(uData[,"Upper_Seq"]), sep = "_")
    
    #Count total PSMs per protein
    uData <- merge(uData, table(uData$Master.Protein.Accessions), by.x = "Master.Protein.Accessions", by.y = "Var1")
    colnames(uData)[ncol(uData)] <- "Prot_totalPSMs"
    
    #Count labeled PSMs per protein
    uData <- merge(uData, table(labeled_psms$Master.Protein.Accessions), by.x = "Master.Protein.Accessions", by.y = "Var1", all = TRUE)
    colnames(uData)[ncol(uData)] <- "Prot_labeledPSMs"
    
    #Count unique peptides per protein
    uData <- merge(uData, table(unique(uData[,c("Master.Protein.Accessions","Upper_Seq")])$Master.Protein.Accessions), by.x = "Master.Protein.Accessions", by.y = "Var1", all = TRUE)
    colnames(uData)[ncol(uData)] <- "Prot_uniquePeptides"
    
    #Count labeled unique peptides per protein
    uData <- merge(uData, table(unique(labeled_psms[,c("Master.Protein.Accessions","Upper_Seq")])$Master.Protein.Accessions), by.x = "Master.Protein.Accessions", by.y = "Var1", all = TRUE)
    colnames(uData)[ncol(uData)] <- "Prot_labeledPeptides"
    
    train_RT_mean <- mean(train_set_forPred$RT_Diff_fromUPL)
    train_RT_sd <- sd(train_set_forPred$RT_Diff_fromUPL)
    
    train_set_forPred$RT_Diff_fromUPL_scaled <- (train_set_forPred$RT_Diff_fromUPL - train_RT_mean)/train_RT_sd
    
    train_scoreDiff_mean <- mean(train_set_forPred$scoreDiff)
    train_scoreDiff_sd <- sd(train_set_forPred$scoreDiff)
    
    train_set_forPred$scoreDiff_scaled <- (train_set_forPred$scoreDiff - train_scoreDiff_mean)/train_scoreDiff_sd
    
    train_numPSMs_mean <- mean(train_set_forPred$numPSMs)
    train_numPSMs_sd <- sd(train_set_forPred$numPSMs)
    
    train_set_forPred$numPSMs_scaled <- (train_set_forPred$numPSMs - train_numPSMs_mean)/train_numPSMs_sd
    
    #Read in previous H/L paired data to train model
    glm.fit <- glm(labeled ~ numPSMs + scoreDiff + RT_Diff_fromUPL_scaled + agreePSMs,
                   data = train_set_forPred,
                   family = binomial)
    
    glm.probs_train <- predict(glm.fit,
                               newdata = train_set_forPred,
                               type = "response")
    
    glm.probs_test <- predict(glm.fit,
                              newdata = uData,
                              type = "response")
    
    glm.pred_train = ifelse(glm.probs_train > 0.5, "Labeled", "Unlabeled")
    
    train_out <- cbind(glm.probs_train, glm.pred_train, train_set_forPred)
    colnames(train_out)[1:2] <- c("glm.probs", "glm.pred")
    train_out["model_set"] <- "train"
    #write.table(train_out, file="train_preds.txt", sep="\t")
    
    table(train_out$glm.pred, train_out$labeled)

    glm.pred_test = ifelse(glm.probs_test > 0.5, "Labeled", "Unlabeled")
    
    test_out <- cbind(glm.probs_test, glm.pred_test, uData)
    colnames(test_out)[1:2] <- c("glm.probs", "glm.pred")
    test_out["model_set"] <- "test"
    #write.table(test_out, file="test_preds.txt",sep="\t")
    
    table(test_out$glm.pred, test_out$labeled)
    
    test_out[,"label_confidence"] <- ifelse(test_out[,"glm.probs"] >= 0.85,"High",
                                             ifelse(test_out[,"glm.probs"] >= 0.50, "Medium", "Low"))
    
    return(test_out)
}
processPeps_TMT <- function(processedPSMs, modIDs, plex = "16plex"){
  
  processedPSMs_sorted <- processedPSMs[with(processedPSMs, order(scanID, Prot_totalPSMs, XCorr)),]
  #write.table(processedPSMs_sorted, file="processedPSMs_sorted_Check.txt", sep = "\t", quote = FALSE, row.names = FALSE)
  
  uniqueScans <- unique(processedPSMs$scanID)
  #write.table(uniqueScans, file="uniqueScan_Check.txt", sep = "\t", quote = FALSE, row.names = FALSE)
  
  uniqueScan_DF <- data.frame(matrix(nrow=length(uniqueScans), ncol = ncol(processedPSMs_sorted)))
  
  for(i in 1:length(uniqueScans)){
    uniqueScan_DF[i,] <- processedPSMs_sorted[grep(paste(c("\\b", uniqueScans[i], "\\b"), collapse=""), processedPSMs_sorted$scanID)[1],]
    print(i)
  }
  
  colnames(uniqueScan_DF) <- colnames(processedPSMs_sorted) 
  write.csv(uniqueScan_DF, file = paste(udata_path, "uniqueScansDF_", Sys.Date(), ".csv", sep = ""), quote = FALSE, row.names = FALSE)
  
  #Read in uniqueScan_DF
  #uniqueScan_DF <- read.table("uniqueScanDF_2021-11-05.txt", header=TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  #Create new DF using unique peptides/FFF label status
  pepDF <- data.frame(unique(uniqueScan_DF$pep_label))
  colnames(pepDF)[1] <- "Unique_ID"
  
  #Extract protein IDs from each unique peptide/FFF label status
  pepDF <- unique(merge(pepDF, uniqueScan_DF[,c("Master.Protein.Accessions","pep_label")], by.x = "Unique_ID", by.y = "pep_label", all.x = TRUE))
  
  colnames(pepDF)[ncol(pepDF)] <- "Uniprot_ID"
  
  #Count the number of unique scans for each peptide
  pepDF <- merge(pepDF, table(uniqueScan_DF$pep_label), by.x = "Unique_ID", by.y = "Var1", all.x = TRUE)
  colnames(pepDF)[ncol(pepDF)] <- "#_Spectra"
  
  #Count the number of PSMs for each peptide
  pepDF <- merge(pepDF, table(processedPSMs$pep_label), by.x = "Unique_ID", by.y = "Var1", all.x = TRUE)
  colnames(pepDF)[ncol(pepDF)] <- "#_PSMs"
  
  ####### Sum TMT Data per peptide #######
  
  if(plex == "10plex"){
    tmtData <- uniqueScan_DF[,c("Abundance.126", "Abundance.127N", "Abundance.127C", "Abundance.128N", "Abundance.128C", 
                                "Abundance.129N", "Abundance.129C", "Abundance.130N", "Abundance.130C", "Abundance.131", "pep_label")]
    
    tmtData[is.na(tmtData)] <- 1
    #write.table(tmtData, "presummed_TMT.txt", quote = FALSE, sep ="\t", row.names = FALSE)
    
    summedTMT <- aggregate(tmtData[,c(1:10)], by = list(Unique_ID = tmtData$pep_label), FUN = sum)
    #write.table(summedTMT, "summed_TMT.txt", quote = FALSE, sep ="\t", row.names = FALSE)
    
  }
  
  if(plex == "16plex"){
    tmtData <- uniqueScan_DF[,c("Abundance.126", "Abundance.127N", "Abundance.127C", "Abundance.128N", "Abundance.128C", 
                                "Abundance.129N", "Abundance.129C", "Abundance.130N", "Abundance.130C", "Abundance.131N", 
                                "Abundance.131C", "Abundance.132N", "Abundance.132C", "Abundance.133N", "Abundance.133C",
                                "Abundance.134N", "pep_label")]
    
    tmtData[is.na(tmtData)] <- 1
    #write.table(tmtData, "presummed_TMT.txt", quote = FALSE, sep ="\t", row.names = FALSE)
    
    summedTMT <- aggregate(tmtData[,c(1:16)], by = list(Unique_ID = tmtData$pep_label), FUN = sum)
    #write.table(summedTMT, "summed_TMT.txt", quote = FALSE, sep ="\t", row.names = FALSE)
  }

  
  pepDF <- merge(pepDF, summedTMT, by = "Unique_ID")
  
  ####### Average Chimeric Spectra Features per peptide #######
  
  # avgXCorr <- aggregate(processedPSMs$XCorr, by = list(Unique_ID = processedPSMs$pep_label), FUN = mean)
  # avgDscore <- aggregate(processedPSMs$scoreDiff, by = list(Unique_ID = processedPSMs$pep_label), FUN = mean)
  # avgNumPSMs <- aggregate(processedPSMs$numPSMs, by = list(Unique_ID = processedPSMs$pep_label), FUN = mean)
  # avgRTshift <- aggregate(processedPSMs$RT_Diff_fromUPL, by = list(Unique_ID = processedPSMs$pep_label), FUN = mean)
  
  feats_toAvg <- c("glm.probs", "XCorr", "scoreDiff", "numPSMs", "RT_Diff_fromUPL")

  averagedFeatures <- aggregate(processedPSMs[,feats_toAvg], by = list(Unique_ID = processedPSMs$pep_label), FUN = mean)
  
  # for MS2 averagedFeatures <- aggregate(processedPSMs[,c(1,37,53,56,60)], by = list(Unique_ID = processedPSMs$pep_label), FUN = mean)
  pepDF <- merge(pepDF, averagedFeatures, by = "Unique_ID")
  
  max_glm <- aggregate(processedPSMs[,"glm.probs"], by = list(Unique_ID = processedPSMs$pep_label), FUN = max)
  
  pepDF <- merge(pepDF, max_glm, by = "Unique_ID")
  
  colnames(pepDF)[ncol(pepDF)] <- "max.glm.probs"
  
  ####### Aggregate Label Localization Info per peptide #######
  for(i in 1:nrow(pepDF)){
    all_sites <- processedPSMs$label_site[which(processedPSMs$pep_label == pepDF$Unique_ID[i])]
    
    sites <- unique(processedPSMs$label_site[which(processedPSMs$pep_label == pepDF$Unique_ID[i])])
    locs <- unique(processedPSMs$label_loc[which(processedPSMs$pep_label == pepDF$Unique_ID[i])])
    aas <- substr(sites, start = 1, stop = 1)
    
    siteDF <- data.frame(sites)
    siteDF <- cbind(siteDF, aas, locs)
    siteDF <- merge(siteDF, table(all_sites), by.x = "sites", by.y = "all_sites")
    
    siteDF <- siteDF[order(siteDF$locs),]
    
    pepDF$sites[i] <- paste(siteDF$sites, collapse = ";")
    pepDF$aas[i] <- paste(siteDF$aas, collapse = ";")
    pepDF$locs[i] <- paste(siteDF$locs, collapse = ";")
    pepDF$site_psms[i] <- paste(siteDF$Freq, collapse = ";")
    
    print(i)
  }
  
  temp <- matrix(unlist(strsplit(pepDF$Unique_ID, "_")), nrow = nrow(pepDF), ncol = length(modIDs)+1, byrow =TRUE)
  colnames(temp) <- c("PEPTIDE.SEQUENCE", modIDs)
  pepDF <- cbind(pepDF, temp)
  
  return(pepDF)
}

matchUniFeatures <- function(processedPeps, featureTab){
  #Change column names
  colnames(featureTab)[which(colnames(featureTab)=="Cross.reference..PDB.")] <- "PDB"
  colnames(featureTab)[which(colnames(featureTab)=="Cross.reference..ChEMBL.")] <- "ChEMBL"
  colnames(featureTab)[which(colnames(featureTab)=="Cross.reference..DrugBank.")] <- "DrugBank"
  
  #Create "Gene" column (first gene in Gene.names column)
  featureTab["Gene"] <- sub(" .*", "", featureTab$Gene.names)
  
  #Replace gene with Uniprot ID if no gene name available
  featureTab[which(featureTab$Gene == ""),]$Gene <- featureTab[which(featureTab$Gene == ""),]$Entry
  
  #Extract features for reviewed proteins and proteins from proteomic data
  featureTab_reviewed <- featureTab[which(featureTab$Status == "reviewed"),]
  featureTab_proteomics <- featureTab[which(featureTab$Entry %in% processedPeps$Uniprot_ID),]
  
  #Match gene name to proteomic data based on Uniprot ID
  proteomic_feature_cols <- c("Entry", "Gene")
  processedPeps <- merge(processedPeps, featureTab_proteomics[,proteomic_feature_cols], by.x="Uniprot_ID", by.y="Entry", all.x=TRUE)
  
  #Match reviewed features to proteomic data based on gene name
  #Note: if no reviewed gene (or no gene name), columns will be blank
  reviewed_feature_cols <- c("Protein.names", "Gene.names", "Active.site", "Binding.site", "Calcium.binding", "Cofactor", "DNA.binding", "Metal.binding", "Nucleotide.binding", "Site", "PDB", "ChEMBL", "DrugBank", "Sequence", "Gene", "AlphaFold", "Variants_LP")              
  processedPeps <- merge(processedPeps, featureTab_reviewed[,reviewed_feature_cols], by.x = "Gene", by.y = "Gene", all.x=TRUE)
  
  processedPeps$Labeled.Peptide <- ""
  processedPeps$TARG_PEPRANGE <- ""
  for(i in 1:nrow(processedPeps)){
    
    start <- regexpr(processedPeps$PEPTIDE.SEQUENCE[i], processedPeps$Sequence[i], ignore.case = TRUE)[[1]]
    end <- start + nchar(processedPeps$PEPTIDE.SEQUENCE[i]) - 1
    peprange <- paste(start, end, sep = "_")
    
    targ_peprange <- paste(processedPeps$Uniprot_ID[i], peprange, sep = "_")
    
    processedPeps$Labeled.Peptide[i] <- peprange
    processedPeps$TARG_PEPRANGE[i] <- targ_peprange
    print(i)
  }
  
  processedPeps[is.na(processedPeps)] <- ""
  
  return(processedPeps)
}

#####################################################################################
#######                                                                       #######
#######                      PSM Text File Processing                         #######
#######                                                                       #######
#####################################################################################

#Set directory to experimental data location
udata_path <- file.path("insert_file_path_here")

#Read in experimental data (.txt export from Proteome Discoverer)
uData <- read.table(paste(udata_path, "insert_file_name_here", sep =""), header=TRUE, sep = "\t", stringsAsFactors = FALSE)

#Set path to metadata (training model data, Uniprot data, AF data, etc.)
sol_metadata_path <- file.path("insert_file_path_here")

#Read in model training data and select H/L paired data for model generation
train_set_input <- read.table(paste(sol_metadata_path, "SoL-7-8_merged_PSMs_040722.txt", sep = ""), header = TRUE, row.names = 1, sep = "\t")
train_set_forPred_all <- train_set_input[which((train_set_input$paired==1)|(train_set_input$labeled==0)),]
train_set_forPred <- train_set_input[which((train_set_input$paired==1)|(train_set_input$labeled==0)),]

#Set Mod IDs (for identifying modified PSMs) and plex size
modIDs <- c("RS1")
#modIDs <- c("-AC-light", "-AC-heavy")

plex_size <- "16plex"

#Process SoL data 
processedPSMs <- processPSMs_TMT(uData, modIDs, train_set_forPred_all, ptmRS = TRUE, plex_size, quant_level <- "MS3")
#processedPSMs <- processPSMs_isotope(uData, modIDs, ptmRS = FALSE)

#Check model prediction (confusion matrix)
table(processedPSMs$glm.pred, processedPSMs$labeled)

#Write processedPSMs to a file
write.csv(processedPSMs, file = paste(udata_path, "insert_file_name_here", Sys.Date(), ".csv", sep = ""), quote = FALSE, row.names = FALSE)

############################## Process SoL Peptides (warning: can take a bit) ##############################

#Filter for labeled peptides
processedPSMs_filtered <- processedPSMs[which(processedPSMs$labeled==1),]

#Process SoL into unique peptide
processedPeps <- processPeps_TMT(processedPSMs_filtered, modIDs, plex_size)
#processedPeps <- processPeps_isotope(processedPSMs_filtered)

#Write processedPeptides to a file
write.csv(processedPeps, file = paste(udata_path, "insert_file_name_here", Sys.Date(), ".csv", sep = ""), quote = FALSE, row.names = FALSE)


############################## Match structural features to processedPeptides ##############################

#Read in structural features from Uniprot
featureTab <- read.table(paste(sol_metadata_path, "Uniprot_FeatureTable_04.26.22.txt", sep =""), header=TRUE, sep = "\t",  quote="", fill=TRUE, stringsAsFactors = FALSE)
alphaFold <- read.table(paste(sol_metadata_path, "alphaFold_fpocketTab_F1only.txt", sep =""), header=TRUE, sep = "\t",  quote="", fill=TRUE, stringsAsFactors = FALSE)

featureTab <- merge(featureTab,alphaFold, by.x="Entry", by.y = "Uniprot", all.x=TRUE)

#Match structural features to processedPeptides
processedPeps_wFeatures <- matchUniFeatures(processedPeps, featureTab)

#Write processedPeps_wFeatures to file (for use with PDB Processing Script)
write.csv(processedPeps_wFeatures, file = paste(udata_path, "insert_file_name_here", Sys.Date(), ".csv", sep = ""),  quote = FALSE, row.names = FALSE)

##############################################################################




