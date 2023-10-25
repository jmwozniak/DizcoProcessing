require(bio3d)

#get the residues around the binding pockets
get_resid <- function(path,
                      pdb_id,
                      pep_seq,
                      isotop_start,
                      isotop_end) {
  #read the pdb file and obtain information of every chain and sequence
  pdb <-
    read.pdb(file.path(dirname(path), tolower(paste(
      pdb_id, ".pdb", sep = ""
    ))), maxlines = -1)
  
  pdb_info <- list(
    chain = character(),
    seq = character(),
    shift = numeric(),
    target = c()
  )
  
  #Get chains from pdb file
  pdb_info$chain <- unique(pdb$atom$chain[pdb$atom$type == "ATOM"])
  
  #Get chains from pdb file
  for (i in 1:length(pdb_info$chain)) {
    #Get sequences from each chain
    pdb_info$seq[i] <- 
      paste(pdbseq(pdb, inds = atom.select(pdb, "calpha", chain = pdb_info$chain[i])), collapse = "")
    
    #Determine which, if any, chain contains pep_seq
    pdb_info$target[i] <-
      ifelse(length(grep(pep_seq, pdb_info$seq[i], fixed = TRUE)) != 0, TRUE, FALSE)
    
    #Determine shift (in # of AAs) of pdb file from actual protein location
    if (pdb_info$target[i] == TRUE) {
      pdb_info$shift[i] <-
        as.numeric(names(pdbseq(
          pdb, inds = atom.select(pdb, "calpha", chain = pdb_info$chain[i])
        )[regexpr(pep_seq, pdb_info$seq[i])])) - isotop_start
    } else{
      pdb_info$shift[i] <- NA
    }
    
  }
  
  #If not found in any chain, return
  if (any(pdb_info$target) == FALSE) {
    return(paste("Isotop Peptide Not Found in ", pdb_id, sep = ""))
  }
  
  #Only consider those chains that contain isotop peptide sequence
  chain_of_target <- pdb_info$chain[pdb_info$target]
  number_of_pockets <- length(list.files(file.path(path, "pockets"))) /
    2
  
  if (number_of_pockets > 0) {
    #Create list to store fpocket residues and name of chain
    resid_by_fpocket <- replicate(length(chain_of_target), character())
    names(resid_by_fpocket) <- chain_of_target
    
    #Get info on all pockets from fpocket analysis
    pocket.info <-
      read.delim(
        file.path(path, tolower(paste(pdb_id, "_info.txt", sep = ""
        ))),
        header = FALSE,
        stringsAsFactors = FALSE,
        strip.white = TRUE
      )
    
    #Extract all pocket volumes and filter by volume cutoff
    pocket_vol <- numeric()
    for (i in 1:number_of_pockets) {
      pocket_vol <- c(pocket_vol, pocket.info$V3[20 * i - 12])
    }
    pockets_ind <- (1:number_of_pockets)[pocket_vol >= vol_cutoff]
    
    
    if (length(pockets_ind) > 0) {
      for (i in pockets_ind) {
        
        # Read in pockets from fpocket output
        pdb <-
          read.pdb(file.path(
            path,
            "pockets",
            paste("pocket", i, "_atm.pdb", sep = "")))
        
        # Find the chains of each pocket
        chain_of_pocket <-
          unique(pdb$atom$chain[pdb$atom$type == "ATOM"])
        
        
        if ("FALSE" %in% as.character(chain_of_pocket)) {
          # If peptide isn't in chain, return F
          chain_of_pocket[as.character(chain_of_pocket) == "FALSE"] <- "F"
          chain_to_go <- intersect(chain_of_target, chain_of_pocket)
          
          # If peptide isn't in chain, return F
          if (length(chain_to_go) != 0) {
            for (j in 1:length(chain_to_go)) {
              if (chain_to_go[j] == "F") {
                res <- pdbseq(pdb, inds = atom.select(pdb, chain = "FALSE", type = "ATOM"))
                res <- res[!duplicated(res)]
                resid_by_fpocket[[chain_to_go[j]]] <-
                  union(resid_by_fpocket[[chain_to_go[j]]],
                        paste(names(res), chain_to_go[j], sep = "."))
              } else{
                res <-
                  pdbseq(pdb,
                         inds = atom.select(pdb, chain = chain_to_go[j], type = "ATOM"))
                res <- res[!duplicated(res)]
                resid_by_fpocket[[chain_to_go[j]]] <-
                  union(resid_by_fpocket[[chain_to_go[j]]],
                        paste(names(res), chain_to_go[j], sep = "."))
              }
            }
          }
        } else{
          chain_to_go <-
            as.character(intersect(chain_of_target, chain_of_pocket))
          if (length(chain_to_go) != 0) {
            for (j in 1:length(chain_to_go)) {
              res <-
                pdbseq(pdb,
                       inds = atom.select(pdb, chain = chain_to_go[j], type = "ATOM"))
              res <- res[!duplicated(res)]
              resid_by_fpocket[[chain_to_go[j]]] <-
                union(resid_by_fpocket[[chain_to_go[j]]],
                      paste(names(res), chain_to_go[j], sep = "."))
            }
          }
        }
      }
    } else{
      return("No Pockets Pass Filter")
    }
  } else{
    return("No Predicted Pockets")
  }
  
  
  overlapped <- replicate(length(resid_by_fpocket), character())
  names(overlapped) <- names(resid_by_fpocket)
  shifted_isotop <- replicate(length(resid_by_fpocket), character())
  names(shifted_isotop) <- names(resid_by_fpocket)
  
  for (i in 1:length(resid_by_fpocket)) {
    
    shifted_isotop[[names(resid_by_fpocket[i])]] <-
      paste((isotop_start:(isotop_end)) + pdb_info$shift[pdb_info$chain == names(resid_by_fpocket[i])],
            names(resid_by_fpocket[i]),
            sep = ".")
    
    overlapped[[names(resid_by_fpocket[i])]] <-
      intersect(shifted_isotop[[names(resid_by_fpocket[i])]], resid_by_fpocket[[i]])
  }
  
  resid_by_fpocket <-
    lapply(Filter(length, resid_by_fpocket), paste, collapse = ";")
  overlapped <- lapply(Filter(length, overlapped), paste, collapse = ";")
  shifted_isotop <-
    lapply(Filter(length, shifted_isotop), paste, collapse = ";")
  
  if (length(overlapped) != 0) {
    return(c(
      paste(paste(resid_by_fpocket, collapse = ";"), sep = " "),
      paste(paste(overlapped, collapse = ";"), sep = " "),
      paste(paste(shifted_isotop, collapse = ";"), sep = " "),
      length(pockets_ind)
    ))
  } else{
    return(c(
      paste(paste(resid_by_fpocket, collapse = ";"), sep = " "),
      "No Overlap",
      paste(paste(shifted_isotop, collapse = ";"), sep = " "),
      length(pockets_ind)
    ))
  }
}

get_pocket <- function(path, resid.overlapped) {
  
  #determine the number of pockets of the pdb file
  number_of_pockets <- length(list.files(file.path(path, "pockets"))) /
    2
  
  pocket.info <-
    read.delim(
      file.path(path, tolower(paste(pdb_id, "_info.txt", sep = ""
      ))),
      header = FALSE,
      stringsAsFactors = FALSE,
      strip.white = TRUE
    )
  
  pocket_vol <- numeric()
  for (i in 1:number_of_pockets) {
    pocket_vol <- c(pocket_vol, pocket.info$V3[20 * i - 12])
  }
  
  pockets_ind <- (1:number_of_pockets)[pocket_vol >= vol_cutoff]
  
  resid.overlapped <- strsplit(resid.overlapped, ";")[[1]]
  overlap_w_pocket <- c()
  
  for (i in pockets_ind) {
    
    pdb <-
      read.pdb(file.path(path, "pockets", paste("pocket", i, "_atm.pdb", sep = "")))
    
    chain_of_pocket <- unique(pdb$atom$chain[pdb$atom$type == "ATOM"])
    resid_by_pocket <- c()
    
    for (j in 1:length(chain_of_pocket)) {
      res <-
        pdbseq(pdb, inds = atom.select(
          pdb,
          chain = as.character(chain_of_pocket[j]),
          type = "ATOM"
        ))
      res <- res[!duplicated(res)]
      if (as.character(chain_of_pocket[j]) == "FALSE") {
        resid_by_pocket <- c(resid_by_pocket, paste(names(res), "F", sep = "."))
      } else{
        resid_by_pocket <-
          c(resid_by_pocket,
            paste(names(res), chain_of_pocket[j], sep = "."))
      }
    }
    if (length(intersect(resid.overlapped, resid_by_pocket)) > 0) {
      overlap_w_pocket <- c(overlap_w_pocket, i)
    }
  }
  return(overlap_w_pocket)
}

#Read in data from PDB filtering script
udata_path <- file.path("insert_file_path_here")
fpocket <- read.csv(paste(udata_path,"insert_file_name_here", sep =""), header=TRUE, stringsAsFactors = FALSE)

#Set paths to location of all PDB/AF files
af_path <- file.path("/Users/Woz_1/UP000005640_9606_HUMAN_v2/")
pdb_path <- file.path("/Users/Woz_1/PDB_Downloads/")

fpocket$PDB <- fpocket$PDB_quick

fpocket$resid_by_fpocket <- ""
fpocket$shifted_isotop <- ""
fpocket$overlapped <- ""
fpocket$isotop.covered <- ""
fpocket$no.filtered.pockets <- ""
fpocket$number.of.pockets <- "-"
fpocket$pocket.overlapped <- "-"
fpocket$pocket_struct_used <- "-"

vol_cutoff <- 0

all.pocket.info <- NULL
total_pockets <- 1

for (k in 1:dim(fpocket)[1]) {
  try(
  if(fpocket$PDB[k] != ""){
    
    path <- file.path(pdb_path, tolower(paste(fpocket$PDB[k], "out", sep = "_")))

    isotop_start <-
      as.numeric(strsplit(fpocket$Labeled.Peptide[k], "_")[[1]][1])
    
    isotop_end <-
      as.numeric(strsplit(fpocket$Labeled.Peptide[k], "_")[[1]][2])
    
    pdb_id <- fpocket$PDB[k]
    pep_seq <- fpocket$PEPTIDE.SEQUENCE[k]
    
    if(is.na(isotop_start) || is.na(isotop_start) || is.na(pep_seq)){
      
      fpocket$resid_by_fpocket[k] <- NA
      fpocket$overlapped[k] <- "No Overlap"
      fpocket$shifted_isotop[k] <- NA
      fpocket$no.filtered.pockets[k] <- NA
      
    }else{
      result <- get_resid(path, pdb_id, pep_seq, isotop_start, isotop_end)
      
      fpocket$resid_by_fpocket[k] <- result[1]
      fpocket$overlapped[k] <- result[2]
      fpocket$shifted_isotop[k] <- result[3]
      fpocket$no.filtered.pockets[k] <- result[4]
      
      if(is.na(fpocket$overlapped[k])) fpocket$overlapped[k] <- "No Overlap"
    }
    
    if (fpocket$overlapped[k] != "No Overlap") {
      overlapped_resid <- strsplit(fpocket$overlapped[k], ";")[[1]]
      
      mean_overlapped <-
        mean(table(gsub("[0-9]+.", "", overlapped_resid)))
      
      isotop_peptide_length <- nchar(pep_seq)
      
      fpocket$isotop.covered[k] <-
        round(mean_overlapped / isotop_peptide_length, digits = 3)
      
    } else{
      fpocket$isotop.covered[k] <- 0
    }
    
    fpocket$number.of.pockets[k] <-
      length(list.files(file.path(path, "pockets"))) / 2
    
    if (fpocket$overlapped[k] != "No Overlap") {
      fpocket$pocket.overlapped[k] <-
        toString(get_pocket(path, fpocket$overlapped[k]))
     
      temp.pocket.info <-
        read.delim(
          file.path(path, paste(pdb_id, "_info.txt", sep = ""
          )),
          header = FALSE,
          stringsAsFactors = FALSE,
          strip.white = TRUE
        )
      temp.pocket.info <- cbind(pdb_id, temp.pocket.info)
      all.pocket.info <- rbind(all.pocket.info, temp.pocket.info)
    }
    
    fpocket$pocket_struct_used[k] <- "PDB"
    
    print(k)
    
    #Use AF part if peptide not found in PDB
    if (grepl("Not Found", fpocket$resid_by_fpocket[k])){
      path <- file.path(af_path, paste(fpocket$AF[k], "out", sep = "_"))
      
      isotop_start <-
        as.numeric(strsplit(fpocket$Labeled.Peptide[k], "_")[[1]][1])
      
      isotop_end <-
        as.numeric(strsplit(fpocket$Labeled.Peptide[k], "_")[[1]][2])
      
      pdb_id <- fpocket$AF[k]
      pep_seq <- fpocket$PEPTIDE.SEQUENCE[k]
      
      if(is.na(isotop_start) || is.na(isotop_start) || is.na(pep_seq)){
        
        fpocket$resid_by_fpocket[k] <- NA
        fpocket$overlapped[k] <- "No Overlap"
        fpocket$shifted_isotop[k] <- NA
        fpocket$no.filtered.pockets[k] <- NA
        
      }else{
        result <- get_resid(path, pdb_id, pep_seq, isotop_start, isotop_end)
        
        fpocket$resid_by_fpocket[k] <- result[1]
        fpocket$overlapped[k] <- result[2]
        fpocket$shifted_isotop[k] <- result[3]
        fpocket$no.filtered.pockets[k] <- result[4]
        
        if(is.na(fpocket$overlapped[k])) fpocket$overlapped[k] <- "No Overlap"
      }
      
      if (fpocket$overlapped[k] != "No Overlap") {
        overlapped_resid <- strsplit(fpocket$overlapped[k], ";")[[1]]
        
        mean_overlapped <-
          mean(table(gsub("[0-9]+.", "", overlapped_resid)))
        
        isotop_peptide_length <- nchar(pep_seq)
        
        fpocket$isotop.covered[k] <-
          round(mean_overlapped / isotop_peptide_length, digits = 3)
        
      } else{
        fpocket$isotop.covered[k] <- 0
      }
      
      fpocket$number.of.pockets[k] <-
        length(list.files(file.path(path, "pockets"))) / 2
      
      if (fpocket$overlapped[i] != "No Overlap") {
        fpocket$pocket.overlapped[i] <-
          toString(get_pocket(path, fpocket$overlapped[i]))
        
        temp.pocket.info <-
          read.delim(
            file.path(path, paste(pdb_id, "_info.txt", sep = ""
            )),
            header = FALSE,
            stringsAsFactors = FALSE,
            strip.white = TRUE
          )
        temp.pocket.info <- cbind(pdb_id, temp.pocket.info)
        all.pocket.info <- rbind(all.pocket.info, temp.pocket.info)
      }
      
      print(k)
      
      fpocket$pocket_struct_used[k] <- "AF"
    }
    
  } else if (fpocket$AF[k] != ""){
    
    ##################################
    #####       AF Section        ####
    ##################################
    
    path <- file.path(af_path, paste(fpocket$AF[k], "out", sep = "_"))
    
    isotop_start <-
      as.numeric(strsplit(fpocket$Labeled.Peptide[k], "_")[[1]][1])
    
    isotop_end <-
      as.numeric(strsplit(fpocket$Labeled.Peptide[k], "_")[[1]][2])
    
    pdb_id <- fpocket$AF[k]
    pep_seq <- fpocket$PEPTIDE.SEQUENCE[k]
    
    if(is.na(isotop_start) || is.na(isotop_start) || is.na(pep_seq)){
      
      fpocket$resid_by_fpocket[k] <- NA
      fpocket$overlapped[k] <- "No Overlap"
      fpocket$shifted_isotop[k] <- NA
      fpocket$no.filtered.pockets[k] <- NA
      
    }else{
      result <- get_resid(path, pdb_id, pep_seq, isotop_start, isotop_end)
      
      fpocket$resid_by_fpocket[k] <- result[1]
      fpocket$overlapped[k] <- result[2]
      fpocket$shifted_isotop[k] <- result[3]
      fpocket$no.filtered.pockets[k] <- result[4]
      
      if(is.na(fpocket$overlapped[k])) fpocket$overlapped[k] <- "No Overlap"
    }
    
    if (fpocket$overlapped[k] != "No Overlap") {
      overlapped_resid <- strsplit(fpocket$overlapped[k], ";")[[1]]
      
      mean_overlapped <-
        mean(table(gsub("[0-9]+.", "", overlapped_resid)))
      
      isotop_peptide_length <- nchar(pep_seq)
      
      fpocket$isotop.covered[k] <-
        round(mean_overlapped / isotop_peptide_length, digits = 3)
      
    } else{
      fpocket$isotop.covered[k] <- 0
    }
    
    fpocket$number.of.pockets[k] <-
      length(list.files(file.path(path, "pockets"))) / 2
    
    if (fpocket$overlapped[k] != "No Overlap") {
      fpocket$pocket.overlapped[k] <-
        toString(get_pocket(path, fpocket$overlapped[k]))
      
      temp.pocket.info <-
        read.delim(
          file.path(path, paste(pdb_id, "_info.txt", sep = ""
          )),
          header = FALSE,
          stringsAsFactors = FALSE,
          strip.white = TRUE
        )
      
      temp.pocket.info <- cbind(pdb_id, temp.pocket.info)
      all.pocket.info <- rbind(all.pocket.info, temp.pocket.info)
    }
    
    print(k)
    
    fpocket$pocket_struct_used[k] <- "AF"
  }
  )
  
  print(k)
}

write.csv(fpocket, file = paste(udata_path, "insert_file_name_here", Sys.Date(), ".csv", sep = ""),  quote = FALSE, row.names = FALSE)
