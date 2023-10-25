require(bio3d)

#Read in 16plex Data
udata_path <- file.path("insert_file_path_here")
sol <- read.csv(paste(udata_path,"insert_file_name_here", sep =""), header=TRUE, stringsAsFactors = FALSE)

#Set paths to location of all PDB/AF files
pdb_path <- file.path("/Users/Woz_1/PDB_Downloads/")
af_path <- file.path("/Users/Woz_1/UP000005640_9606_HUMAN_v2/")

#Adjust AF naming
sol$AF <- gsub("\\..*","",sol$AlphaFold)
sol$AF[which(grepl("F1-", sol$AF)==FALSE)]<-""

sol$PDB <- sol$PDB_quick
sol$concatenated_features <- paste(sol$Active.site, sol$Binding.site, sol$Calcium.binding, sol$DNA.binding, sol$Metal.binding, sep = "; ")

sol$active_site <- ""

# Process Uniprot active site info to be compatible with previous developed script
for(i in 1:nrow(sol)){
  
  temp <- NULL
  
  act_indices <- gregexpr("ACT_SITE(.*?);", sol$concatenated_features[i])
  if(act_indices[[1]][1] != -1){
  for(k in 1:length(act_indices[[1]])){
    start <- act_indices[[1]][k]
    mid <- act_indices[[1]][k] + 9
    end <- act_indices[[1]][k] + attr(act_indices[[1]],"match.length")[k] - 2
    
    temp2 <- paste(substr(sol$concatenated_features[i], start, end), substr(sol$concatenated_features[i], mid, end), sep = " ")
    temp <- paste(temp, temp2, sep="; ")
  }
  }
  
  bind_indices <- gregexpr("BINDING(.*?);", sol$concatenated_features[i])
  if(bind_indices[[1]][1] != -1){
  for(k in 1:length(bind_indices[[1]])){
    start <- bind_indices[[1]][k]
    mid <- bind_indices[[1]][k] + 8
    end <- bind_indices[[1]][k] + attr(bind_indices[[1]],"match.length")[k] - 2
    
    temp2 <- paste(substr(sol$concatenated_features[i], start, end), substr(sol$concatenated_features[i], mid, end), sep = " ")
    temp <- paste(temp, temp2, sep="; ")
  }
  }
  
  metal_indices <- gregexpr("METAL(.*?);", sol$concatenated_features[i])
  if(metal_indices[[1]][1] != -1){
  for(k in 1:length(metal_indices[[1]])){
    start <- metal_indices[[1]][k]
    mid <- metal_indices[[1]][k] + 6
    end <- metal_indices[[1]][k] + attr(metal_indices[[1]],"match.length")[k] - 2
    
    temp2 <- paste(substr(sol$concatenated_features[i], start, end), substr(sol$concatenated_features[i], mid, end), sep = " ")
    temp <- paste(temp, temp2, sep="; ")
  }
  }
  
  dna_indices <- gregexpr("DNA_BIND(.*?);", sol$concatenated_features[i])
  if(dna_indices[[1]][1] != -1){
  for(k in 1:length(dna_indices[[1]])){
    start <- dna_indices[[1]][k]
    end <- dna_indices[[1]][k] + attr(dna_indices[[1]],"match.length")[k] - 2
    
    temp2 <- gsub("..", " ", substr(sol$concatenated_features[i], start, end), fixed = TRUE)
    temp <- paste(temp, temp2, sep="; ")
  }
  }
  
  np_indices <- gregexpr("NP_BIND(.*?);", sol$concatenated_features[i])
  if(np_indices[[1]][1] != -1){
  for(k in 1:length(np_indices[[1]])){
    start <- np_indices[[1]][k]
    end <- np_indices[[1]][k] + attr(np_indices[[1]],"match.length")[k] - 2
    
    temp2 <- gsub("..", " ", substr(sol$concatenated_features[i], start, end), fixed = TRUE)
    temp <- paste(temp, temp2, sep="; ")
  }
  }

  if(!is.null(temp)){
    sol$active_site[i] <- substr(temp, 3, nchar(temp))
  }
}


sol$spatial_distance <- ""
sol$closest_site <- ""
sol$shift <- ""
sol$site_total <- ""
sol$site_in_pdb <- ""
sol$site_near_isotop <- ""
sol$site_struct_used <- ""

# Map distances from SoL peptide to active sites
for (i in 1:dim(sol)[1]) {
  max_shift <- 0
  
  try(if ((sol$active_site[i] != "") && (sol$PDB[i] != "")) {
    pdb <- read.pdb(file.path(pdb_path, paste(
      tolower(sol$PDB[i]), ".pdb", sep = ""
    )))
    
    isotop <- sol$PEPTIDE.SEQUENCE[i]
    
    #Position in reference seq
    isotop_start <-
      as.numeric(strsplit(sol$TARG_PEPRANGE[i], "_")[[1]][2])
    
    ref_seq <- strsplit(sol$Sequence[i], "")[[1]]
    
    site_list <-
      strsplit(unlist(strsplit(sol$active_site[i], "; ")), " ")
    site <- character()
    for (j in 1:length(site_list)) {
      site <-
        c(site, as.numeric(site_list[[j]][2]):as.numeric(site_list[[j]][3]))
    }
    site <- as.numeric(site)
    
    pdb_chain <- unique(pdb$atom$chain[pdb$atom$type == "ATOM"])
    min_distance <- character()
    site_distance <- vector(mode = "list", length = length(site))
    names(site_distance) <- site
    
    for (j in 1:length(pdb_chain)) {
      #Read in sequence from PDB file
      chain_seq <-
        pdbseq(pdb, inds = atom.select(pdb, "calpha", chain = pdb_chain[j]))
      
      #Remove negative values from chain sequence if additional AAs are added to protein structure
      try(if(names(chain_seq)[1] < 1){
        chain_start <- which(names(chain_seq)==1)
        chain_end <- as.numeric(names(chain_seq)[length(chain_seq)])
        
        chain_seq <- chain_seq[chain_start:chain_end]
      })
     
      #Position of peptide in PDB chain
      pos <- regexpr(isotop, paste(chain_seq, collapse = ""))
      
      #Determine how different the reference seq vs PDB seq are (need to re-code for proper analysis)
      if (pos != -1)
      {
        #Sequence/position in PDB chain
        seq_pos <- chain_seq[pos:(pos + nchar(isotop) - 1)]
        
        isotop_atoms <-
          atom.select(pdb,
                      "protein",
                      chain = pdb_chain[j],
                      resno = as.numeric(names(seq_pos)))
        
        #Calculate shift in residues from reference seq vs. PDB seq
        shift <- as.numeric(names(seq_pos)[1]) - isotop_start
        
        max_shift <-
          ifelse(abs(shift) >= abs(max_shift),
                 abs(shift),
                 abs(max_shift))
        
        mis_alignment <-
          tryCatch(
            sum(!chain_seq == na.omit(ref_seq[as.numeric(names(chain_seq)) - shift])),
            error = function(e) {
              9999
            }
          )
        
        #show the mismatch residue
        #which(!chain_seq==ref_seq[as.numeric(names(chain_seq))-shift])
        
        if (ifelse(mis_alignment >= 4, TRUE, FALSE)) {
          min_distance[j] <-
            paste("Mismatch", sum(!chain_seq == ref_seq[as.numeric(names(chain_seq)) -
                                                          shift]), sep = ": ")
        }
        else{
          site_atoms <-
            atom.select(pdb,
                        "protein",
                        chain = pdb_chain[j],
                        resno = site + shift)
          if (length(site_atoms$atom) == 0) {
            min_distance[j] <- "Active Site Not Found"
          }
          else{
            min_distance[j] <-
              tryCatch(
                min(dist.xyz(
                  pdb$xyz[isotop_atoms$xyz], pdb$xyz[site_atoms$xyz], all.pairs = TRUE
                )),
                error = function(e) {
                  "ERROR"
                }
              )
            for (k in 1:length(site)) {
              site_sle <-
                atom.select(pdb,
                            "protein",
                            chain = pdb_chain[j],
                            resno = site[k] + shift)
              site_distance[[k]][j] <-
                tryCatch(
                  min(dist.xyz(
                    pdb$xyz[isotop_atoms$xyz], pdb$xyz[site_sle$xyz], all.pairs = TRUE
                  )),
                  error = function(e) {
                    9999
                  }
                )
            }
          }
        }
      } else{
        min_distance[j] <- "isoTOP Peptide Not Found"
      }
    }
    
    site_distance <- lapply(lapply(site_distance, function(x) x[!is.na(x)]), min)
    
    if (any(min_distance == "ERROR")) {
      sol$spatial_distance[i] <- "Something's Wrong!"
    } else if (any(!is.na(as.numeric(min_distance)))) {
      sol$spatial_distance[i] <-
        round(min(as.numeric(min_distance), na.rm = TRUE), digits = 3)
      sol$closest_site[i] <- names(which.min(site_distance)[1])
      sol$shift[i] <- max_shift
      sol$site_total[i] <- length(site)
      sol$site_in_pdb[i] <- sum(site_distance != 9999)
      sol$site_near_isotop[i] <- sum(site_distance <= 6)
      sol$site_struct_used[i] <- "PDB"
    } else{
      sol$spatial_distance[i] <- paste(min_distance, collapse = ";")
      sol$shift[i] <- max_shift
      sol$site_struct_used[i] <- "PDB"
    }
    
  } else if ((sol$active_site[i] != "") && (sol$AF[i] != "")) {
    
    ##################################
    #####       AF Section        ####
    ##################################
    
    pdb <- read.pdb(file.path(af_path, paste(
      sol$AF[i], ".pdb", sep = ""
    )))
    
    isotop <- sol$PEPTIDE.SEQUENCE[i]
    
    #Position in reference seq
    isotop_start <-
      as.numeric(strsplit(sol$TARG_PEPRANGE[i], "_")[[1]][2])
    
    ref_seq <- strsplit(sol$Sequence[i], "")[[1]]
    
    site_list <-
      strsplit(unlist(strsplit(sol$active_site[i], "; ")), " ")
    
    site <- character()
    
    for (j in 1:length(site_list)) {
      site <-
        c(site, as.numeric(site_list[[j]][2]):as.numeric(site_list[[j]][3]))
    }
    
    site <- as.numeric(site)
    
    pdb_chain <- unique(pdb$atom$chain[pdb$atom$type == "ATOM"])
    min_distance <- character()
    site_distance <- vector(mode = "list", length = length(site))
    names(site_distance) <- site
    
    for (j in 1:length(pdb_chain)) {
      #Read in sequence from PDB file
      chain_seq <-
        pdbseq(pdb, inds = atom.select(pdb, "calpha", chain = pdb_chain[j]))
      
      #Remove negative values from chain sequence if additional AAs are added to protein structure
      if(names(chain_seq)[1] < 1){
        chain_start <- which(names(chain_seq)==1)
        chain_end <- as.numeric(names(chain_seq)[length(chain_seq)])
        
        chain_seq <- chain_seq[chain_start:chain_end]
      }
      
      #Position of peptide in PDB chain
      pos <- regexpr(isotop, paste(chain_seq, collapse = ""))
      
      #Determine how different the reference seq vs PDB seq are (need to re-code for proper analysis)
      if (pos != -1)
      {
        #Sequence/position in PDB chain
        seq_pos <- chain_seq[pos:(pos + nchar(isotop) - 1)]
        
        isotop_atoms <-
          atom.select(pdb,
                      "protein",
                      chain = pdb_chain[j],
                      resno = as.numeric(names(seq_pos)))
        
        #Calculate shift in residues from reference seq vs. PDB seq
        shift <- as.numeric(names(seq_pos)[1]) - isotop_start
        
        max_shift <-
          ifelse(abs(shift) >= abs(max_shift),
                 abs(shift),
                 abs(max_shift))
        
        mis_alignment <-
          tryCatch(
            sum(!chain_seq == na.omit(ref_seq[as.numeric(names(chain_seq)) - shift])),
            error = function(e) {
              9999
            }
          )
        
        #show the mismatch residue
        #which(!chain_seq==ref_seq[as.numeric(names(chain_seq))-shift])
        
        if (ifelse(mis_alignment >= 4, TRUE, FALSE)) {
          min_distance[j] <-
            paste("Mismatch", sum(!chain_seq == ref_seq[as.numeric(names(chain_seq)) -
                                                          shift]), sep = ": ")
        }
        else{
          site_atoms <-
            atom.select(pdb,
                        "protein",
                        chain = pdb_chain[j],
                        resno = site + shift)
          if (length(site_atoms$atom) == 0) {
            min_distance[j] <- "Active Site Not Found"
          }
          else{
            min_distance[j] <-
              tryCatch(
                min(dist.xyz(
                  pdb$xyz[isotop_atoms$xyz], pdb$xyz[site_atoms$xyz], all.pairs = TRUE
                )),
                error = function(e) {
                  "ERROR"
                }
              )
            for (k in 1:length(site)) {
              site_sle <-
                atom.select(pdb,
                            "protein",
                            chain = pdb_chain[j],
                            resno = site[k] + shift)
              
              site_distance[[k]][j] <-
                tryCatch(
                  min(dist.xyz(
                    pdb$xyz[isotop_atoms$xyz], pdb$xyz[site_sle$xyz], all.pairs = TRUE
                  )),
                  error = function(e) {
                    9999
                  }
                )
            }
          }
        }
      }
      else{
        min_distance[j] <- "isoTOP Peptide Not Found"
      }
    }
    
    site_distance <- lapply(site_distance, min)
    
    if (any(min_distance == "ERROR")) {
      sol$spatial_distance[i] <- "Something's Wrong!"
    } else if (any(!is.na(as.numeric(min_distance)))) {
      sol$spatial_distance[i] <-
        round(min(as.numeric(min_distance), na.rm = TRUE), digits = 3)
      sol$closest_site[i] <- names(which.min(site_distance)[1])
      sol$shift[i] <- max_shift
      sol$site_total[i] <- length(site)
      sol$site_in_pdb[i] <- sum(site_distance != 9999)
      sol$site_near_isotop[i] <- sum(site_distance <= 6)
      sol$site_struct_used[i] <- "AF"
    } else{
      sol$spatial_distance[i] <- paste(min_distance, collapse = ";")
      sol$shift[i] <- max_shift
      sol$site_struct_used[i] <- "AF"
    }
    
  }
  )
  
  #Rerun AF part if active site or peptide not found in PDB
  try(if ((grepl("Not Found", sol$spatial_distance[i]))){
    pdb <- read.pdb(file.path(af_path, paste(
      sol$AF[i], ".pdb", sep = ""
    )))
    
    isotop <- sol$PEPTIDE.SEQUENCE[i]
    
    #Position in reference seq
    isotop_start <-
      as.numeric(strsplit(sol$TARG_PEPRANGE[i], "_")[[1]][2])
    
    ref_seq <- strsplit(sol$Sequence[i], "")[[1]]
    
    site_list <-
      strsplit(unlist(strsplit(sol$active_site[i], "; ")), " ")
    
    site <- character()
    
    for (j in 1:length(site_list)) {
      site <-
        c(site, as.numeric(site_list[[j]][2]):as.numeric(site_list[[j]][3]))
    }
    
    site <- as.numeric(site)
    
    pdb_chain <- unique(pdb$atom$chain[pdb$atom$type == "ATOM"])
    min_distance <- character()
    site_distance <- vector(mode = "list", length = length(site))
    names(site_distance) <- site
    
    for (j in 1:length(pdb_chain)) {
      #Read in sequence from PDB file
      chain_seq <-
        pdbseq(pdb, inds = atom.select(pdb, "calpha", chain = pdb_chain[j]))
      
      #Remove negative values from chain sequence if additional AAs are added to protein structure
      if(names(chain_seq)[1] < 1){
        chain_start <- which(names(chain_seq)==1)
        chain_end <- as.numeric(names(chain_seq)[length(chain_seq)])
        
        chain_seq <- chain_seq[chain_start:chain_end]
      }
      
      #Position of peptide in PDB chain
      pos <- regexpr(isotop, paste(chain_seq, collapse = ""))
      
      #Determine how different the reference seq vs PDB seq are (need to re-code for proper analysis)
      if (pos != -1)
      {
        #Sequence/position in PDB chain
        seq_pos <- chain_seq[pos:(pos + nchar(isotop) - 1)]
        
        isotop_atoms <-
          atom.select(pdb,
                      "protein",
                      chain = pdb_chain[j],
                      resno = as.numeric(names(seq_pos)))
        
        #Calculate shift in residues from reference seq vs. PDB seq
        shift <- as.numeric(names(seq_pos)[1]) - isotop_start
        
        max_shift <-
          ifelse(abs(shift) >= abs(max_shift),
                 abs(shift),
                 abs(max_shift))
        
        mis_alignment <-
          tryCatch(
            sum(!chain_seq == na.omit(ref_seq[as.numeric(names(chain_seq)) - shift])),
            error = function(e) {
              9999
            }
          )
        
        #show the mismatch residue
        #which(!chain_seq==ref_seq[as.numeric(names(chain_seq))-shift])
        
        if (ifelse(mis_alignment >= 4, TRUE, FALSE)) {
          min_distance[j] <-
            paste("Mismatch", sum(!chain_seq == ref_seq[as.numeric(names(chain_seq)) -
                                                          shift]), sep = ": ")
        }
        else{
          site_atoms <-
            atom.select(pdb,
                        "protein",
                        chain = pdb_chain[j],
                        resno = site + shift)
          if (length(site_atoms$atom) == 0) {
            min_distance[j] <- "Active Site Not Found"
          }
          else{
            min_distance[j] <-
              tryCatch(
                min(dist.xyz(
                  pdb$xyz[isotop_atoms$xyz], pdb$xyz[site_atoms$xyz], all.pairs = TRUE
                )),
                error = function(e) {
                  "ERROR"
                }
              )
            for (k in 1:length(site)) {
              site_sle <-
                atom.select(pdb,
                            "protein",
                            chain = pdb_chain[j],
                            resno = site[k] + shift)
              site_distance[[k]][j] <-
                tryCatch(
                  min(dist.xyz(
                    pdb$xyz[isotop_atoms$xyz], pdb$xyz[site_sle$xyz], all.pairs = TRUE
                  )),
                  error = function(e) {
                    9999
                  }
                )
            }
          }
        }
      }
      else{
        min_distance[j] <- "isoTOP Peptide Not Found"
      }
    }
    
    site_distance <- lapply(site_distance, min)
    
    if (any(min_distance == "ERROR")) {
      sol$spatial_distance[i] <- "Something's Wrong!"
    } else if (any(!is.na(as.numeric(min_distance)))) {
      sol$spatial_distance[i] <-
        round(min(as.numeric(min_distance), na.rm = TRUE), digits = 3)
      sol$closest_site[i] <- names(which.min(site_distance)[1])
      sol$shift[i] <- max_shift
      sol$site_total[i] <- length(site)
      sol$site_in_pdb[i] <- sum(site_distance != 9999)
      sol$site_near_isotop[i] <- sum(site_distance <= 6)
      sol$site_struct_used[i] <- "AF"
    } else{
      sol$spatial_distance[i] <- paste(min_distance, collapse = ";")
      sol$shift[i] <- max_shift
      sol$site_struct_used[i] <- "AF"
    }
  })
  
  print(i)
}

write.csv(sol, file = paste(udata_path, "insert_file_name_here", Sys.Date(), ".csv", sep = ""), row.names = FALSE)


