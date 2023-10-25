require(bio3d)

#Read in data from SoL Processing Script
udata_path <- file.path("insert_file_path_here")
sol <- read.csv(paste(udata_path,"insert_file_name_here", sep =""), header=TRUE, stringsAsFactors = FALSE)

#Set path to location of all PDB files
pdb_path <- file.path("insert_file_path_here")

sol$PDB_original <- sol$PDB
sol$PDB_clean <- ""

i <- 1
while (i <= dim(sol)[1]){
  
  #those have multiple PDB
  if (nchar(sol$PDB[i]) > 4) {
    #get isoTOP peptides
    pep_seq <- sol[sol$Uniprot == sol$Uniprot[i], ]$PEPTIDE.SEQUENCE
    
    #get PDBs associated with this Uniprot ID
    pdbs_to_check <- strsplit(sol$PDB[i], ";")[[1]]
    number_of_peptide_in_PDB <- numeric()
    
    for (j in 1:length(pdbs_to_check)) {
      #find the PDB that contains most of isoTOP peptide
      number_of_peptide_in_PDB[j] <- tryCatch({
        
        #read PDB file, get its chains and sequences
        pdb <- read.pdb(file.path(pdb_path, paste(
          tolower(pdbs_to_check[j]), ".pdb", sep = ""
        )))
        
        pdb_info <- list(chain = character(),
                         seq = character())
        pdb_info$chain <-
          unique(pdb$atom$chain[pdb$atom$type == "ATOM"])
        for (k in 1:length(pdb_info$chain)) {
          pdb_info$seq[k] <-
            paste(pdbseq(pdb, inds = atom.select(pdb, "calpha", chain = pdb_info$chain[k])), collapse = "")
        }
        
        #determine whether or not each isoTOP peptide is contained in this PDB
        in_pdb <- logical()
        for (n in 1:length(pep_seq)) {
          in_pdb[n] <-
            ifelse(length(grep(pep_seq[n], pdb_info$seq, fixed = TRUE)) > 0, TRUE, FALSE)
        }
        sum(in_pdb)
      }, error = function(e) {
        0
      })
    }
    
    if(sum(number_of_peptide_in_PDB)==0){
      sol[sol$Uniprot == sol$Uniprot[i], ]$PDB_clean <- ""
    }else{
      sol[sol$Uniprot == sol$Uniprot[i], ]$PDB_clean <-
        paste(pdbs_to_check[number_of_peptide_in_PDB == max(number_of_peptide_in_PDB)], collapse = " ")
    }
    
    i <- i + length(pep_seq)
  }else{
    sol$PDB_clean[i] <- sol$PDB[i]
    i <- i + 1
  }
  
  print(i)
}

sol$PDB_quick <- gsub(" .*", "", sol$PDB_clean)

write.csv(sol, file = paste(udata_path, "insert_file_name_here", Sys.Date(), ".csv", sep = ""), row.names = FALSE)