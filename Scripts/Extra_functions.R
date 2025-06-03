#                                             Extra functions ####

#not needed this function but left it here bc its pretty cool

#colnames(hgnc_symbols)<- c('ensembl_id', 'hgnc_symbol')
replace_ensembl_with_hgnc <- function(data, mapping) {
  # Ensure that 'mapping' has the expected columns
  if(!("ensembl_id" %in% names(mapping)) | !("hgnc_symbol" %in% names(mapping))) {
    stop("The mapping data frame must contain 'ensembl_id' and 'hgnc_symbol' columns.")
  }
  
  # Create a named vector for lookup
  lookup_vector <- setNames(mapping$hgnc_symbol, mapping$ensembl_id)
  
  # Replace rownames in 'data' using the lookup vector
  new_rownames <- lookup_vector[rownames(data)]
  
  # Identify unmatched Ensembl IDs
  unmatched_ids <- rownames(data)[is.na(new_rownames)]
  
  if(length(unmatched_ids) > 0) {
    cat("Number of unmatched Ensembl IDs:", length(unmatched_ids), "\n")
    cat("Unmatched Ensembl IDs are:\n")
    print(unmatched_ids)
  }
  
  # Replace NAs with original Ensembl IDs (keeps original IDs if no match is found)
  new_rownames[is.na(new_rownames)] <- unmatched_ids
  
  # Set the new rownames
  rownames(data) <- new_rownames
  
  return(data)
}


#sig_cryo_fresh_symbol<-replace_ensembl_with_hgnc(sig_results_cryo, hgnc_symbols)

#Number of unmatched Ensembl IDs: 74 
#Unmatched Ensembl IDs are:
#  [1] "ENSG00000226380" "ENSG00000278806" "ENSG00000255823" "ENSG00000130723"
#[5] "ENSG00000281657" "ENSG00000275149" "ENSG00000269028" "ENSG00000279645"
#[9] "ENSG00000274152" "ENSG00000259747" "ENSG00000150526" "ENSG00000279906"
#[13] "ENSG00000226598" "ENSG00000278937" "ENSG00000188206" "ENSG00000243587"
#[17] "ENSG00000253379" "ENSG00000271430" "ENSG00000204365" "ENSG00000236673"
#[21] "ENSG00000231429" "ENSG00000261550" "ENSG00000231459" "ENSG00000243444"
#[25] "ENSG00000280156" "ENSG00000254561" "ENSG00000253426" "ENSG00000242349"
#[29] "ENSG00000228614" "ENSG00000130489" "ENSG00000166748" "ENSG00000232559"
#[33] "ENSG00000273038" "ENSG00000233878" "ENSG00000239467" "ENSG00000236269"
#[37] "ENSG00000281881" "ENSG00000185087" "ENSG00000224638" "ENSG00000228651"
#[41] "ENSG00000237647" "ENSG00000263570" "ENSG00000231575" "ENSG00000225986"
#[45] "ENSG00000186354" "ENSG00000273301" "ENSG00000281041" "ENSG00000237404"
#[49] "ENSG00000272367" "ENSG00000226751" "ENSG00000237764" "ENSG00000255090"
#[53] "ENSG00000261662" "ENSG00000225520" "ENSG00000236166" "ENSG00000274956"
#[57] "ENSG00000250599" "ENSG00000250061" "ENSG00000261438" "ENSG00000276182"
#[61] "ENSG00000236740" "ENSG00000266411" "ENSG00000251331" "ENSG00000279184"
#[65] "ENSG00000273237" "ENSG00000261295" "ENSG00000223414" "ENSG00000280551"
#[69] "ENSG00000257794" "ENSG00000269814" "ENSG00000227869" "ENSG00000256863"
#[73] "ENSG00000278309" "ENSG00000275879"

