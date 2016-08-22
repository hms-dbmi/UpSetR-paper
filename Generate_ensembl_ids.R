#This script is for outputting a file with all of the ensembl gene IDs that belong to a specific intersection
#As well as an option to filter by a specific mutation (e.g. ">-" for deletions, or ">G" for a substitution to G)

#Provide a vector of project names that make up a specific intersection, as well as the mutation result 
#(put "none" to get every ensembl gene ID from intersection)
#Will write a file containing a list of ensembl gene IDs related to the mutations in that intersection with a specific mutation
#Results can be uploaded to DAVID
ensembl_id_for_intersection <- function(projects, mutation) {
  intersections <- paste(projects, collapse = " + ")
  ens_id <- myData[which(myData$project == intersections), ]
  if(mutation != "none"){
  ens_id <- ens_id[which(grepl(mutation, ens_id$mutation)), ]
  }
  ens_id <- as.data.frame(unique(as.character(ens_id$ensembl_id)))
  write.table(
    as.data.frame(ens_id),
    "ensemble_id.txt",
    sep = "",
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
  )
}

#Calling the function with these inputs will grab all of the mutations only present in the THCA-SA project
#And then filter down further to grab all of the genes where the mutation was a deletion
ensembl_id_for_intersection("THCA-SA", ">-")