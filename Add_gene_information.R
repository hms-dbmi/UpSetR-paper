#This script will use biomaRt to map our ensembl gene ids affected by the mutations to a gene name
#As well as a description of the genes function

library("biomaRt")
#set up a biomart for human genes
human <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "www.ensembl.org")
#grab all of the ensemb IDs contained in the data pulled from ICGC
ensembl_genes <- as.character(myData$ensembl_id)

#Using the ensembl IDs from the mutations pulled via the ICGC API, find their corresponding gene name
#As well as a description of the gene function
gene_data <- getBM(
  filters= "ensembl_gene_id", 
  attributes= c("ensembl_gene_id", "external_gene_name", "description"),
  values= ensembl_genes,
  mart= human)

#All of the mutations that affected a gene with an ensembl gene ID not recognized in the human biomart
no_info <- myData[which(!myData$ensembl_id %in% gene_data$ensembl_gene_id), ]

#All recognized ensembl IDs
myData <- myData[which(myData$ensembl_id %in% gene_data$ensembl_gene_id), ]

#Prep columns on gene name and description information to be added on to myData
blanks <- as.data.frame(matrix(rep("", nrow(myData)*2), ncol = 2))
names(blanks) <- c("gene_name", "description")
myData$gene_name <- as.character(myData$gene_name)
myData$description <- as.character(myData$description)
myData <- cbind(myData, blanks)

num_of_genes_recognized <- nrow(gene_data)

#This loop will take a few minutes
for(i in seq(num_of_genes_recognized)){
  print(paste("Getting info on gene", i, "of", num_of_genes_recognized))
  myData[which(myData$ensembl_id == gene_data$ensembl_gene_id[i]), ]$gene_name <- gene_data$external_gene_name[i]
  myData[which(myData$ensembl_id == gene_data$ensembl_gene_id[i]), ]$description <- gene_data$description[i]
}