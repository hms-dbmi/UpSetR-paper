require(jsonlite)
require(curl)
require(plyr)

#An example of the function input. 
#Specify the project name, fields, and number of entries(size) to pull in each list
#This version of pulling from the ICGC REST API specifcally checks whether the mutations affected a coding region in the transcripts field\
#If the mutation affects a non-coding region, it is filtered out of the data frame
data <- list(list(project = "THCA-US", fields = c("id", "mutation", "chromosome", "start", "end", "transcripts"), size = 6659),
             list(project = "THCA-SA", fields = c("id", "mutation", "chromosome", "start", "end", "transcripts"), size = 45126),
             list(project = "LUSC-US", fields = c("id", "mutation", "chromosome", "start", "end", "transcripts"), size = 65063),
             list(project = "LUSC-KR", fields = c("id", "mutation", "chromosome", "start", "end", "transcripts"), size = 64671),
             list(project = "LUSC-CN", fields = c("id", "mutation", "chromosome", "start", "end", "transcripts"), size = 419),
             list(project = "LAML-KR", fields = c("id", "mutation", "chromosome", "start", "end", "transcripts"), size = 42977))

#Function to pull data from ICGC Rest API and convert to UpSetR data frame

icgcData <- function(data){
  aggregateData <- data.frame()
  geneIdCol <- data.frame(ensembl_id = character())
  projectNameCol <- data.frame()
  IDs <- list()
  
  for(i in 1:length(data)){
    IDs[[i]] <- c(0)
    size <- data[[i]]$size
    remainder <- size%%100
    size <- floor(size/100)
    
    if(remainder != 0){
      size <- size + 1
    }
    
    for(j in seq(size)){
      if(j == size && remainder != 0){
        count <- remainder
      }
      else{
        count <- 100
      }
      from <- ((j-1)*100)+1
      partURL <- ""
      
      for(k in 1:length(data[[i]]$fields)){
        field <- data[[i]]$fields[k]
        if(k != 1){
          field <- paste0("%2C", field)
        }
        partURL <- paste0(partURL, field)
      }
      
      url <- paste0("https://dcc.icgc.org:443/api/v1/projects/",
                    data[[i]]$project, "/mutations?field=",
                    partURL, "&&&from=", as.character(from),
                    "&size=",as.character(count),"&&order=desc")
      
      temp <- fromJSON(url)$hits
      geneID_vector <- c()
      length(geneID_vector) <- length(temp$transcripts)
      for(k in 1:length(temp$transcripts)){
        if("protein_coding" %in% temp$transcripts[[k]]$type){
          if (!is.null(temp$transcripts[[k]]$consequence$geneAffectedId)) {
          geneID <- count(temp$transcripts[[k]]$consequence$geneAffectedId)
          names(geneID) <- c("ensembl_id", "freq")
          geneID <- geneID[order(geneID$freq, decreasing = T), ]
          geneID$ensembl_id <- as.character(geneID$ensembl_id)
          geneID <- geneID[1, 1]
          geneID_vector[k] <- geneID
          }
          else{
            geneID_vector[k] <- NA
          }
          temp$transcripts[[k]] <- "protein_coding"
        }
        else {
          if (!is.null(temp$transcripts[[k]]$consequence$geneAffectedId)) {
            geneID <- count(temp$transcripts[[k]]$consequence$geneAffectedId)
            names(geneID) <- c("ensembl_id", "freq")
            geneID <- geneID[order(geneID$freq, decreasing = T),]
            geneID$ensembl_id <- as.character(geneID$ensembl_id)
            geneID <- geneID[1, 1]
            geneID_vector[k] <- geneID
          }
          else{
            geneID_vector[k] <- NA
          }
          temp$transcripts[[k]] <- "non_coding"
        }
      }
      
      geneID_vector <- as.data.frame(matrix(geneID_vector, ncol = 1))
      colnames(geneID_vector) <- "ensembl_id"
      geneIdCol <- rbind(geneIdCol, geneID_vector)
      aggregateData <- rbind(aggregateData, temp)
      IDs[[i]] <- c(IDs[[i]], as.vector(temp$id))
      
      print(paste("On page", as.character(j), "of project", data[[i]]$project))
      
    }
    
    addcolumn <- as.data.frame(matrix(data = rep(data[[i]]$project, data[[i]]$size), nrow = data[[i]]$size, ncol = 1))
    names(addcolumn) <- "project"
    projectNameCol <- rbind(projectNameCol, addcolumn)
  }
  
  aggregateData <- cbind(aggregateData, geneIdCol)
  aggregateData <- cbind(aggregateData, projectNameCol)
  
  num <- nrow(aggregateData)
  projects <- c()
  
  for(i in 1:length(data)){
    name <- data[[i]]$project
    
    setCol <- rep(0, num)
    setCol <- as.data.frame(setCol)
    names(setCol) <- name
    aggregateData <- cbind(aggregateData, setCol)
    projects[i] <- name
    
  }

  aggregateData <- aggregateData[-which(duplicated(aggregateData$id)), ]
  
  for(i in 1:length(data)){
    IDs[[i]] <- IDs[[i]][-1]
    aggregateData[which(aggregateData$id %in% IDs[[i]]), projects[i]] <- 1 
  }
  projects <- apply(aggregateData, 1, function(x){x <- x[-c(1:6)]; x <- names(x[which(x ==  "1")]); x <- paste(unlist(x), collapse = " + ")})
  projects <- as.character(projects)
  aggregateData$project <- projects
  names(aggregateData$project) <- "project"
  return(aggregateData)
}


#Run the function to generate a data frame compatible with UpSetR
myData <- icgcData(data)

#Change transcripts column from list to vector
myData$transcripts <- unlist(myData$transcripts)

#Only grab mutations that affect protein coding regions
myData <- myData[which(myData$transcripts == "protein_coding"), ]

setdata <- data.frame(
  projects= c("THCA-US", "THCA-SA", "ORCA-IN", "BLCA-US",
              "BLCA-CN", "LUSC-CN", "LUSC-US", "LUSC-KR",
              "KIRC-US", "KIRP-US", "LAML-KR"),
  Donors = c(507, 15, 119, 412, 103, 10, 504, 36, 537, 291, 171),
  Site = c("Head & Neck", "Head & Neck", "Head & Neck",
           "Bladder", "Bladder", "Lung", "Lung", "Lung", "Kidney",
           "Kidney", "Blood"),
  Country = c("United States", "Saudi Arabia", "India", "United States",
              "China", "China", "United States", "South Korea",
              "United States", "United States", "South Korea")
)

setdata$Country <- as.character(setdata$Country)
setdata$Site <- as.character(setdata$Site)