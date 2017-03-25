require(jsonlite)
require(curl)
require(plyr)

#An example of the function input. 
#Specify the project name, fields, and number of entries(size) to pull in each list

data <- list(list(project = "THCA-US", fields = c("id", "mutation", "chromosome", "start", "end"), size = 6659),
             list(project = "ORCA-IN", fields = c("id", "mutation", "chromosome", "start", "end"), size = 13626),
             list(project = "THCA-SA", fields = c("id", "mutation", "chromosome", "start", "end"), size = 45126),
             list(project = "LUSC-US", fields = c("id", "mutation", "chromosome", "start", "end"), size = 65063),
             list(project = "LUSC-KR", fields = c("id", "mutation", "chromosome", "start", "end"), size = 64671),
             list(project = "LUSC-CN", fields = c("id", "mutation", "chromosome", "start", "end"), size = 419),
             list(project = "KIRC-US", fields = c("id", "mutation", "chromosome", "start", "end"), size = 26371),
             list(project = "LAML-KR", fields = c("id", "mutation", "chromosome", "start", "end"), size = 42977))

#Function to pull data from ICGC Rest API and convert to UpSetR data frame

icgcData <- function(data){
  aggregateData <- data.frame()
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
      aggregateData <- rbind(aggregateData, temp)
      IDs[[i]] <- c(IDs[[i]], as.vector(temp$id))
      
      print(paste("On page", as.character(j), "of project", data[[i]]$project))
      
    }
    
    addcolumn <- as.data.frame(matrix(data = rep(data[[i]]$project, data[[i]]$size), nrow = data[[i]]$size, ncol = 1))
    names(addcolumn) <- "project"
    projectNameCol <- rbind(projectNameCol, addcolumn)
  }
  
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
