myData <- read.table("icgcData.txt", header = T, sep = "\t", check.names = F)

library(UpSetR)

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

myplot <- function(data, colour){
  data <- data[which(data$color == colour), ]
  plot_title <- as.character(unique(data$project))
  data <- count(data["mutation"])
  data$freq <- as.numeric(data$freq)
  data$mutation <- as.character(data$mutation)
  data <- data[which(nchar(data$mutation) == 3), ]
  data <- data[order(data$mutation), ]
  bases <- strsplit(data$mutation, ">")
  original <- unlist(lapply(bases, function(x){x <- x[1]}))
  mutation <- unlist(lapply(bases, function(x){x <- x[2]}))
  data$original <- original
  data$mutation <- mutation
  
  if(length(data[which(data$original == "-"), ]$original)!=0){
    data[which(data$original == "-"), ]$original <- "__"
  }
  if(length(data[which(data$mutation == "-"), ]$mutation)!= 0){
    data[which(data$mutation == "-"), ]$mutation <- "__"
  }
  check <- c("__", "G", "C", "T", "A")
  missing <- check[!check %in% data$mutation]
  add <- c()
  if(length(missing) != 0){
    for(i in 1:length(missing)){
      add[i] <- missing[i]
    }
    add_original <- rep("A", length(missing))
    add_freq <- rep(NA, length(missing))
    add_data <- data.frame(mutation = add, freq = add_freq, original = add_original)
    data <- rbind(data, add_data)
  }
  add <- c()
  add_freq <- c()
  missing <- check[!check %in% data$original]
  if(length(missing) != 0){
    for(i in 1:length(missing)){
      add[i] <- missing[i]
    }
    add_mutation <- rep("A", length(missing))
    add_freq <- rep(NA, length(missing))
    add_data <- data.frame(mutation = add_mutation, freq = add_freq, original = add)
    data <- rbind(data, add_data)
    data[which(is.na(data$freq)), ]$freq <- 0
  }
  
  maxvalue <- max(data$freq)
  
  plot <- (ggplot(data, aes_string(x="mutation", y = "original")) + theme(panel.background = element_rect("white"),
                                                                          axis.ticks = element_blank(),
                                                                          plot.margin = unit(c(0,1.5,0.3,1.5), "cm"),
                                                                          plot.title = element_text(color=colour, size = 10, hjust = 0.5),
                                                                          axis.title.x = element_text(size = 8.3),
                                                                          axis.title.y = element_text(size = 8.3),
                                                                          axis.text.x = element_text(size = 7),
                                                                          axis.text.y = element_text(size = 7),
                                                                          legend.key.size = unit(0.4, "cm"),
                                                                          legend.text = element_text(size = 7),
                                                                          legend.title = element_text(size = 8.3),
                                                                          legend.key.height = unit(0.4, "cm"))
           + geom_tile(aes(fill = freq)) + scale_fill_gradient(low="grey86", high = "grey16", na.value = "white", limits = c(0, maxvalue))
           + ggtitle(plot_title))
  
}

pdf(file = "supfig1.pdf", width = 8, height = 5)
# sup fig 1
upset(myData, nsets = 8, nintersects = 30, order.by = "freq")
dev.off()

pdf(file = "supfig4.pdf", width = 8, height = 5)
# sup fig 4
upset(myData, nsets = 8,
      sets.x.label = "Mutations", mainbar.y.label = "Shared Mutations",
      order.by = "freq", nintersects = 30,
      queries = list(list(query = intersects, params = list("LUSC-US"), color = "#56b4e9"),
                     list(query = intersects, params = list("THCA-SA"), color = "#cc79a9"),
                     list(query = intersects, params = list("LUSC-KR", "LAML-KR"), color = "#009e73")))
dev.off()

# sup fig 5
typeOfMutation <- as.character(myData$mutation)
typeOfMutation <- strsplit(typeOfMutation, split = ">")
typeOfMutation <- unlist(lapply(typeOfMutation, function(x) {
  if(x[2] == "-") {
    return("deletion")
  }
  if(x[1] == "-") {
    return ("insertion")
  }
  else {
    return("substitution")
  }
}))
myData$mutation_type <- typeOfMutation

pdf(file = "supfig5.pdf", width = 8, height = 5)
upset(myData, nsets = 8, mb.ratio = c(0.65, 0.35),
      sets.x.label = "Mutations", mainbar.y.label = "Shared Mutations",
      order.by = "freq", nintersects = 30,
      queries = list(list(query = intersects, params = list("LUSC-US"), color = "#56b4e9"),
                     list(query = intersects, params = list("THCA-SA"), color = "#cc79a7"),
                     list(query = intersects, params = list("LUSC-KR", "LAML-KR"), color = "#009e73"),
                     list(query = elements, params = list("mutation_type", "deletion"), color = "#e69f00", active = T)))
dev.off()

mutationTypeHistogram <- function(data, project) {
  countData <- count(data$mutation_type)
  histData <- data.frame(x = as.character(countData$x), count = countData$freq,
                         check.names = F, stringsAsFactors = F)
  SaudiData <- data[which(data$project == project), ]
  SaudiCountData <- count(SaudiData$mutation_type)
  SaudiHistData <- data.frame(x = as.character(SaudiCountData$x), count = SaudiCountData$freq,
                              check.names = F, stringsAsFactors = F)
  plot <- (ggplot(histData, aes_string(x="x", y = "count")) + theme(panel.background = element_rect("white"),
                                                                    axis.ticks = element_blank(),
                                                                    plot.margin = unit(c(0,1.5,0.3,1.5), "cm"),
                                                                    plot.title = element_text(size = 10, hjust = 0.5),
                                                                    axis.title.x = element_text(size = 8.3),
                                                                    axis.title.y = element_text(size = 8.3),
                                                                    axis.text.x = element_text(size = 7),
                                                                    axis.text.y = element_text(size = 7),
                                                                    legend.key.size = unit(0.4, "cm"),
                                                                    legend.text = element_text(size = 7),
                                                                    legend.title = element_text(size = 8.3),
                                                                    legend.key.height = unit(0.4, "cm"))
           + geom_bar(stat = "identity", width = 0.6, position = "identity")
           + geom_bar(data = SaudiHistData,
                      aes_string(x="x", y = "count"),
                      fill = "#cc79a7",
                      stat = "identity", position = "identity", width = 0.6)
           + xlab("mutation type")
           + ylab("count")
           + ggtitle("Counts per mutation type"))
}

pdf(file = "fig1.pdf", width = 14, height = 6.2)
#Figure 1
upset(myData, nsets = 8,
      sets.x.label = "Mutations", mainbar.y.label = "Shared Mutations", mb.ratio = c(0.6,0.4),
      order.by = "freq", nintersects = 30,
      queries = list(list(query = intersects, params = list("LUSC-US"), color = "#56b4e9"),
                     list(query = intersects, params = list("THCA-SA"), color = "#cc79a7"),
                     list(query = intersects, params = list("LUSC-KR", "LAML-KR"), color = "#009e73"),
                     list(query = elements, params = list("mutation_type", "deletion"), color = "#e69f00", active = T)),
      attribute.plots = list(gridrows = 40,
                             plots = list(list(plot = myplot, x = "#56b4e9", queries = T),
                                          list(plot = myplot, x = "#cc79a7", queries = T),
                                          list(plot = myplot, x = "#009e73", queries = T),
                                          list(plot = mutationTypeHistogram, x = "THCA-SA", queries = F)),
                             ncols = 4),
      set.metadata = list(data = setdata,
                          plots = list(list(type="hist", column = "Donors", assign =15),
                                       list(type="text", column = "Site", assign = 8))))
dev.off()
