#A function to generate SNP heat maps
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
                                                                          plot.title = element_text(color=colour, size = 10),
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

#Call the upset function and apply intersection queries
#Then use the info from the intersection queries to plot the SNP heat maps as attribute plots
upset(myData, nsets = 6,
      sets.x.label = "Mutations", mainbar.y.label = "Shared Mutations", mb.ratio = c(0.5,0.5),
      order.by = "freq", nintersects = 30,
      queries = list(list(query = intersects, params = list("LUSC-US"), color = "#E69F00"),
                     list(query = intersects, params = list("LUSC-KR"), color = "#56B4E9"),
                     list(query = intersects, params = list("LAML-KR"), color = "#009E73"),
                     list(query = intersects, params = list("THCA-SA"), color = "#F0E442"),
                     list(query = intersects, params = list("THCA-SA", "LAML-KR"), color = "#0072B2"),
                     list(query = intersects, params = list("LUSC-KR", "THCA-SA"), color = "#D55E00")),
      attribute.plots = list(gridrows = 80,
                             plots = list(list(plot = myplot, x = "#E69F00", queries = T),
                                          list(plot = myplot, x = "#56B4E9", queries = T),
                                          list(plot = myplot, x = "#009E73", queries = T),
                                          list(plot = myplot, x = "#F0E442", queries = T),
                                          list(plot = myplot, x = "#0072B2", queries = T),
                                          list(plot = myplot, x = "#D55E00", queries = T)), ncols = 3),
      set.metadata = list(data = setdata,
                          plots = list(list(type="hist", column = "Donors", assign =15),
                                       list(type="text", column = "Site", assign = 8))))