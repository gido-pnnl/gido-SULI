bNTI <- read.csv(file.choose(), row.names = 1) #load in bNTI values

list.bNTI <- c(NULL) #create dummy list

for(i in 1:ncol(bNTI)){ #list all bNTI values as one column
  list.bNTI <- c(list.bNTI, bNTI[,i])
}

list.bNTI <- as.data.frame(list.bNTI) #store as data frame
list.bNTI <- na.omit(list.bNTI) #omit NAs

mean.bNTI <- mean(list.bNTI$list.bNTI) #mean

r <- ggplot(data = list.bNTI, aes(x = list.bNTI)) + #create histogram with density curve
  geom_histogram(aes(y = ..density..), bins = 200, color = "red", fill = "white") +
  geom_density(alpha = 0.5, size = 1) +
  geom_vline(xintercept = mean(list.bNTI$list.bNTI), size = 1) +
  xlab("bNTI") +
  ylab("Density") +
  ggtitle("bNTI Distribution for 300A Site", )
r
