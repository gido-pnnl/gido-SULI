##matching bNTI data with geochem data - adding rows for replicates##

geochem.data <- read.csv(file.choose()) #load in geochemical characteristics of samples
bNTI <- read.csv(file.choose(), row.names = 1) #load in bNTI values
solenoid <- read.csv(file.choose())

rep.names <- rownames(bNTI) #stores sample names with replicate ID
rep.names <- as.data.frame(as.character(rep.names))
colnames(rep.names) <- "rep.names"

rep.names$Sample = gsub("\\..*","",rep.names$rep.names) #stores sample names without replicate ID for matching

#this step fixes a problem in the FTIRC-MS sample names; some samples only had 5 digits instead of 6 which did not match geochem data
rep.names[117:131,2] <- c("A000001", "A000001", "A000002", "A000003", "A000004", "A000005", "A000006", "A000007", "A000007", "A000008", "A000009", "A000009", "A000009", "A000010", "A000010")

geochem.data <- merge(rep.names, geochem.data, all.x = TRUE, by = "Sample") #matches geochem data with sample replicates -- essentially creates extra rows for sample replicates

#this step fixes another problem with the FTIRC-MS data set: it is organized starting with the 104th sample then repeats back to the 1st sample later.. this is how bNTI output is structured and is all but unavoidable
geochem.data <- geochem.data[c(188:303,1:187),] 
rownames(geochem.data) <- geochem.data$rep.names #resetting row names to avoid confusion with numbers

#this step merges in solenoid data
geochem.data$order <- 1:nrow(geochem.data) 
geochem.data <- merge(solenoid, geochem.data, by = "Solenoid", all.x = FALSE, all.y = TRUE)
geochem.data <- geochem.data[order(geochem.data$order),]
rownames(geochem.data) <- geochem.data$rep.names

#parsing geochemical data to numerical variables
geochem.num <- geochem.data[, c(2,3,11,12,13,22:30,32:42)]
geochem.num <- apply(geochem.num, 2, function(x) gsub("<", NA, x))
geochem.num <- apply(geochem.num, 2, as.numeric)
geochem.num <- as.data.frame(geochem.num)

geochem.diff <- lapply(geochem.num, function(x) abs(outer(x,x,'-'))) #calculates difference matrix for each geochemical variable

##mantel test

mant <- matrix(nrow = 25, ncol = 2) #create storage dummy
rownames(mant) <- colnames(geochem.diff)
colnames(mant) <- c("Mantel Stat", "Significance")

for(i in 1:length(geochem.diff)){ #calculate mantel values for each geochemical variable
  temp <- geochem.diff[[i]]
  m <- mantel(as.dist(temp), as.dist(bNTI), permutations = 999, na.rm = T, method = "spearman")
  mant[i,1] <- m[[3]]
  mant[i,2] <- m[[4]]}

##creating graphs of average bNTI vs geochem values

bNTI <- as.matrix(bNTI) #create full bNTI matrix
bNTI <- reflect_triangle(bNTI, from = "lower")

avg.bNTI <- apply(bNTI, 1, mean, na.rm = T) #calculate averages per samples
avg.bNTI <- as.data.frame(avg.bNTI)
geochem.num$bNTI <- avg.bNTI$avg.bNTI

for(i in 1:25){ #save all plots of avg bNTI vs geochem data
  formula <- geochem.num[,i]~geochem.num$bNTI
  print(ggplot(data = geochem.num, aes(x = geochem.num[,i], y = bNTI)) + geom_point() +
    geom_smooth(method = "lm") +
    stat_poly_eq(formula = formula, 
                 aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
                 parse = TRUE) +
    xlab(colnames(geochem.num)[i]))
  ggsave(paste0("avg_bNTI_", colnames(geochem.num)[i], ".pdf", sep = ""))
}

##creating same graphs but with outliers (avg. bNTI > 0) removed

geochem.num.no.outliers <- geochem.num[-which(geochem.num$bNTI > 0),]

for(i in 1:25){ #save all plots of avg bNTI vs geochem data
  formula <- geochem.num.no.outliers[,i]~geochem.num.no.outliers$bNTI
  print(ggplot(data = geochem.num.no.outliers, aes(x = geochem.num.no.outliers[,i], y = bNTI)) + geom_point() +
          geom_smooth(method = "lm") +
          stat_poly_eq(formula = formula, 
                       aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
                       parse = TRUE) +
          xlab(colnames(geochem.num.no.outliers)[i]))
  ggsave(paste0("avg_bNTI_", colnames(geochem.num.no.outliers)[i], "_no_outliers.pdf", sep = ""))
}

##creating box plots of categorical variables and average bNTI

geochem.data$bNTI <- avg.bNTI$avg.bNTI
geochem.data.no.outliers <- geochem.data[-which(geochem.data$bNTI > 0),]
geochem.data.no.outliers$Solenoid <- as.factor(geochem.data.no.outliers$Solenoid)
geochem.data.no.outliers$Date_time <- gsub("\\ .*", "", geochem.data.no.outliers$Date_time)

ggplot(data = geochem.data.no.outliers, aes(x = Elevation.1, y = bNTI)) + #for terrain
  geom_boxplot() +
  geom_jitter()
ggsave("avg_bNTI_terrain_plot.pdf")

ggplot(data = geochem.data.no.outliers, aes(x = Solenoid, y = bNTI)) + #for solenoid
  geom_boxplot() +
  geom_jitter()
ggsave("avg_bNTI_solenoid_plot.pdf")

ggplot(data = geochem.data.no.outliers, aes(x = Date_time, y = bNTI)) + #for solenoid
  geom_boxplot() +
  geom_jitter()
ggsave("avg_bNTI_date_plot.pdf")
