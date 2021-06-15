##matching bNTI data with geochem data - adding rows for replicates##

geochem.data <- read.csv(file.choose()) #load in geochemical characteristics of samples
bNTI <- read.csv(file.choose(), row.names = 1) #load in bNTI values

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

##pairwise comparison of sample geochemical characteristics##

anions <- as.data.frame(geochem.data$Anions_eq.per.L) #find sample values
rownames(anions) <- rownames(geochem.data)

anion.diff <- matrix(0, nrow(anions), nrow(anions)) #make dummy matrix

for(i in 1:nrow(anions)){ #comparing each value to all other values
  anion.diff[,i] <- anions[i,1] - anions[,1]
}

anion.diff <- as.data.frame(anion.diff) #organizing matrix
rownames(anion.diff) <- rownames(anions)
colnames(anion.diff) <- rownames(anions)
anion.diff <- abs(anion.diff)

plot(anion.diff[lower.tri(anion.diff)], bNTI[lower.tri(bNTI)]) #plot bNTI values vs sample values