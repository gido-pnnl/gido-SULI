setwd("C:\\Users\\gido282\\OneDrive - PNNL\\Documents\\Data") #set working directory

data <- read.csv(file.choose(), row.names = 1) #load in peak data
bad <- read.csv(file.choose(), header = FALSE) #load in peaks that need removed
tree <- read.tree(file.choose()) #load in tree

peak.data <- data[,!(colnames(data) %in% bad[,1])] #removing bad samples from data
w <- which(rowSums(peak.data) == 0) #finding indices where peaks no longer exist
peak.data <- peak.data[-w,] #remove non-existent peaks from peak data

peak.data[peak.data > 1] = 1 #turns peak data into present/absent

library(picante)
library(reshape2)
library(abind)

tree_type <- "broad_tree"
Sample_Name <- "300A_site"

acomb <- function(...) abind(..., along = 3) #binds matrices along 3 dimensions

match <- match.phylo.data(tree, peak.data)

coph <- cophenetic(match$phy)

bMNTD <- as.matrix(comdistnt(t(match$data), coph, abundance.weighted = F, exclude.conspecifics = F)) #calculates experimental bMNTD values

files <- list.files(path = paste(tree_type, "_Null_Results/", sep = "") #lists all nulls
                   , pattern = "bMNTD_rep", full.names = TRUE) 

rand.bMNTD = NULL # Dummy object

for(curr.file in files){ #merges nulls into a list 
  temp = as.data.frame(read.csv(curr.file, row.names = 1))
  rand.bMNTD = c(rand.bMNTD, list(temp))
} 

rand.bMNTD = do.call(acomb, rand.bMNTD) #merges nulls in list into 3D array
rm("curr.file")

bNTI = matrix(c(NA), nrow = ncol(rand.bMNTD[,,1]), ncol = ncol(rand.bMNTD[,,1])) #creates empty matrix

for(i in 1:(ncol(rand.bMNTD[,,1])-1)){
  for(j in (i+1):ncol(rand.bMNTD[,,1])){
    m = rand.bMNTD[j,i,] # Just setting all the randomizations for a given comparison to a matrix
    bNTI[j,i] = ((bMNTD[j,i]-mean(m))/sd(m)) # The bNTI calculation
  }
}

dimnames(bNTI) = dimnames(rand.bMNTD[,,1]) #renaming dimensions 
rm("m", "j")

write.csv(bNTI, paste(Sample_Name, "_", tree_type, "_bNTI_", length(files), ".csv", sep = ""), quote = F)
