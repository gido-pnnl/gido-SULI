##preparing data

setwd("C:\\Users\\gido282\\OneDrive - PNNL\\Documents\\Data") #set working directory

data <- read.csv(file.choose(), row.names = 1) #load in peak data
bad <- read.csv(file.choose(), header = FALSE) #load in peaks that need removed
tree <- read.tree(file.choose()) #load in tree

peak.data <- data[,!(colnames(data) %in% bad[,1])] #removing bad samples from data
w <- which(rowSums(peak.data) == 0) #finding indices where peaks no longer exist
peak.data <- peak.data[-w,] #remove non-existent peaks from peak data

peak.data[peak.data > 1] = 1 #turns peak data into present/absent

tree_type <- "tree_type" #name of tree
Sample_Name <- "sample_name" #name of data set

if(!dir.exists(paste0(tree_type, "_Null_Results"))){ #creates directory for null bMNTD matrix storage
  dir.create(paste0(tree_type, "_Null_Results"))
}

##bMNTD process

matched <- match.phylo.data(tree, peak.data) #matches tree to peak data, stores as list with "data" and "phy"
coph <- cophenetic(matched$phy) #creates cophenetic distance matrix relating species

range = 1:999 #number of iterations

for(i in range){ #in this iteration, data is structured with peaks as rows and samples as columns - returns distances between samples
  bMNTD.rand = as.matrix(comdistnt(t(matched$data), taxaShuffle(coph), abundance.weighted = F, exclude.conspecifics = F)) #calculates bMNTD relating samples with randomized distance matrix
  write.csv(bMNTD.rand, paste(tree_type, "_Null_Results/FTICR_", Sample_Name, "_", tree_type, "_bMNTD_rep", i, ".csv", sep = ""), quote = F) #writes matrix as .csv
  rm("bMNTD.rand")
  
  print(c(date(),i)) #prints date/time
}
