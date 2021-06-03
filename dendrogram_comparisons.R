# ##################################### #
# REMOVING BAD PEAKS AND PREPARING DATA #
# ##################################### #

setwd("C:\\Users\\gido282\\OneDrive - PNNL\\Documents\\Data") #set working directory

data <- read.csv(file.choose(), row.names = 1) #load in peak data
mol <- read.csv(file.choose(), row.names = 1) #load in molecular characteristics
bad <- read.csv(file.choose(), header = FALSE) #load in peaks that need removed

peak.data <- data[,!(colnames(data) %in% bad[,1])] #removing bad samples from data

w <- which(rowSums(peak.data) == 0) #finding indices where peaks no longer exist

peak.data <- peak.data[-w,] #removing non-existent peaks from data frame

mol.data <- mol[-w,] #removing corresponding peaks from mol. data frame

#check to see if row names match
if(!identical(rownames(mol.data), rownames(peak.data))){
  stop("the row names between mol_data and peak_data don't match")
} else {
  print("row names of mol_data and peak_data match")
}

#transposing peak data - necessary for pd function
peak.data.t <- t(peak.data)

# ############## #
# CREATING TREES #
# ############## #

tree1 <- mol.data[-which(mol.data$MolForm %in% NA),] #removing peaks with no formula assignment
tree1 <- tree1[,c("C", "H", "O", "N", "S", "P"), drop = F] #choose variables for clustering
tree1 <- as.data.frame(apply(tree1, 2, scale), row.names = row.names(tree1)) #scaling values
tree1 <- as.phylo(hclust(vegdist(tree1, "euclidean")), "average") #creating distance matrix -> storing clusters as tree

tree2 <- mol.data[-which(mol.data$MolForm %in% NA),]
tree2 <- tree2[,c("OtoC_ratio", "HtoC_ratio", "NtoC_ratio", "PtoC_ratio"), drop = F] #choose different variables for clustering
tree2 <- as.data.frame(apply(tree2, 2, scale), row.names = row.names(tree2))
tree2 <- as.phylo(hclust(vegdist(tree2, "euclidean"), "average"))

# ################################## #
# CALCULATING PHYLOGENETIC DIVERSITY #
# ################################## #

tree.results <- matrix(NA, nrow = nrow(peak.data.t), ncol = 2) #making empty matrix for pd results
colnames(tree.results) <- c("Tree_1_PD", "Tree_2_PD") #naming columns 
rownames(tree.results) <- rownames(peak.data.t) #naming rows

for(i in 1:nrow(peak.data.t)){
  tree.results[i,1] <- pd(t(peak.data.t[i,]), tree1)[,1] #stores PD for each sample in respective column of tree.results
  tree.results[i,2] <- pd(t(peak.data.t[i,]), tree2)[,1]
}

tree.results <- as.data.frame(tree.results) #storing results as data frame

# ############## #
# VISUALIZE DATA #
# ############## #

ggplot(tree.results, aes(x = Tree_1_PD, y = Tree_2_PD)) + #sets tree 1 PD on x-axis and tree 2 PD on y-axis 
  geom_point() + #scatter plot
  geom_abline(slope = 1) + #y = x reference line
  coord_cartesian(xlim = c(0,NA), ylim = c(0,NA)) + #graph starts at (0,0)
  xlab("Phylogenetic Diversity - Number of Atoms") + #x-axis (tree 1) label
  ylab("Phylogenetic Diversity - Atomic Ratios") #y-axis (tree 2) label

#done :D