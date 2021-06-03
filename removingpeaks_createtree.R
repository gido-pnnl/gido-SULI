# ################## #
# REMOVING BAD PEAKS #
# ################## #

setwd("C:\\Users\\gido282\\OneDrive - PNNL\\Documents\\Data") #set working directory

data <- read.csv(file.choose(), row.names = 1) #load in peak data
mol <- read.csv(file.choose(), row.names = 1) #load in molecular characteristics
bad <- read.csv(file.choose(), header = FALSE) #load in peaks that need removed

data2 <- data[,!(colnames(data) %in% bad[,1])] #removing bad samples from data

w <- which(rowSums(data2) == 0) #finding indices where peaks no longer exist

peak_data <- data2[-w,] #removing non-existant peaks from data frame

mol_data <- mol[-w,] #removing corresponding peaks from mol. dataframe


#check to see if row names match
if(!identical(rownames(mol_data), rownames(peak_data))){
  stop("the row names between mol_data and peak_data don't match")
} else {
  print("row names of mol_data and peak_data match")
}

# #### #
# TREE #
# #### #

mol.info <- mol_data[-which(mol_data$MolForm %in% NA),] #removing peaks with no formula assignment
mol.info <- mol.info[,c("C", "H", "O", "N", "S", "P", "DBE", "AI_Mod", "NeutralMass"), drop = F] #choosing variables for clustering
mol.info <- as.data.frame(apply(mol.info, 2, scale), row.names = row.names(mol.info)) #scaling values
tree <- as.phylo(hclust(vegdist(mol.info, "euclidean")), "average") #creating distance matrix -> storing clusters as tree

mol_data <- cbind(rownames(mol_data), mol_data) #for tree aesthetics to work

#visualizing tree
col = colorRampPalette(c("dodgerblue4", "goldenrod3", "firebrick3"))(length(unique(mol_data$bs1_class)))
the_tree <- ggtree(tree, layout = "circular") %<+% mol_data + geom_tippoint(aes(color = bs1_class), na.rm = T) +
  scale_color_manual(values = col)+
  theme(legend.position = "right")
