#Feb 3, 2025 - Katja C. Seltmann

rm(list=ls())
library(tidyverse)
library("ape")
library("phytools")
library('picante')

#### tree and taxon name list comparisons. First read in generic and species trees

#### Henriquez Piskulich, Patricia Andrea; Hugall, Andrew F.; Stuart-Fox, Devi (2023). A supermatrix phylogeny of the world’s bees (Hymenoptera: Anthophila) [Dataset]. Dryad. https://doi.org/10.5061/dryad.80gb5mkw1

#############genus tree##############################
#genus tree
mytreegenus <- read.tree('../data/BEE_mat7gen_p8pmAa_fst.nwk')

# Extract the genus from the tip labels
# Format in data file: "genus_species~~family~~subfamily~tribe"
mytreegenus$tip.label <- sub("_.*", "", mytreegenus$tip.label)
mytreegenus$tip.label <- sub("~.*", "", mytreegenus$tip.label)

#unrooted tree, includes branch lengths, 433 tips and 431 internal nodes
mytreegenus

#tip labels
mytreegenus$tip.label

writeLines(mytreegenus$tip.label, con = "../output/output-genus-tree-default-list.txt")


#############species tree##############################
#species tree
mytreespecies <- read.tree('../data/BEE_mat7_fulltree.nwk')

# Extract the species from the tip labels
#imported format: Andrena_asteris_n4m1_Andrenidae

#format like: "Protandrena solidaginis, Andrenidae"
mytreespecies$tip.label <- sub("^(\\w+)_([^_]+)_n..._(\\w+)$", "\\1 \\2, \\3", mytreespecies$tip.label)

#just species names
mytreespecies$tip.label <- sub("^(\\w+)_([^_]+)_n..._(\\w+)$", "\\1 \\2", mytreespecies$tip.label)

#print to a file
writeLines(mytreespecies$tip.label, con = "../output/output-species-tree-default-list.txt")

#######################compare with occurrence data###########
df <- read.csv("../data/species.csv", header = FALSE, stringsAsFactors = FALSE)

#occurrence data generic names
genus_names <- sub("^(\\S+).*", "\\1", df$V1)

#occurrence data species names
df

#missing in genus tree
missing_in_tree <- setdiff(df$Genus, mytreegenus$tip.label)

# Convert the list to a vector and collapse elements into a single string with new lines
result <- paste(unlist(missing_in_tree), collapse = "\n")

# Print the result so that each element is on a new line
cat(result)

#missing genera
#Xenoglossodes
#Peponapis
#Micralictoides
#Xeroheriades
#Agapanthinus

#missing in genus tree

species <- data.frame(tip.label = mytreespecies$tip.label, stringsAsFactors = FALSE)

#rename column so they are the same in both dataframes
names(df)[names(df) == "V1"] <- "tip.label"

missing_in_tree <- setdiff(df, species)

# Convert the list to a vector and collapse elements into a single string with new lines
result <- paste(unlist(missing_in_tree), collapse = "\n")

row_count <- nrow(missing_in_tree)
print(row_count)


#########read.table()#################OLD CODE BELOW#########################
### read in the globi bee data
globi = read_csv('modeling_data/globi_allNamesUpdated_Henriquez_Piskulich.csv') %>%
  mutate(bee_genus = sub(" .*","",scientificName))
globi_genera = unique(globi$bee_genus)

###root with apoid wasp outgroups
#Pulverro (Ammoplanidae), Bembix (Bembicidae), Cerceris (Philanthidae), Philanthus (Philanthidae) and Tachysphex (Crabronidae) 

outgroups <- c("Philanthus","Pulverro","Bembix","Tachysphex","Cerceris")
workingtree=root(mytree,outgroup=outgroups, resolve.root = TRUE)
workingtree=as.phylo(mytree)

##Make ultrametric
workingtree=chronos(workingtree, lambda = 0.001, control = chronos.control(maxit = 10000000, eval.max = 100000000))
is.ultrametric(workingtree)

plot.phylo(workingtree, cex = .1, main="Henriquez Piskulich et. al")
p <- ggtree(workingtree) + geom_tiplab()
print(p)

# Find the intersection of the extracted genus names and your list of names
tree_tips <- workingtree$tip.label
matching_names <- tree_tips[!tree_tips %in% globi_genera]

# Trim the tree to include only the matching names
trimmed_tree <- drop.tip(workingtree, matching_names)

#check labels
trimmed_tree$tip.label

#check if ultrametric
is.ultrametric(trimmed_tree)

#check branch lengths
branch_lengths <- trimmed_tree$edge.length

##Make ultrametric
#The output of chronos is an ultrametric phylogenetic tree where branch lengths represent time 

#forced ultrametric
#forcedTree<- force.ultrametric(trimmed_tree)
#phylo_dist = cophenetic.phylo(forcedTree)

#1.0 is the default lambda, lacks differences betwen genera
#calculatedTree=chronos(trimmed_tree, lambda = 0.000001)
#plot.phylo(calculatedTree, cex = .5, main="Henriquez Piskulich et. al")


#check if ultrametric
is.ultrametric(trimmed_tree)

#mytree_edited <- read.tree('final_data/modified_tree_Henriquez_Piskulich.nwk')
#plot.phylo(mytree_edited, cex = .5, main="Henriquez Piskulich et. al")

#plot and write new tree
#plot.phylo(calculatedTree, cex = .5, main="Henriquez Piskulich et. al")
#write.tree(trimmed_tree, file="modified_tree_Henriquez_Piskulich.nwk")

#double check tree looks okay for a few genera in each fam
keep = c("Osmia","Hesperapis", 'Bombus',"Andrena","Melitta","Perdita",'Colletes','Megachile','Lasioglossum')
all_gen = new_bee_tree$tip.label
#(rm_genera = all_gen[!all_gen %in% keep])
rm_genera = trimmed_tree$tip.label[!trimmed_tree$tip.label %in% keep]
very_pruned = drop.tip(trimmed_tree,rm_genera)

#plot and write new tree
plot.phylo(very_pruned, cex = .7, main="HP")

#calculate phylogenetic distances
phylo_dist = cophenetic.phylo(trimmed_tree)
my_pcoa <- stats:::cmdscale(phylo_dist,eig=T)
plot(my_pcoa$points) # looks pretty weird

# Print the distance matrix with high precision
print(phylo_dist, digits = 10)

# Check for repeated distances
unique_distances <- unique(as.vector(phylo_dist))
print(unique_distances)

my_pcoa$points
my_pcoa$eig[1:6]
bee_fam = data.frame(bee_genus = row.names(my_pcoa$points),eigen1 = my_pcoa$points[,1],eigen2 = my_pcoa$points[,2]) %>%
  left_join(globi %>% distinct(bee_genus,bee_family))

labels <- data.frame(bee_genus = trimmed_tree$tip.label)

#add pairwise phylo_dist
phylo_df = as.data.frame(phylo_dist)
phylo_df$bee_genus = row.names(phylo_dist)
row.names(phylo_df) <- NULL

 bee_fam2 = labels %>%
   left_join(bee_fam %>% left_join(phylo_df))
 color_pal=RColorBrewer::brewer.pal(8,'Set3')[c(1,3:8)]
 my_cols = adjustcolor(color_pal[as.factor(bee_fam$bee_family)],.4)
 my_cols2 = adjustcolor(color_pal[as.factor(bee_fam$bee_family)],.7)

with(bee_fam,plot(eigen1,eigen2,pch=16,col = my_cols))
legend("bottomleft",unique(bee_fam$bee_family),col=unique(my_cols2),pch=16)


write_csv(bee_fam2 %>% select(bee_genus,everything()),'modeling_data/bee_phylogenetic_data_Henriquez_Piskulich_tree.csv')


# Comparing output trees
Hedtke <- read.tree('modified_tree_Hedtke.nwk')
Hedtke_binary_tree <- multi2di(Hedtke)
Henriquez_Piskulich <- read.tree('modified_tree_Henriquez_Piskulich.nwk')

# Compare Tree Topologies
dend1 <- as.dendrogram(Hedtke_binary_tree)
dend2 <- as.dendrogram(Henriquez_Piskulich)
topology_comparison <- all.equal(dend1, dend2)
cat("Topology Comparison:", topology_comparison, "\n")

# Are the trees rooted?
is.rooted(Hedtke)
is.rooted(Henriquez_Piskulich)

# Differences in tips
tips_Hedtke <- Hedtke$tip.label
tips_Henriquez_Piskulich <- Henriquez_Piskulich$tip.label #has more tips
setdiff(tips_Hedtke, tips_Henriquez_Piskulich)
setdiff(tips_Henriquez_Piskulich, tips_Hedtke)



