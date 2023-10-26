library("knitr")
library("rmarkdown")
library("ggplot2")
library("reshape2")
library("ggrepel")
library("RColorBrewer")
library("viridis")
library("vegan")
library("betapart") ##
library("multcomp")
library("tidyr")
library("broom")
library("dplyr")
library("ztable") ##
library("cowplot")
library("VennDiagram")
library("ggsci")
library("gplots")
library("eulerr") ##
library("agricolae")
library("labdsv")
library("colorspace")
library('pairwiseAdonis') ## PACKAGE NOT AVAILABLE
library('phyloseq')
library('energy')
library('colorspace')
library('scales')
library('spaa') ##
library('igraph')

setwd("/projectnb2/talbot-lab-data/zrwerbin/Naylor_functional_groups/FunctionalModule_DataFiles_forGitHub")

no_meta.16s <- import_biom("OTU.biom", "OTU.tree")

metadatatable.16s <- import_qiime_sample_data("OTU_metadata.txt")
metadatatable.16s$X <- NULL # The import function gives an extra column that should be rownames.
# The above function doesn't give you the option to specify rownames, so you have to edit those out manually.
# Merge these to create one phyloseq object.
unrar.16s <- merge_phyloseq(no_meta.16s, metadatatable.16s)
###########################################################################################
### Making some minor modifications to this dataset to exclude some extraneous samples. ###
###########################################################################################
# 'Early benomyl' samples were extraneous to this dataset. We only need the late samples.
unrar.16s <- subset_samples(unrar.16s, !grepl("BEN_EARLY", rownames(sample_data(unrar.16s))))
# 'Tetracycline' had too low of sample counts to be usable. Our rule was if any samples have 3 or more
# replicates below the rarefaction threshold, we'd get rid of them, as we want at least 3 of the 5 total
# replicates to be usable in the downstream analysis.
unrar.16s <- subset_samples(unrar.16s, !grepl("tetracycline", rownames(sample_data(unrar.16s))))
# Similarly, 'glyoxylate' had 3 samples below the rarefaction threshold, so I'm getting rid of them here.
unrar.16s <- subset_samples(unrar.16s, !grepl("glyoxylate", rownames(sample_data(unrar.16s))))
#################################################################################################
### Subsetting to include up to 5 replicates with the highest values for Shannon's diversity. ###
#################################################################################################
# Note: for a few modules, we prepared 6 or more replicates to have extras just in case, as we 
# wanted to have at least 5 replicates for as many modules as possible, and it might have
# transpired that some replicates would have been excluded due to poor quality or too low of 
# readcounts per sample. What this means is that we have uneven replicate numbers per module -
# so, we are going to make sure we have a maximum of 5 replicates for each module, choosing only
# the 5 replicates that have the highest values for Shannon's diversity.
# Make a dummy phyloseq object - we have to rarefy for Shannon's diversity to work correctly.
testrar.16s <- rarefy_even_depth(unrar.16s, rngseed = 711, sample.size = 8000)
df.topdiv <- cbind(dplyr::select(as(sample_data(testrar.16s), "data.frame"), ModuleName, ModuleNameRep),
                   Shannons = estimate_richness(testrar.16s, measures = c("Shannon")))
df.topdiv <- data.table::data.table(df.topdiv, key = "ModuleName", keep.rownames = TRUE)
# Order by ModuleName then by Shannon in descending order.
df.topdiv <- df.topdiv[order(ModuleName, -Shannon),]
# Create a dataframe with just the top 5 values for Shannon by ModuleName.
df.top5div <- df.topdiv[, head(.SD, 5), by = ModuleName]
# Testing to see how many modules are left with fewer than 5 reps.
table(df.top5div$ModuleName) %>% as.data.frame() %>% subset(., Freq < 5)
# ^ There are 4 modules with 4 reps and one with 3. Given that this is based on the rarefied
# dataset, this indicates that after rarefaction the vast majority of our modules should still have
# 5 replicates.
# Now subset the phyloseq object so it just has the top 5 replicates for Shannon's diversity.
unrar.16s <- subset_samples(unrar.16s, ModuleNameRep %in% df.top5div$ModuleNameRep)
######################################################
### Processing the sample data for downstream use. ###
######################################################
# Converting all variables to character variables (even Replicate).
sample_data(unrar.16s)[, ] <- lapply(sample_data(unrar.16s)[, ], as.character)
# Adding in factor levels.
sample_data(unrar.16s)$Replicate <- factor(sample_data(unrar.16s)$Replicate, levels = sort(unique(sample_data(unrar.16s)$Replicate)))
sample_data(unrar.16s)$Category <- factor(sample_data(unrar.16s)$Category, levels = c("Soil", "Simple Substrates", "Antibiotics", "Polysaccharides",
                                                                                      "Anaerobic", "Stresses"))
sample_data(unrar.16s)$Stress <- factor(sample_data(unrar.16s)$Stress, levels = c("control", "2,4-D", "late", "PEG", "pH6", "pH8", "light",
                                                                                  "10%C", "heat", "salt"))
# Adding in taxonomic levels.
colnames(tax_table(unrar.16s)) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
# Rarefy samples.
rardepth <- 8000 # Note: we'll be referring back to this number later, so keep it in mind.
rar.16s <- rarefy_even_depth(unrar.16s, rngseed = 711, sample.size = rardepth)
rar.16s # <- This is the full dataset. We're rarefying to 8000, which will end up eliminating some reps.
## We are left with a phyloseq data object encompassing 324 samples and 2542 taxa.



###################################################################################################
### Function that takes a phyloseq object, subsets it down to a module, then calculates a core. ###
###################################################################################################
# The criteria we are using for our core that a taxa is present in at least 40% of replicates (at least 2 reps),
# and has > 0.01% of cumulative readcounts for that core.

get_core_otus.table <- function(physeq, modulename, pct = 0.4){
  
  # 1. Get your subset.
  modulename <- as.character(modulename)
  subset.16s <- prune_samples(sample_data(physeq)$ModuleName == modulename, physeq)
  
  # 2. Make your OTU dataframe.
  subset.otus <- data.frame(otu_table(subset.16s))
  sample_data(subset.16s)$ModuleRep <- with(sample_data(subset.16s), interaction(ModuleName, Replicate))
  colnames(subset.otus) <- sample_data(subset.16s)$ModuleRep
  subset.otus <- subset.otus[rowSums(subset.otus)>0,]
  
  # 3. Get the full value for full readcount abundance for all OTUs in the subset.otus dataframe.
  cumulative <- sum(colSums(subset.otus))
  
  # 4. Subset the dataframe to just those present in at least [amount you specify as pct] % of samples.
  cutoff <- ncol(subset.otus) * pct
  
  notzero <- data.frame(apply(subset.otus, 1, function(c)sum(c!=0)))
  greaterthanpct <- data.frame(apply(notzero, 1, function(c)(c>=cutoff))) 
  
  colnames(greaterthanpct) <- "Persistence"
  greaterthanpct$Rownames <- rownames(greaterthanpct)
  greaterthanpct <- greaterthanpct[(greaterthanpct$Persistence == TRUE),]
  subset.otus <- subset(subset.otus, rownames(subset.otus) %in% rownames(greaterthanpct))
  
  # 5. Subset the dataframe to those greater than 0.01% of the cumulative abundance.
  core <- subset.otus[rowSums(subset.otus) > (cumulative * 0.0001),]
  core$sums <- rowSums(core)
  
  # 6. Get the taxonomic information, including Phylum, Class, Genus, and OTU.
  core$Phylum <- tax_table(subset.16s)[rownames(tax_table(subset.16s)) %in% rownames(core),
                                       colnames(tax_table(subset.16s)) == "Phylum"] %>% gsub("p__", "", .)
  
  core$Class <- tax_table(subset.16s)[rownames(tax_table(subset.16s)) %in% rownames(core),
                                      colnames(tax_table(subset.16s)) == "Class"] %>% gsub("c__", "", .)
  
  core$Genus <- tax_table(subset.16s)[rownames(tax_table(subset.16s)) %in% rownames(core),
                                      colnames(tax_table(subset.16s)) == "Genus"] %>% gsub("g__", "", .)
  core$OTU <- rownames(core)
  
  # Create one column that has all of the taxonomic information combined.
  core$Taxa <- with(core, interaction(OTU, Phylum, Class, Genus, sep = "."))
  
  return(core)
}
########################################################################################################
### Function that takes a core table and reduces it down to unique taxa, each having a count of '1'. ###
########################################################################################################
# We will need this for consolidating all of the cores together.
get_core_counts.table <- function(table, name){
  # original
  #x <- table %>% group_by(Taxa) %>% summarize(count = n()) %>% arrange(-count) %>% data.frame()
  # keeping genus cols
  x <- table %>% group_by(Taxa) %>% mutate(count = n()) %>% arrange(-count) %>% 
    distinct(Taxa, .keep_all=T) %>% dplyr::select(Taxa, count, Phylum, Class, Genus, OTU) %>% data.frame()
  x$Taxa <- as.character(x$Taxa)
  x$ModuleName <- name
  return(x)
}





# Create a list of all the modules included in this analysis.
module_list <- unique(sample_data(rar.16s)$ModuleName)
full_metadata <- data.frame(sample_data(rar.16s))
# Initialize a dataframe to which we will add the core information.
stepwise.core <- as.data.frame(matrix(nrow = 0, ncol = 6))
colnames(stepwise.core) <- c("Taxa", "OTU", "Phylum", "Genus", "count", "ModuleName")
#######################################################
### Generate a table encompassing all of the cores. ###
#######################################################
# Note: this takes a couple minutes to run, longer depending on how complex the dataset is.
for(module_name in module_list){
  print(module_name) # Tracking how far we've gone.
  module_core <- get_core_otus.table(physeq = rar.16s, modulename = module_name, pct = 0.4)
  module_table <- get_core_counts.table(table = module_core, name = module_name)
  stepwise.core <- rbind(stepwise.core, module_table)
}
stepwise.core.duplicate <- stepwise.core
##########################################################
### Add in the appropriate metadata to the core table. ###
##########################################################
full.cat.mod <- dplyr::select(full_metadata, ModuleName, Category) # Making this dataframe as an easy way to reassign 'Category' in our core table.
full.cat.mod <- full.cat.mod[!duplicated(full.cat.mod),]
full.cat.mod$Category <- as.character(full.cat.mod$Category)
for(i in 1:nrow(stepwise.core)){
  cat <- as.character(full.cat.mod[full.cat.mod$ModuleName == stepwise.core$ModuleName[i], colnames(full.cat.mod) == "Category"])
  stepwise.core$Category[i] <- cat
}


classification_df <- stepwise.core %>% distinct(Class, Genus, ModuleName, Category)

saveRDS(classification_df, "/projectnb/talbot-lab-data/zrwerbin/Naylor_functional_groups/functional_module_df.r")

# Extract names.
module_names <- unique(classification_df$ModuleName) 
categories <- unique(classification_df$Category) 
columns <- c(categories, module_names)

# Create empty matrix for classifications

# Read in phyloseq data 
ps.list <- readRDS("/projectnb/talbot-lab-data/zrwerbin/decomposition/ps_allstudy_silva.rds")
ps_master <- merge_phyloseq(ps.list[[1]], ps.list[[2]], ps.list[[3]], ps.list[[4]], ps.list[[5]], ps.list[[6]])
tax <- as(tax_table(ps_master), "matrix")
tax_df <- as.data.frame(tax)
tax_df[,categories] <- "other"
tax_df[,module_names] <- "other"

# Loop through and assign all taxa to categories or modules
for (c in categories){
  # Subset to taxa in category
  cat_df <- classification_df %>% filter(Category == !!c)
  cat_modules <- cat_df %>% distinct(ModuleName) %>% unlist()
  # Assign a 1 if the genus is in this category of modules
  tax_df[tax_df$Genus %in% cat_df$Genus,c] <- c
  
  for (m in cat_modules){
    # Subset to taxa in module
    module_df <- cat_df %>% filter(ModuleName == !!m)
    # Assign a 1 if the genus is in this specific modules
    tax_df[tax_df$Genus %in% module_df$Genus,m] <- m
  }
}
tax_table(ps_master) <- as(tax_df, "matrix")
colnames(tax_table(ps_master)) <- make.names(colnames(tax_table(ps_master)))

source( "/projectnb2/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/binTaxGroups.r")
colnames(sample_data(ps_master))[[1]] <- "sampleID"
category_abuns <- get_tax_level_abun(ps_master, tax_rank_list = make.names(categories))
module_abuns <- get_tax_level_abun(ps_master, tax_rank_list = make.names(module_names))

saveRDS(list(categories = category_abuns, modules = module_abuns), "/projectnb/talbot-lab-data/zrwerbin/decomposition/functional_module_abundances.rds")


# combine and visualize!
cat_module_ref <- classification_df %>% distinct(Category, ModuleName) %>% mutate(ModuleName = make.names(ModuleName), Category = make.names(Category))

combined <- do.call(cbind, unname(lapply(category_abuns, "[[", 2)))
combined$sampleID <- rownames(combined)
combined <- combined %>% dplyr::select(-ends_with("other"))
df <- data.frame(sample_data(ps_master))
category_to_plot <- merge(df, combined, by = "sampleID")


combined <- do.call(cbind, unname(lapply(module_abuns, "[[", 2)))
combined$sampleID <- rownames(combined)
combined <- combined %>% dplyr::select(-ends_with("other"))
df <- data.frame(sample_data(ps_master))
modules_to_plot <- merge(df, combined, by = "sampleID")
modules_long <- modules_to_plot %>% pivot_longer(cols = 8:73, names_to="module")
modules_long$category <- cat_module_ref[match(modules_long$module, cat_module_ref$ModuleName),]$Category
ggplot(modules_long[modules_long$category=="Simple.Substrates",], aes(x = pct_mass_remaining, y = value, color = study)) + geom_point() + scale_x_reverse() + ylab("xylan") + geom_smooth() + facet_grid(~module)

ggplot(modules_long[modules_long$category=="Polysaccharides",], aes(x = pct_mass_remaining, y = value, color = study)) + geom_point() + scale_x_reverse() + ylab("xylan") + geom_smooth() + facet_grid(~module)

ggplot(modules_long[modules_long$category=="Polysaccharides",], aes(x = pct_mass_remaining, y = value, color = plant)) + geom_point() + scale_x_reverse() + ylab("xylan") + geom_smooth() + facet_grid(~module)


