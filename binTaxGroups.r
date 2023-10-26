# get abundances by taxon rank.

get_tax_level_abun <- function(ps, tax_rank_list = c("phylum","class","order","family","genus"), min_seq_depth = 1000) {
require("phyloseq")
	require("speedyseq")
	require("dplyr")
	require("data.table")
  
  if(taxa_are_rows(ps)){otu_table(ps) <- t(otu_table(ps))}
	ps_orig <- prune_samples(rowSums(otu_table(ps)) > min_seq_depth, ps)
	ps_orig <- prune_samples(!duplicated(sample_data(ps_orig)$sampleID), ps_orig)
	out.list <- list()
for (r in 1:length(tax_rank_list)) {
  ps <- ps_orig
  tax_rank <- tax_rank_list[[r]]
  cat(paste("\nEvaluating at rank:", tax_rank))
  # get sequence total
  seq_total <- rowSums(otu_table(ps))
  
  if (!tax_rank %in% c("phylum","class","order","family","genus","Phylum","Class","Order","Family","Genus")) {
    tax_table(ps) <-  tax_table(ps)[,tax_rank]
  }
  
  glom <- speedyseq::tax_glom(ps, taxrank=tax_rank)
  glom_melt <- speedyseq::psmelt(glom)
  
  if("dnaSampleID" %in% colnames(glom_melt)){
  	form <- as.formula(paste0("dnaSampleID ~ ", tax_rank))
  	glom_wide <- reshape2::dcast(glom_melt, form, value.var = "Abundance", fun.aggregate = sum)
  	out_abun <- transform(glom_wide, row.names=dnaSampleID, dnaSampleID=NULL)
  	out_abun <- out_abun[sample_data(ps)$dnaSampleID,]
  	
  } else {
form <- as.formula(paste0("sampleID ~ ", tax_rank))
glom_wide <- reshape2::dcast(glom_melt, form, value.var = "Abundance", fun.aggregate = sum)
out_abun <- transform(glom_wide, row.names=sampleID, sampleID=NULL)
out_abun <- out_abun[sample_data(ps)$sampleID,]
}
out_abun$other <- NULL
out_abun$other <- seq_total - rowSums(out_abun)

# turn into relative abundances
out_rel <- out_abun/rowSums(out_abun)

# Make a prevalence (frequency) table too
# Compute prevalence of each feature, store as data.frame
prevdf <- apply(X = otu_table(glom), 2, function(x){sum(x > 0)})
N.SVs <- data.frame(table(tax_table(ps)[, tax_rank], exclude = NULL))
colnames(N.SVs)[2] <- "N.SVs"
# Add taxonomy and total read counts to this data.frame
prevdf <- data.frame(prevalence = prevdf/nsamples(glom),
                    totalAbundance = taxa_sums(glom),
                    tax_table(glom)[,tax_rank])
prevdf <- merge(prevdf, N.SVs, by.x = colnames(prevdf)[3], by.y  = "Var1", all.x=T)

out <- list(out_abun, out_rel, seq_total, prevdf)
names(out) <- c("abundances","rel.abundances","seq.total","prevalence")
out.list[[r]] <- out
}
names(out.list) <- tax_rank_list
return(out.list)
}
# 
# setwd("/projectnb/dietzelab/zrwerbin/NEON_soil_microbe_processing")
# library(phyloseq)
# library(speedyseq)
# ps <- readRDS("./data/NEON_16S_phyloseq.rds")
# ps_prune <- phyloseq::filter_taxa(ps, function(x) sum(x > 0) > 3, TRUE)
# saveRDS(ps_prune, "./data/NEON_16S_phyloseq_subset.rds")
# 
# 
# sample_data(ps_prune)$sampleID <- sample_names(ps_prune)
# seq_total <- rowSums(otu_table(ps_prune))
# rel <- ps %>% transform_sample_counts(~ . / seq_total)
# 
# out <- get_tax_level_abun(ps_prune, tax_rank_list = c("Phylum","Class","Order","Family","Genus"), min_seq_depth = 5000)
# saveRDS(out, "./data/abundances_16S")


# #### 16S #####
# 
# library(dplyr)
# recent_ps <- readRDS("/projectnb/dietzelab/zrwerbin/NEON_soil_microbe_processing/data/NEON_16S_phyloseq_subset.rds")
# legacy_ps <- readRDS("/projectnb/dietzelab/zrwerbin/NEON_soil_microbe_processing/data/NEON_16S_phyloseq_legacy.rds")
# source("/projectnb/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/helperFunctions.r")
# new_sample_dat <- parseNEONsampleIDs(as.character(sample_data(recent_ps)$dnaSampleID))
# rownames(new_sample_dat) <- rownames(sample_data(recent_ps))
# sample_data(recent_ps) <- new_sample_dat
# 
# master_ps <- merge_phyloseq(legacy_ps, recent_ps)
# colnames(tax_table(master_ps)) <- tolower(colnames(tax_table(master_ps)))
# 
# out <- get_tax_level_abun(master_ps, tax_rank_list = c("Phylum","Class","Order","Family","Genus"), min_seq_depth = 5000)
# saveRDS(out, "/projectnb/dietzelab/zrwerbin/NEON_soil_microbe_processing/data/master_abundances_16S")

