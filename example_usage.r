library(phyloseq)
library(dplyr)
source("./helperFunctions.r")
source("./createTaxFunction.r")
source("./addBacterialFunction.r")
source("./binTaxGroups.r")

# RECENT DATA #
# Load combined sequence table and taxonomic table
seqtab_orig <- readRDS("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/raw/MCC_otu_16S.rds")
taxa_orig <- read.csv("/projectnb/microbiome/zrwerbin/NEON_amplicon/16S/NEON_16S_taxonomy/final_tax_table.csv", header=F, sep="\t")

# Prep for phyloseq
new_tax <- do.call(rbind, (lapply(taxa_orig$V2, parse_taxonomy_qiime)))
rownames(new_tax) <- taxa_orig$V1
colnames(new_tax)[1:5] <- c("Kingdom","Phylum","Class","Order","Family")
sample_dat <- parseNEONsampleIDs(rownames(seqtab_orig))

ps <- phyloseq(otu_table(seqtab_orig, taxa_are_rows = F),
							 tax_table(new_tax), sample_data(sample_dat))


# Assign functional groups
tax_ref <- createTaxFunction(ref.path = "./reference_data/bacteria_func_groups.csv",
														 N.path = "./reference_data/Npathways_Albright2018.csv",
														 C.path = "./reference_data/cellulolytic_Berlemont.csv",
														 Naylor.path = "./reference_data/reference_data/functional_module_df.rds",
														 out.path = "./tax_function_ref.csv")

tax_df = as.data.frame(as(tax_table(ps), "matrix"))
new_tax <- addBacterialFunction(tax = tax_df[,!colnames(tax_df) %in% "Kingdom"],
																tax_fun_ref = tax_ref)
tax_table(ps) <- new_tax


# store the DNA sequences of our ASVs in the refseq slot of the phyloseq object,
# and then rename our taxa to a short string
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))



# Reshape into abundance and summary dataframes, drop low-quality samples
# Takes a few minutes to run.
out <- get_tax_level_abun(ps,
													tax_rank_list = colnames(tax_table(ps)),
													min_seq_depth = 5000)

# Now go through ranks to get top abundances


n.taxa <- 10

rank.df.bac <- list()
val.out.bac <- list()
for (tax_rank in names(out)){
	print(tax_rank)
	rank_abun <- out[[tax_rank]]$rel.abundances
	prev_top <- out[[tax_rank]]$prevalence[order(out[[tax_rank]]$prevalence$prevalence, decreasing = T),]
	most_abundant_taxa <- prev_top[,1]
	most_abundant_taxa <- most_abundant_taxa[!grepl("unassigned|other|genus|family|order|^bacteria",most_abundant_taxa)][1:n.taxa]
	most_abundant_taxa <- gsub("\\-|\\(|\\)| ", "\\.", most_abundant_taxa)
	out_top10 <- rank_abun[,colnames(rank_abun) %in% most_abundant_taxa, drop=F]
	seqDepth <- rowSums(rank_abun)
	out_top10$other <- 1-rowSums(out_top10)
	# Remove samples with a sum of above one (not sure why they exist)
	out_top10 <- out_top10[which(!rowSums(out_top10) > 1),]


	ps.rank.filt <- prune_samples(sample_names(ps) %in% rownames(out_top10), ps)
	rank.df <- cbind(sample_data(ps.rank.filt)[,c("siteID","plotID","dateID","sampleID","dates","plot_date")], out_top10)

	# organize by date
	rank.df$dates <- as.Date(as.character(rank.df$dates), "%Y%m%d")
	rank.df <- rank.df[order(rank.df$dates),]

	rank.df.bac[[tax_rank]]	<- rank.df
}


saveRDS(rank.df.bac, "./groupAbundances.rds")



