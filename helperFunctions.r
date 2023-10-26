# helper functions
#if (!require(dada2)) stop("Please install dada2 according to the directions at: https://benjjneb.github.io/dada2/dada-installation.html")
require(ShortRead)
require(Biostrings)
require(RCurl)
require(phyloseq)


# from the dada2 tutorial
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna),
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}


# check for primers
checkPrimers <- function(fwd, rev, FWD.orients, REV.orients){
  print(
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fwd),
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = rev),
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fwd),
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = rev))
)
}


# tic/toc functions from Colin Averill:
#' Two clock functions.
#' Place tic() at the line in the code where you want to start timing.
#' Place toc() at the position in the code where you want to stop timing and report.
tic <- function() {assign("timer", Sys.time(), envir=.GlobalEnv)}
toc <- function() print(Sys.time()-timer)




#normalize otu table. function from Colin Averill:
pro.function <- function(otu){
  for(i in 1:ncol(otu)){
    otu[,i] <- otu[,i] / sum(otu[,i])
  }
  return(otu)
}

#' Set truncation length
#'
#' Decides on truncation length for trimmed reads based on quality score means. Default cutoff is a score of 30.
#' Warns you if the beginning (first 10) bases are low-quality, but returns the first low-quality base after the 10th base.
#' If no bases are below score, returns the last base.
#'
#'  fl : input files (prints suggested length for each; if used in script, could take the first input, or the lowest)
#'  qscore : default = 30.
#'  n : default = 5e+05 don't know exactly why this matters, but it's part of how 'qa' is called within the dada2 scripts...
#'
#'
get_truncation_length <-function (fl, qscore = 30, n = 5e+05, quiet = TRUE){
  trunc_lengths <- data.frame(file = character(0), early_lowqual = numeric(0),
                              trunc_length = numeric(0))
  for (f in fl) {
    srqa <- ShortRead::qa(f, n = n)
    df <- srqa[["perCycle"]]$quality
    means <- rowsum(df$Score * df$Count, df$Cycle)/rowsum(df$Count,
                                                          df$Cycle)
    lowqual <- which(means < qscore)
    lowqual_in_first_20 <- length(which(lowqual <= 20))
    lowqual_after_20 <- lowqual[lowqual > 20][1]

    trunc_lengths <- rbind(trunc_lengths, data.frame(file = f, early_lowqual = lowqual_in_first_20,
                                                     trunc_length = ifelse(is.na(lowqual_after_20), length(means), lowqual_after_20)))

    if (quiet == FALSE) {
      if(lowqual_in_first_20 > 0){
        cat(paste0(basename(f),': ', lowqual_in_first_20, ' bases within first 20 bases are below your quality score. Consider trimming left.\n'))
      }
      if (is.na(lowqual_after_20)){
        cat(paste0(basename(f), ': After first 20 bases, no bases with a mean under your quality score found. Truncate at end of read, base: ', length(means),'\n'))

      } else if (!is.na(lowqual_after_20)){
        cat(paste0(basename(f),': After first 20 bases, first mean below your quality score is base: ',lowqual_after_20,'\n'))
        #return(lowqual_after_20)

      } else "Something's wrong. Inspect function/reads."
    } # end printing to console
  } # end loop
  return(trunc_lengths$trunc_length)
} # end function


checkPrimers.wide <- function(fwd, rev, FWD.orients, REV.orients){
    df <-  rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fwd),
          FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = rev),
          REV.ForwardReads = sapply(REV.orients, primerHits, fn = fwd),
          REV.ReverseReads = sapply(REV.orients, primerHits, fn = rev))
    df.melted <- reshape2::melt(df)
    df.melted$cat <- paste(df.melted$Var1, df.melted$Var2, sep="-")
    df.qa <- as.data.frame(t(df.melted[,c(3)]))
    names(df.qa) <- t(df.melted[,c(4)])
    rownames(df.qa) <- gsub("_R1", "", basename(fwd))
    return(df.qa)
}


## "tree" command by JennyBC
myfile <- RCurl::getURL("https://gist.githubusercontent.com/jennybc/2bf1dbe6eb1f261dfe60/raw/c53fba8a861f82f90d895458e941a1056c8118f5/twee.R", ssl.verifypeer=F)
eval(parse(text = myfile))



# Create data directory
create_data_directory <- function(path = ".", amplicon = c("ITS", "16S"), ...) {
  # twee vignette here: https://gist.github.com/jennybc/2bf1dbe6eb1f261dfe60
  if (length(amplicon) == 1) {
    cmd <- paste0("mkdir -p data/{filt_seqs/",amplicon,",raw_seqs/",amplicon,",trimmed_seqs/",amplicon,",seq_tables/",amplicon,",output_files/",amplicon,"}")
  } else {
    cmd <- "mkdir -p data/{filt_seqs/{ITS,16S},raw_seqs/{ITS,16S},trimmed_seqs/{ITS,16S},seq_tables/{ITS,16S},output_files/{ITS,16S}}"
  }
  system(cmd)
  twee("data/", level = 2)
}



rbind.named.dfs <- function(df.list){
  # solution from https://stackoverflow.com/questions/15162197/combine-rbind-data-frames-and-create-column-with-name-of-original-data-frames
  dfs <- df.list[sapply(df.list, function(x) !is.null(dim(x)))]
  all.out <- cbind.data.frame(do.call(rbind,dfs),
                              seqRun = rep(names(dfs), vapply(dfs, nrow, numeric(1))))
  return(all.out)
}


# create sample information data.frame from NEON sample names
parseNEONsampleIDs <- function(sampleID){
df <- data.frame(siteID = substr(sampleID, 1, 4), sampleID = sampleID, stringsAsFactors = F) %>%
  mutate(sample = sapply(strsplit(sampleID, "-GEN|-gen"),  "[[" , 1)) %>%
  mutate(sampleID = sapply(strsplit(sampleID, "-gen.fastq"),  "[[" , 1)) %>%
  mutate(dates = sapply(strsplit(sample, "-"), function(x) x[grep("[2]\\d\\d\\d\\d\\d\\d\\d", x)])) %>%
  mutate(dates = ifelse(dates == "21040514", "20140514", dates)) %>%
  mutate(asDate = as.Date(as.character(dates), "%Y%m%d")) %>%
  mutate(dateID = substr(as.character(dates), 1, 6)) %>%
  mutate(plotID = substr(sample, 1, 8)) %>%
  mutate(site_date = substr(sample, 1, 8)) %>%
  mutate(horizon = ifelse(grepl("-M-", sample), "M", "O")) %>%
  mutate(without_horizon = gsub("-[M|O]-", "-", sample)) %>%
  mutate(plot_date = paste0(plotID, "-", dateID)) %>%
  as.data.frame()
return(df)
}
