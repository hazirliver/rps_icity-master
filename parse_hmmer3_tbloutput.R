library(data.table)
library(dplyr)
library(rhmmer)
library(tidyr)

args <- commandArgs(TRUE)
suffix = as.character(args[1])
fname <- paste0("hmm_out", suffix,".txt")

hmm <- read_tblout(fname)%>% 
  select_at(c(1, 4, 5)) %>% 
  separate(domain_name, into = c("gi0", "gi"), sep = "\\|", extra = "drop") %>% 
  select(-gi0)

sorted_hmm <- hmm %>%  
  group_by(gi) %>% 
  dplyr::filter(sequence_evalue == min(sequence_evalue))

gi_pfam_named <- sorted_hmm[,c(1,2)]
colnames(gi_pfam_named) <- c("gi", "pfam_cog_id")


fwrite(sorted_hmm, paste0("sorted_hmm", suffix, ".tsv"), quote = F, sep = "\t")
fwrite(gi_pfam_named, paste0("gi_pfam_named", suffix,".tsv"), quote = F, sep = "\t")
