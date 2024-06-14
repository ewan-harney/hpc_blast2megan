# get path
path<-getwd()

## read in the files
fasta <- read.table(file = paste(path, "/working_data/06_ASV_seqs.fasta", sep=""))
counts <- read.table(file = paste(path, "/working_data/06_ASV_counts.tsv", sep=""))
taxsum <- read.table(file = paste(path, "/blast_out/megan_summary_out.tsv", sep=""))
taxpath <- read.table(file = paste(path, "/blast_out/megan_taxonpath_out.tsv", sep=""))

## 1). Turn the fasta into a table
names_TRUE <- grep(">", fasta$V1, value = FALSE)
seqs_TRUE <- grep(">", fasta$V1, value = FALSE, invert = TRUE)
names_vec <- fasta[names_TRUE, ]
seqs_vec <- fasta[seqs_TRUE, ]
fastatab <-cbind(gsub('>','', names_vec),gsub('>ASV_','', names_vec),seqs_vec)

## 2). Merge the fasta table by the counts table
mi_tab_1 <- merge(fastatab,counts, by.x="seqs_vec", by.y = "row.names", all=TRUE)

## 3). Merge taxon path with the sequence + counts
mi_tab_path <- merge(taxpath,mi_tab_1, by ="V1", all=TRUE)

## 4). Rename columns in taxon path and order by ASV number
colnames(mi_tab_path)[1] <- "asv" 
colnames(mi_tab_path)[2] <- "domain"
colnames(mi_tab_path)[3] <- "kingdom"
colnames(mi_tab_path)[4] <- "phylum"
colnames(mi_tab_path)[5] <- "class"
colnames(mi_tab_path)[6] <- "order"
colnames(mi_tab_path)[7] <- "family"
colnames(mi_tab_path)[8] <- "genus"
colnames(mi_tab_path)[9] <- "species"
colnames(mi_tab_path)[10] <- "subspp"
colnames(mi_tab_path)[12] <- "asv_number"

## 5). Rename columns in taxa summary file
colnames(taxsum)[1] <- "asv" 
colnames(taxsum)[2] <- "lca_taxa_rank"
colnames(taxsum)[3] <- "lca_taxa_name"

## 6). Merge summary with path + seq + counts and order by asv number
mi_tab_final <- merge(taxsum,mi_tab_path, by ="asv", all=TRUE)
mi_tab_final$asv_number<-as.numeric(mi_tab_final$asv_number)
mi_tab_final  <- mi_tab_final [order(mi_tab_final$asv_number),]

## 5). Save the file, removing asv_number
write.table(mi_tab_final [,-c(14)], "working_data/ASV_taxa_seq_counts.tsv", sep = "\t", quote=F, col.names=T, row.names = F)
