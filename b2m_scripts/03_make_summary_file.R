# get path
path<-getwd()

## read in the files
fasta <- read.table(file = paste(path, "/working_data/06_ASV_seqs.fasta", sep=""))
counts <- read.table(file = paste(path, "/working_data/06_ASV_counts.tsv", sep=""))
taxa <- read.table(file = paste(path, "/blast_out/megan_sum_out.tsv", sep=""))

## 1). Turn the fasta into a table
names_TRUE <- grep(">", fasta$V1, value = FALSE)
seqs_TRUE <- grep(">", fasta$V1, value = FALSE, invert = TRUE)
names_vec <- fasta[names_TRUE, ]
seqs_vec<- fasta[seqs_TRUE, ]
fastatab<-cbind(gsub('>','', names_vec),gsub('>ASV_','', names_vec),seqs_vec)

## 2). Merge the fasta table by the counts table
mi_tab_1 <- merge(fastatab,counts, by.x="seqs_vec", by.y = "row.names", all=TRUE)

## 3). Label NCBI as unknown in taxa results and then merge with the previous
taxa$V2<-gsub('NCBI','unknown', taxa$V2)
mi_tab_2 <- merge(taxa,mi_tab_1, by ="V1", all=TRUE)

## 4). Rename some columns and order by ASV number
colnames(mi_tab_2)[1] <- "asv" 
colnames(mi_tab_2)[2] <- "taxa_rank"
colnames(mi_tab_2)[3] <- "taxa"
colnames(mi_tab_2)[5] <- "num_rank"
mi_tab_2$num_rank<-as.numeric(mi_tab_2$num_rank)
mi_tab_2  <- mi_tab_2 [order(mi_tab_2$num_rank),]

## 5). Save the file
write.table(mi_tab_2 [,-c(5)], "blast_out/ASV_taxa-summary_counts.tsv", sep = "\t", quote=F, col.names=T, row.names = F)
