#Remove all Ns at the beginning and end of the seq to prep files for prokka annotation 

rm(list=ls()); 
library(Biostrings)

options(echo=TRUE)
args <- commandArgs(trailingOnly = TRUE)

sequences<-readDNAStringSet(args[1],format='fasta')

sequencestrimmed<-DNAStringSet(gsub("N*N$",'',gsub("^N*",'',as.character(sequences))));
names(sequencestrimmed)<-substring(names(sequences),1,20); 
#prokka needs contig name to be <=20 chars long
fastafilename<-as.character(names(sequencestrimmed))
fastafilename<-paste(fastafilename, ".fasta")

writeXStringSet(sequencestrimmed,file=fastafilename,format='fasta');
file.remove(args[1])
quit(save="no")
