#Takes sam file and matches to ensembl biomart to find fpkm per each gene. Writes info to csv file
library(Biostrings)
library(Rsamtools)
library(GenomicAlignments)
library(biomaRt)
library(edgeR)

args = commandArgs(trailingOnly=TRUE)
filename<-args[1]
datasetname<-args[2]
#imports bam file
my.reads <- readGAlignments(file=filename,use.names = TRUE)


#write.csv(as.list(alignedreads), 'namesadsfa.txt')
#takes aligned read IDs (ensemble transcript IDs) and puts into data frame
alignedreads<-runValue(seqnames(my.reads))
names<-as.data.frame(runValue(seqnames(my.reads)))
#converts all values from binary to char
nameslist<- data.frame(lapply(names, as.character), stringsAsFactors=FALSE)
colnames(nameslist)[1]<-'Reads'
#trims off last two characters (don't know what those are for but they aren't in the ensemble transcript database)
#nameslist$Reads<-substr(nameslist$Reads,1,nchar(nameslist$Reads)-2)
nameslist[,2]<-0
colnames(nameslist)[2]<-'Genes'
uniquenames<-unique(nameslist)

#biomart pulls gene ids, transcript length, and descriptions corresponding to transcript ids 
mart <- useMart(biomart = "ensembl", dataset = datasetname)
results<-getBM(attributes = c("ensembl_transcript_id_version", "external_gene_name","description","transcript_length"),
               filters = "ensembl_transcript_id_version", values = uniquenames$Reads,
               mart = mart)
results[,5]<-0
colnames(results)[5]<-'read_count'
#matches read counts to the transcript IDs 
#necessary because biomart creates dataframe for unique values and throws out duplicates
for (i in 1:length(results$ensembl_transcript_id_version)){
  results$read_count[i]=sum(results$ensembl_transcript_id[i]==nameslist$Reads)
  print(i)
}
##24353672 reads
##5168876 aligned
results[,6]<-0
#fpkm calculation from https://www.biostars.org/p/273537/
colnames(results)[6]<-'FPKM'
# for (i in 1:length(results$ensembl_transcript_id)) { 
#   results$FPKM[i]<-(results$read_count*(1000*1000000))/(results$transcript_length[i])
#   print(i)
#   }
for (i in 1:length(results$ensembl_transcript_id)) { 
  results$FPKM[i]<-(results$read_count/((results$transcript_length[i])/1000))/(5168876/1000000)
  print(i)
}

csv_filename<-paste0(filename,".csv")
write.csv(results, csv_filename)








