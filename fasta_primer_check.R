####Pulls gene from fasta files using Gff3 information
### requires aligned whole genomes rooted to AF157706.1 
### requires AF157706.1 gff3 file
options(echo=TRUE)
args <- commandArgs(trailingOnly = TRUE)
#args[1] is gene name as Uxx
#args[2] forward or reverse as f or r 
#args[3] is primer sequence

library (ape)
library(seqinr)
library(Biostrings)
library(rtracklayer)

fastafile<-readDNAStringSet(file='hhv6balignedgenomes.fasta')
gfffile<-readGFF('sequence.gff3')
geneinfo<-grep(args[1],gfffile$Name)

primer<-args[3]
#deletes whitespace
primer<-gsub(" ", "", primer, fixed = TRUE)
primer<-DNAString(primer)
#if reverse, primer converts to reverse complement
if (args[2]=='r'){ 
  primer<- reverseComplement(primer)
}

#checks if gene exists in gff3 
if(isEmpty(geneinfo)==TRUE){ 
  print('gene not found')
  stop()
}

#trims fasta headers
names(fastafile)<-gsub('([A-z]+) .*', '', names(fastafile))
names(fastafile)<-trimws(names(fastafile),which='right')

#print(names(fastafile))

#for csv files
totalmatch<-c()
nomatch<-c()

#checks for primer inside of fasta files
cwd <- getwd()        
newdir <- paste0(args[3])
dir.create(newdir)      
setwd(newdir) 
for (i in 1:length(fastafile)){
  if(vcountPattern(primer,fastafile[i])==1){
    totalmatch[[i]]<-names(fastafile[i])
  }
  else{
    if((vcountPattern(primer,fastafile[i])>1)){
      print('multiple sites found for sequences:')
      print(names(fastafile[i]))
    }
  }
  if(vcountPattern(primer,fastafile[i])==0){
    for (j in 1:length(geneinfo))
      nomatch[[i]]<-names(fastafile[i])
    z<-subseq(fastafile[i], start=gfffile$start[geneinfo[1]], end=gfffile$end[geneinfo[1]], width=NA)
    q<-names(fastafile[i])
    q<-sub('\\..*', '.fasta', q)
    writeXStringSet(z, q, format='fasta')
  }
}

#omits NAs 
totalmatch<-totalmatch[!is.na(totalmatch)]
nomatch<-nomatch[!is.na(nomatch)]
#writes CSVs for match/nomatch
write.csv(totalmatch, file=paste0((args[1]),'completematch.csv'))
write.csv(nomatch, file=paste0((args[1]),'zeromatch.csv'))

system('cat *.fasta > combined.fasta')

setwd(cwd)






