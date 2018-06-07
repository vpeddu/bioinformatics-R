#Takes bam file, start position, end position, and creats fasta file
library (Rsamtools)
library(GenomicAlignments)

options(echo=TRUE)
args <- commandArgs(trailingOnly = TRUE)
#args[1]=bamfilename
#args[2]=start position 
#args[3]= end position 


#sets parameters for scanbam
params<-ScanBamParam(flag=scanBamFlag(isUnmappedQuery=FALSE), what=c('rname','qname','strand','qwidth','cigar','seq'))

#defines index for bamfile 
baifname<-indexBam(args[1])

#scans bam file 
bamfile<- scanBam(args[1],param=params, index=baifname)
startpos<-as.numeric(args[2])
endpos<-as.numeric(args[3])
stackparam<-RangesList(IRanges(start=startpos, end=endpos))

names(stackparam)<-unique(bamfile[[1]]$rname)
stackedbam<-stackStringsFromBam(args[1],index=baifname, param=ScanBamParam(which=stackparam))
names(stackedbam)<-bamfile[[1]]$qname
fastaname<-paste0(args[1],'.fasta')
writeXStringSet(stackedbam, fastaname, format='fasta')

quit(save="no")



