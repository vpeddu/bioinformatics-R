#takes sam file, makes sorted bam, runs through samtools depth, and saves smoothed depth plot
#samname is sam file name
args <- commandArgs(TRUE)
samname<-args[1]
print(samname)
system(paste0('samtools view -Sb ', samname, ' > ', samname, '.bam ',sep=""))
print ('sam->bam done')
system(paste0('samtools sort ', samname, '.bam ', '-o ',samname, 'sorted.bam',sep=""))
print ('bam sorting done')
system(paste0('samtools depth ', samname,'sorted.bam',  '> ',samname,'.txt' ,sep=""))
print ('depth done')
#print(paste0('samtools view -Sb ', samname, ' > ', samname, '.bam ',sep=""))
coveragetable<-read.table(paste0(samname,'.txt'),header=FALSE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
coveragetable$V1<-NULL
colnames(coveragetable)[1]<-'position'
colnames(coveragetable)[2]<-'depth'

lowesssmoothed<-lowess(coveragetable$position, coveragetable$depth, f=.1, iter=3, delta=.01*diff(range(coveragetable$position)))
lowesssmootheddf<-as.data.frame(lowesssmoothed)
colnames(lowesssmootheddf)[1]<-'position'
colnames(lowesssmootheddf)[2]<-'depth'
png(filename=paste0(samname,'_coverage_plot.png'))
plotted<-plot(lowesssmootheddf$position,lowesssmootheddf$depth, xlab='position', ylab='depth', main=samname, ylim = c(0, 55))
dev.off()


