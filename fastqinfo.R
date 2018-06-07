#given file and path makes csv of fastq stats
file<-args[1]
filepath<-args[2]
system("cd HSV_WGS
        cd sequencing_runs") 
filenames<-list.files(file)

for (files in filenames) { 
x<-paste(filepath,sep = '')
assign(paste0(files),(list.files(x)))
}
all<-c()
for(i in 1:length(filenames)) {
  all[[i]]<-get(filenames[i])
}
all<-unlist(all)
namedf<-as.data.frame(all)
names(namedf)[1]<-'Filenames'

#checks for VRC
namedf[2]<-grepl("*(\\d{4})-", namedf$Filenames)
names(namedf)[2]<-'VRC'

#checks for duplicates
namedf[3]<-duplicated(namedf$Filenames)
names(namedf)[3]<-'Duplicates'

#finds file size
filesize<-c()
reads<-c()
for(i in 1:nrow(namedf)) {
  for (j in 1:length(filenames)) {
    x<-paste(filepath,namedf$Filenames[i],sep='')
    if(file.exists(x)==TRUE) {
      filesize[i]<-file.size(x)
      #how to get output from console into vector?
      #reads[i]<-system2(system(paste("awk '{s++}END{print s/4}'",x)),stdout = '')
      #bash commands for file size
      #print(system(paste("du -h",x)))
    }
  }
}
#filesize gives size in bytes
namedf[4]<-(filesize/1000000)
names(namedf)[4]<-'Size (MB)'

#writes to csv
write.csv(namedf, file="fastq_info.csv")






