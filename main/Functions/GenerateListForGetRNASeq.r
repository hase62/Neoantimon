data<-t(sapply(scan(commandArgs(TRUE)[1], "character", sep="\n"), 
		function(x) strsplit(x, "\t")[[1]]))
data<-cbind(data[,3],as.numeric(data[,12]) - 2,as.numeric(data[,12]) + 2)
rownames(data)<-NULL
colnames(data)<-NULL

write.table(data, commandArgs(TRUE)[2], 
      row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")