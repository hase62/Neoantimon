data<-t(sapply(scan(commandArgs(TRUE)[1], "character", sep="\n"), 
	function(x) strsplit(x, "\t")[[1]]))
colnames(data)<-NULL
rownames(data)<-NULL

if(is.na(commandArgs(TRUE)[2])){
   ratio_matrix<-matrix(nrow=nrow(data), ncol=6, NA)
   write.table(cbind(data, ratio_matrix), commandArgs(TRUE)[1],
      row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
   q("no")
}

#name, mean, min, max
data_2<-t(sapply(scan(commandArgs(TRUE)[2], "character", sep="\n", skip=1), 
	function(x) strsplit(x, "\t")[[1]]))
colnames(data_2)<-NULL
rownames(data_2)<-NULL
if(length(data_2)==0){
   ratio_matrix<-matrix(nrow=nrow(data), ncol=6, NA)
   write.table(cbind(data, ratio_matrix), commandArgs(TRUE)[1],
      row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
   q("no")
}

tmp<-matrix(nrow=nrow(data), ncol=ncol(data_2), -1)
tmp[match(as.numeric(sapply(data_2[,1], function(x) strsplit(x, ";")[[1]][2])), 1:nrow(data)),]<-data_2

if(nrow(data)==1){
   write.table(t(c(data, tmp[,c(1,2,3,4)])), paste(commandArgs(TRUE)[1],sep=""), 
      row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
}else{
   write.table(cbind(data, tmp[,c(1,2,3,4)]), paste(commandArgs(TRUE)[1],sep=""), 
      row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
}

