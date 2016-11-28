index<-as.numeric(commandArgs(TRUE)[1])
num<-as.numeric(commandArgs(TRUE)[2])

data<-t(sapply(scan("RNAExpression.txt","character",sep="\n"), function(x) strsplit(x,"\t")[[1]]))
unq<-unique(data[,1])
ran<-c(floor(length(unq)/num*(index-1))+1, floor(length(unq)/num*index))
ran[2]<-ifelse(ran[2] > length(unq), length(unq), ran[2])
u<-unq[ran[1]:ran[2]]
res<-sapply(u, function(x) max(as.numeric(data[!is.na(match(data[,1], x)),2])))

write.table(res, "RNAExpression.max.txt", row.names=TRUE, col.names=FALSE,quote=FALSE,sep="\t",append=TRUE)
