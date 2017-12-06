#'Get refMrna file
#'
#'@param url Url for getting the corresponding refMrna.fa.gz. 
#'
#'@param file The file path to refMrna.fa.gz or refMrna.fa. 
#'
#'@param export_dir Export directory when using url (Default="lib").
#'
#'@return void
#'
#'@export
InstallRefMrnaFile<-function(url = NULL, 
                             file = NULL, 
                             export_dir = "lib"){
  if(is.null(url) & is.null(file)) print("Please specify url or file.")

  twd<-getw()
  if(!is.null(url)){
    if(!dir.exists(export_dir)) dir.create(export_dir)
    setwd(export_dir)
    if(is.null(file)) file<-rev(strsplit(url, "/")[[1]])[1]
    download.file(url = url, destfile = file)
  }
  if(!is.na(grep(".gz", file))) {
    system(paste("gunzip", file))
    file<-strsplit(file, ".gz")[[1]][1]
  }  
  data<-scan(file, "character", "\n")
  ID<-grep(">", data)
  Sp_var<-ID + 1
  Seq_start<-ID + 2
  Seq_end<-c(ID[-1] - 1, length(data))
  merge<-cbind(gsub(">", "", data[ID]), 
        data[Sp_var], 
        sapply(1:length(Seq_start), function(x) paste(data[Seq_start[x]:Seq_end[x]], collapse = ""))
  )
  write.table(merge, paste(file, "merge.fa", sep = "."), row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
}
