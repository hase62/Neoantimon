#'Get refFlat file
#'
#'@param url Url for getting the corresponding refFlat.txt.gz (Default = "http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refFlat.txt.gz"). 
#'
#'@param export_dir Export directory (Default = "lib").
#'
#'@return void
#'
#'@export
InstallRefFlat<-function(url = "http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refFlat.txt.gz", 
                         export_dir = "lib"){
  twd<-getwd()
  if(!dir.exists(export_dir)) dir.create(export_dir)
  setwd(export_dir)
  
  file<-rev(strsplit(url, "/")[[1]])[1]
  download.file(url, destfile = file)
  system(paste("gunzip", file))
  
  setwd(twd)
}