#'Get refMrna file
#'
#'@param url Url for getting the corresponding refMrna.fa.gz (Default = "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/refMrna.fa.gz"). 
#'
#'@param export_dir Export directory (Default = "lib").
#'
#'@return void
#'
#'@export
InstallRefMrnaFile<-function(url = "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/refMrna.fa.gz", 
                             export_dir = "lib"){
  twd<-getwd()
  if(!dir.exists(export_dir)) dir.create(export_dir)
  setwd(export_dir)

  file<-rev(strsplit(url, "/")[[1]])[1]
  download.file(url, destfile = file)
  system(paste("gunzip", file))
  
  setwd(twd)
}
