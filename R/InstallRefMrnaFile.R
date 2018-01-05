#'Get refMrna file
#'
#'@param url Url for getting the corresponding refMrna.fa.gz
#'
#' (Default = "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/refMrna.fa.gz").
#'
#'@param export_dir Export directory (Default = "lib").
#'
#'@return void
#'
#'@export
InstallRefMrnaFile<-function(url = "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/refMrna.fa.gz",
                             export_dir = "lib"){
  file<-rev(strsplit(url, "/")[[1]])[1]
  if(file.exists(paste(export_dir, strsplit(file, ".gz")[[1]][1], sep="/"))){
    print(paste(export_dir, "/", strsplit(file, ".gz")[[1]][1], " exists.", sep=""))
    return()
  }
  twd<-getwd()
  if(!dir.exists(export_dir)) dir.create(export_dir)
  setwd(export_dir)

  download.file(url, destfile = file)
  system(paste("gunzip", file))



  setwd(twd)
}
