#'Get Sample Files for Neoantimon
#'
#'@param url Url for getting the corresponding refFlat.txt.gz. 
#'
#'@param export_dir Export directory for using Sample Files (Default="lib").
#'
#'@return void
#'
#'@export
InstallRefFlat<-function(url, export_dir = "lib"){
  twd<-getwd()
  if(!dir.exists(export_dir)) dir.create(export_dir)
  setwd(export_dir)
  
  #refFlatfile
  file<-rev(strsplit(url, "/")[[1]])[1]
  download.file(url, destfile = file)
  system(paste("gunzip", file))
  
  setwd(twd)
}