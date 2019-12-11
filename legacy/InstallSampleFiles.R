#'Get Sample Files for Neoantimon
#'
#'
#'
#'
#'
#'@param export_dir Export directory (Default = "lib").
#'
#'@return void
#'
#'@export
InstallSampleFiles<-function(export_dir = "lib"){
  url <- "https://github.com/hase62/Neoantimon/raw/master/lib/data.zip"
  file<-rev(strsplit(url, "/")[[1]])[1]
  if(dir.exists(paste(export_dir, strsplit(file, ".zip")[[1]][1], sep="/"))){
    print(paste(export_dir, "/", strsplit(file, ".zip")[[1]][1], " exists.", sep=""))
    return()
  }
  twd<-getwd()
  if(!dir.exists(export_dir)) dir.create(export_dir)
  setwd(export_dir)

  download.file(url, destfile = file)
  system(paste("unzip", file))
  file.remove(file)


  setwd(twd)
}
