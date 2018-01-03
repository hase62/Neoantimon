#'Get Sample Files for Neoantimon
#'
#'@param export_dir Export directory for using Sample Files (Default="lib").
#'
#'@return void
#'
#'@export
InstallSampleFiles<-function(export_dir = "lib"){
  twd<-getw()
  if(!dir.exists(export_dir)) dir.create(export_dir)
  setwd(export_dir)

  #sample files
  download.file("https://github.com/hase62/Neoantimon/raw/master/lib/data.zip", 
                destfile = "data.zip")
  system("unzip data.zip")

  setwd(twd)
}