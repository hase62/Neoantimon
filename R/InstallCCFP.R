#'Get CCFP.jar to Calculate Cancer Cell Fraction. 
#'
#'@param export_dir Export directory for using Sample Files (Default="lib").
#'
#'@return void
#'
#'@export
GetCCFP<-function(export_dir = "lib"){
  twd<-getwd()
  if(!dir.exists(export_dir)) dir.create(export_dir)
  setwd(export_dir)
     
  #ccfp
  download.file("https://github.com/hase62/Neoantimon/raw/master/lib/ccfp.jar", 
                destfile = "ccfp.jar")

  setwd(twd)
}