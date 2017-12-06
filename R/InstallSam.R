#'Install Samtools
#'
#'@param export_dir Export directory for using Sample Files (Default="lib").
#'
#'@return void
#'
#'@export
InstallSam<-function(export_dir = "lib"){
  twd<-getwd()
  if(!dir.exists(export_dir)) dir.create(export_dir)
  setwd(export_dir)
  
  #samtools
  download.file("https://github.com/hase62/Neoantimon/raw/master/lib/samtools-1.6.tar.bz2", 
                destfile = "samtools-1.6.tar.bz2")
  system("tar jxf samtools-1.6.tar.bz2")
  setwd("samtools-1.6")
  system("./configure")
  system("make")
  system("make install")
  setwd(twd)
  file.remove("samtools-1.6.tar.bz2")
  
  #samtools
  download.file("https://github.com/hase62/Neoantimon/raw/master/lib/samtools-0.1.19.tar.bz2", 
                destfile = "samtools-0.1.19.tar.bz2")
  system("tar jxf samtools-0.1.19.tar.bz2")
  setwd("samtools-0.1.19")
  system("make")
  setwd(twd)
  file.remove("samtools-0.1.19.tar.bz2")
}