#'Install Samtools
#'
#'@param url Url for getting samtools (Default = "https://github.com/hase62/Neoantimon/raw/master/lib/samtools-0.1.19.tar.bz2").
#'
#'@param export_dir Export directory (Default = "lib").
#'
#'@return void
#'
#'@export
InstallSamtools<-function(url = "https://github.com/hase62/Neoantimon/raw/master/lib/samtools-0.1.19.tar.bz2",
                          export_dir = "lib"){
  file<-rev(strsplit(url, "/")[[1]])[1]
  if(dir.exists(paste(export_dir, strsplit(file, ".tar")[[1]][1], sep="/"))){
    print(paste(export_dir, "/", strsplit(file, ".tar")[[1]][1], " exists.", sep=""))
    return()
  }
  twd<-getwd()
  if(!dir.exists(export_dir)) dir.create(export_dir)
  setwd(export_dir)

  download.file(url, destfile = file)
  system(paste("tar jxf", file))
  file.remove(file)
  setwd(strsplit(file, ".tar")[[1]][1])
  system("make")
  setwd(twd)
}
