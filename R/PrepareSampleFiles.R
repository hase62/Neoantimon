dir.create("lib")
setwd("lib")

urls<-c(,
        , 
        "https://github.com/hase62/Neoantimon/raw/master/lib/samtools-0.1.19.tar.bz2")

#samtools
download.file("https://github.com/hase62/Neoantimon/raw/master/lib/samtools-0.1.19.tar.bz2", 
              destfile = "samtools-0.1.19.tar.bz2")
system("tar jxf samtools-0.1.19.tar.bz2")
setwd("samtools-0.1.19")
system("make")
setwd("./../")

#ccfp
download.file("https://github.com/hase62/Neoantimon/raw/master/lib/ccfp.jar", 
              destfile = "ccfp.jar")
system("tar jxf samtools-0.1.19.tar.bz2")
setwd("samtools-0.1.19")
system("make")
setwd("./../")


download.file(, destfile = "temp")
download.file(, destfile = "temp")

wget 