tryCatch2<-function(g){
    tryCatch(g,
             error=function(e){
              message(e)
              cat("sasasa")
                   }
             )
}