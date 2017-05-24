tryCatch2<-function(g){
    tryCatch(g,
             error=function(e){
              message(e)
              cat("sasasa")
                   }
             )
}

apply2<-function(x, y, f){
  if(is.null(nrow(x))){
    return(f(x))
  } else {
    apply(x, y, function(a) f(a))
  }
}
