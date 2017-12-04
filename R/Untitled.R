install.packages("rJava")
library(rJava)
dyn.load("/Library/Java/JavaVirtualMachines/jdk1.8.0_144.jdk/Contents/Home/jre/lib/server/libjvm.dylib")

SET_CCFP <- list(
  CLASS_JAR_PATH = stringr::str_c(getwd(), "../lib/ccfp.jar")
)

ccfp <- function (str, set_mode = NULL, jar_path) {
  if (is.null(str) || str == "") {
    return (NULL)
  }
  
  CLASS_JAR_PATH = stringr::str_c(getwd(), "/../lib/ccfp.jar")
  rJava::.jinit(classpath = "/Users/takanorihasegawa/Git/Neoantimon/lib/ccfp.jar", force.init = TRUE)
  #rJava::.jaddClassPath(path = )
  rJava::J("java.lang.Math")
  rJava::J("MUT.clone.CCFP")
  
  mode <- rJava::J(class = "MUT.clone.CCFP")
  builder <- rJava::J(class = "java.clone.CCFP", method = "builder")
  if (!is.null(set_mode)) {
    if (set_mode == "SEARCH") {
      builder <- builder$mode(mode$SEARCH)
    } else if (set_mode == "EXTENDED") {
      builder <- builder$mode(mode$EXTENDED)      
    } else {
      set_mode <- "NORMAL"
    }
  } else {
    set_mode <- "NORMAL"
  }
  
  tokenizer <- builder$build()
  tokened <- tokenizer$tokenize(str)

  return (
    dplyr::bind_rows(
      lapply(X = as.list(tokened), FUN = function (token) {
        return(
          dplyr::data_frame(
            surface = token$getSurfaceForm(), feature = token$getAllFeatures(),
            is_know = token$isKnown(), is_unk = token$isUnknown(), is_user = token$isUser(),
            mode = set_mode
          )
        )
      })
    )
  )
}