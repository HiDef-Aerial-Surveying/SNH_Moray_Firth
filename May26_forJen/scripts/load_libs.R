####################################################
### Script libraries function                    ###
### a function for loading libraries from a list ###
####################################################

Load_Libs <- function(LIBS){
  for(l in LIBS){
    if(require(l,character.only=TRUE,quietly=TRUE)){
      print(paste(l, "Loaded"))
    }else{
      if(l == "INLA"){
        install.packages("INLA", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
      }
      install.packages(l)
      require(l,character.only=TRUE,quietly=TRUE)
    }
  }
}

