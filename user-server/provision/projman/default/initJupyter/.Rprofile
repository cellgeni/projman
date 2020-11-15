options(repos=structure(c(CRAN="https://cran.ma.imperial.ac.uk/")))
options(menu.graphics=FALSE)
options(stringsAsFactors=FALSE)
options(nwarnings=1000)
#Have warnings printed straight away
options(warn=1)
source("~/.initR/init.R")
#Should we load any of the module defined in init?
import(core,all=TRUE)
if(interactive()){
  .Last = function() try(utils::savehistory("~/.Rhistory"))
  .First = function() try(utils::loadhistory("~/.Rhistory"))
  library(colorout)
  #Load logging function
  #library(RLogBook)
  #initLogBook()
}
