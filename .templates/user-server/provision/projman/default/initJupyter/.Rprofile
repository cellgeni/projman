###################
# Standard options
options(repos=structure(c(CRAN="https://cran.ma.imperial.ac.uk/")))
options(menu.graphics=FALSE)
options(stringsAsFactors=FALSE)
options(nwarnings=1000)
#Have warnings printed straight away
options(warn=1)
#Try and persist history
if(interactive()){
  .Last = function() try(utils::savehistory("~/.Rhistory"))
  .First = function() try(utils::loadhistory("~/.Rhistory"))
  #Uncomment this for glorious technicolour
  #library(colorout)
}
############################
# Define the import function
#Paths to search for importing "modules"
assign('.modPaths',c('.','~/trueHome/Projects/Common/Code/','~/trueHome/Projects/Common/Code/INITR','~/.initR/'),envir = .GlobalEnv)
#' Imports a python style module.  These are just files with .R/r extensions.
#' 
#' If src is an R file, it is used directly.  Otherwise the directories in searchDirs are searched for an R file named as \code{src}.  The file is loaded into an environment that is stored in the callers environment.
#'
#' @param src File to load or name of module in search path.
#' @param as What to name the loaded module.
#' @param searchDirs Directories to search for module.
#' @param all Load all entries into parent namespace.  Use with caution.
#' @param reattach If module is already loaded with all=TRUE, re-attach it.
import = function(src,as=NULL,searchDirs=.modPaths,all=FALSE,reattach=TRUE){
  #Force it to be a character
  src = as.character(substitute(src))
  #Get the calling environment
  #penv = parent.env(environment())
  #Make a new environment for the src file
  module = new.env()
  #Get the thing
  if(!file.exists(src)){
    #Try and find the thing in each of the search paths
    found=FALSE
    for(searchDir in searchDirs){
      if(file.exists(file.path(searchDir,sprintf('%s.r',src)))){
        src = file.path(searchDir,sprintf('%s.r',src))
        found=TRUE
        break
      }else if(file.exists(file.path(searchDir,sprintf('%s.R',src)))){
        src = file.path(searchDir,sprintf('%s.R',src))
        found=TRUE
        break
      }
    }
    if(!found)
      stop(sprintf("Could not find module %s on search path.",src))
  }
  sys.source(src,envir=module,keep.source=interactive())
  #Work out what to call it
  if(is.null(as))
    as=gsub('\\.[Rr]$','',basename(src))
  #Now do the assignment
  if(all){
    nom = sprintf("module:%s",as)
    if(nom %in% search()){
      if(reattach){
        detach(pos=match(nom,search()))
      }else{
        stop(sprintf("A module named %s is already attached.",as))
      }
    }
    attach(module,pos=2,name=nom)
  }else{
    penv = parent.frame()
    assign(as,module,envir=penv)
  }
  invisible(NULL)
}

################
# Core functions
#' Close all display devices
dev.off.all = function() {
  done = tryCatch(dev.off(),error = function(e) FALSE)
  for(i in seq(100)){
    done = tryCatch(dev.off(),error = function(e) FALSE)
    if(done){
      break
    }
  }
}

#' Memory usage of objects
#' 
#' @param what The names of the objects to check.  Usually ls()
usedMem = function(what) {
  x = lapply(what, function(e) object.size(get(e)))
  o = order(sapply(x,as.integer),decreasing=TRUE)
  x = sapply(x,format,unit='auto')
  names(x)=what
  x[o]
}

#' Replace variables in string
#' @param string The string to replace variables in.
#' @param ... What to replace with what.
varString = function(string,...) {
  out = string
  vars = list(...)
  for(i in seq_along(vars)){
    out = gsub(paste0('%',names(vars)[i],'%'),vars[[i]],out)
  }
  return(out)
}

#' An enhancement to the standard table function 
#'
#' @param levels These are entries that are always in the table output 
#' @param dropExtra Should we drop any value not specified?
#' @param defValue Default value
ltable = function(...,levels=NULL,dropExtra=TRUE,defValue=0) {
  if(is.null(levels)){
    table(...)
  }else{
    x = table(...)
    tmp = rep(defValue,length(levels))
    names(tmp)=levels
    if(dropExtra){
      o = which(names(x)%in%levels)
      tmp[names(x[o])] = x[o]
    }else{
      tmp[names(x)] = x
    }
    return(tmp)
  }
}

#' Determines how many lines of comments in a file
#'
#' Reads the file a few lines at a time until a non-comment line is reached.
#'
#' @param fnom The file name of the file.
#' @param comment.char The character that indicates a commented line
#' @return The number of commented lines at the start of the file.
numComments = function(fnom,comment.char='#'){
  guess=10
  repeat{
    dat  = substr(readLines(fnom,guess),1,1)!=comment.char
    x = which(dat)
    #Did we find a non comment
    if(length(x)>0)
      break
    #Did we read the whole file
    if(length(dat)<guess)
      return(length(dat))
    guess=guess+10
  }
  return(x[1]-1)
}
