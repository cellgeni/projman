#Core functions and variables that we almost always want to load
#Some libraries are so universal, we may as well just load them at the start

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

#' Re-set the DISPLAY environmental variable so plotting will work again
#' If n is not given, will loop through possible options until plotting stops throwing an error.
#' @param n The display number to set.
fixDisp = function(n) {
  if(missing(n)){
    #Close everything
    dev.off.all()
    #Set display until it works
    working=FALSE
    for(k in 1:20){
      fixDisp(k)
      working = suppressWarnings(is.null(tryCatch(plot(1),error = function(e) FALSE)))
      dev.off.all()
      if(working){
        break
      }
    }
  }else{
    curr = Sys.getenv('DISPLAY')
    Sys.setenv(DISPLAY= gsub(':([0-9]+)\\.',sprintf(':%d.',n),curr))
  }
}

#' Memory usage of objects
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
