#!/usr/bin/R
#This is a master script that creates a bunch of "import" functions that can be used to load various convenience functions

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

##Core utilities and quality of life improvement functions
#importCore = function(src='~/Projects/Common/Code/INITR/core.R'){
#  source(src)
#}
#
##Functions for working with the Sanger infrastructure
#importSanger = function(src='~/Projects/Common/Code/INITR/sanger.R'){
#  source(src)
#}
#
##Function for working with genomic objects
#importGenome = function(src='~/Projects/Common/Code/INITR/genome.R'){
#  source(src)
#}
#
##Some customs plotting functions that are commonly useful
#importPlot = function(src='~/Projects/Common/Code/INITR/plot.R'){
#  source(src)
#}


  

##DEPRECATED: Now part of TCN function.  
###' Uses the method of Speed and ? to calculate rates at which we expect reads to occur based on GC content.  Returns predicted counts for the regions given in tgts.
###' For now this is intended to be used with targeted sequencing data and so the mapability of sequence is ignored.
###' To handle an entire genome, set tgts to the entire genome, minus unmappable locations e.g. tgts = setdiff(getGenome('hg19'),unmappableLocationsGRanges)
###' @param bfile The path to a BAM file containing the mapped reads.
###' @param tgts A GRanges object giving the locations of the genome to consider.  
###' @param genome The BSgenome object to use.
###' @param lstrip The number of bases to trim from the 5' end of region.
###' @param rstrip The number of bases to trim from the 3' end of region.
###' @param gc_smooth Smooth the predicted rates in bins of this width.  If less than zero interpreted as fraction of insert size.
##GCnorm = function(bfile,tgts,genome=Hsapiens,lstrip=3,rstrip=3,gc_smooth = .03,maxReads=10000000) {
##  require(GenomicFiles)
##  #Make a non-overlapping version of targets for 
##  rtgts = reduce(tgts)
##  #Filter out bad reads
##  #Keep only the first mate of a read, discard anything even mildly suspicious
##  flags = scanBamFlag(isPaired = TRUE, isProperPair = TRUE, isUnmappedQuery = FALSE,
##         hasUnmappedMate = FALSE, isMinusStrand = NA, isMateMinusStrand = NA,
##         isFirstMateRead = TRUE, isSecondMateRead = FALSE, 
##         isSecondaryAlignment = FALSE, isNotPassingQualityControls = FALSE,
##         isDuplicate = NA)
##  bf = open(BamFile(bfile,yieldSize=maxReads))
##  #Get the positions of any read in the target region
##  param = ScanBamParam(flags,what=c('rname','pos','mpos','qwidth','isize'))
##  yield = function(e) scanBam(e,param=param)
##  map = function(e) {
##      pos = unlist(lapply(e,`[[`,'pos'),use.names=FALSE)
##      mpos = unlist(lapply(e,`[[`,'mpos'),use.names=FALSE)
##      qwidth = unlist(lapply(e,`[[`,'qwidth'),use.names=FALSE)
##      strand = pos < mpos
##      tmp = GRanges(unlist(lapply(e,`[[`,'rname')),
##              IRanges(ifelse(strand,pos,pos+qwidth-1),width=1),
##              strand=ifelse(strand,'+','-'),
##              readLen=qwidth,
##              insert=unlist(lapply(e,`[[`,'isize'),use.names=FALSE)
##              )
##      subsetByOverlaps(tmp,rtgts)
##  }
##  #This is just a slightly edited REDUCEsampler
##  sampleSize = maxReads
##  reducer = function(x, y, ...) {
##    #Haven't gotten anything useful yet...
##    if(length(x)==0L & length(y)==0L){
##      cat('Skipping...\n')
##      return(x)
##    }
##    #print(x)
##    yld_n = length(y)
##    yld_for_n = sum(strand(y)=='+')
##    yld_rev_n = sum(strand(y)=='-')
##    #print(c(tot,length(x),length(y)))
##    #If it's the first time, store length of x
##    if(!('tot' %in% colnames(mcols(x)))) {
##      x$tot = length(x)
##      x$nFor = sum(strand(x)=='+')
##      x$nRev = sum(strand(x)=='-')
##      #If we've got options, subsample
##      if(length(x) >= sampleSize){
##        x = x[sample(length(x),sampleSize)]
##      }
##    }
##    #Establish counters in y
##    y$tot = numeric(length(y))+x$tot[1]
##    y$nFor = numeric(length(y))+x$nFor[1]
##    y$nRev = numeric(length(y))+x$nRev[1]
##    #Fill up x, if it's not already filled up
##    if(unique(x$tot) < sampleSize) {
##      #x is established 
##      #First fill up x as much as we can
##      x = append(x,y)
##      #If we've gone too far, trim off the fat as y
##      if(length(x)>sampleSize){
##        y = x[seq(sampleSize+1,length(x))]
##        x = x[seq(sampleSize)]
##      }else{
##        #Got nothing to add then...
##        y=GRanges(tot=numeric(),nFor=numeric(),nRev=numeric())
##      }
##    }
##    #Update the counters
##    x$tot = x$tot + yld_n
##    x$nFor = x$nFor + yld_for_n
##    x$nRev = x$nRev + yld_rev_n
##    #Randomly add in some new ones from y
##    keep = rbinom(1L, min(sampleSize,length(y),length(x)),min(length(x),length(y))/unique(x$tot))
##    #Print
##    cat(sprintf("On target reads (+,-) = %g (%g,%g)\n", unique(x$tot),unique(x$nFor),unique(x$nRev)))
##    i = sample(length(x),keep)
##    j = sample(length(y),keep)
##    x[i] = y[j]
##    x
##  }
##  done = function(e) length(e[[1]][[1]]) == 0L
##  #Actually load a subset of the reads
##  cat('Loading reads.\n')
##  reads = reduceByYield(bf,yield,map,reducer,done)
##  #Close the file
##  close(bf)
##  #If it hasn't gone through the reducer, we have all the reads
##  if(length(reads)<sampleSize){
##    reads = reducer(reads,GRanges())
##  }
##  #How many reads did we get?
##  Nfor = reads$nFor[1]
##  Nrev = reads$nRev[1]
##  Ntot = reads$tot[1]
##  #Now get a fragment size distribution
##  isize = median(abs(reads$insert))
##  readLen = median(reads$readLen)
##  wlen = isize-lstrip-rstrip
##  #Sanitise reads and regions to create a list
##  reads = sanitiseRanges(reads,genome,isize-rstrip)
##  rtgts = sanitiseRanges(rtgts,genome,isize-rstrip)
##  reg_list = list(readSubset=reads,targetGenome=rtgts)
##  #Get GC counts at locations we care about
##  cnts = get_GC_counts(reg_list,isize,lstrip,rstrip,genome=genome,summarise='Region')
##  #Finally calculate the GC normalisation rates
##  #Smooth the result and scale appropriately
##  if(gc_smooth<1){
##    gc_smooth = max(1,ceiling(wlen*gc_smooth))
##  }
##  gc_smooth = gc_smooth + !(gc_smooth%%2)
##  pad = rep(0,gc_smooth %/% 2)
##  #Add a rates column to the output
##  noms = c(rownames(cnts$forw),'rate')
##  cnts$forw = rbind(cnts$forw,rep(0,ncol(cnts$forw)))
##  rownames(cnts$forw) = noms
##  #Smooth the result
##  cnts$forw['rate',] = as.numeric(runsum(Rle(c(pad,cnts$forw['readSubset',],pad)),gc_smooth)) / as.numeric(runsum(Rle(c(pad,cnts$forw['targetGenome',],pad)),gc_smooth))
##  #Dodgy values should be 0
##  cnts$forw['rate',][!is.finite(cnts$forw['rate',])]=0
##  #Scale it appropriately
##  cnts$forw['rate',] = cnts$forw['rate',] * Nfor/sum(cnts$forw['rate',] * cnts$forw['targetGenome',])
##  #Same for reverse
##  noms = c(rownames(cnts$revr),'rate')
##  cnts$revr = rbind(cnts$revr,rep(0,ncol(cnts$revr)))
##  rownames(cnts$revr) = noms
##  cnts$revr['rate',] = as.numeric(runsum(Rle(c(pad,cnts$revr['readSubset',],pad)),gc_smooth)) / as.numeric(runsum(Rle(c(pad,cnts$revr['targetGenome',],pad)),gc_smooth))
##  cnts$revr['rate',][!is.finite(cnts$revr['rate',])]=0
##  cnts$revr['rate',] = cnts$revr['rate',] * Nrev/sum(cnts$revr['rate',] * cnts$revr['targetGenome',])
##  #Return everything we need to make use of this
##  return(list(forw=cnts$forw,revr=cnts$revr,isize=isize,lstrip=lstrip,rstrip=rstrip,readLen=readLen,genome=genome,bfile=bfile,gc_smooth=gc_smooth))
##}



##DEPRECATED: Now part of TCN function.  
###' Truncates or drops regions that are within buf of the start/end of chromosomes
##sanitiseRanges = function(regions,genome,buf){
##  clens = seqlengths(genome)
##  names(clens) = fix_chr(names(clens))
##  #Find ones that are too close to start
##  o = which(start(regions) < buf)
##  #Drop anything that's entirely too close
##  oo = end(regions[o]) < buf
##  to_drop = o[oo]
##  #Edit the start of the others
##  start(regions[o[!oo]]) = buf
##  #Find those too close to the end
##  o = which(end(regions) > clens[fix_chr(as.character(seqnames(regions)))] - buf+1)
##  #Drop anything that's no good
##  oo = start(regions[o]) > clens[fix_chr(as.character(seqnames(regions[o])))] - buf+1
##  to_drop = c(to_drop,o[oo])
##  #Fix the other
##  end(regions[o[!oo]]) = clens[fix_chr(as.character(seqnames(regions[o[!oo]])))] -buf +1
##  #Drop the duds and return
##  if(length(to_drop)>0){
##    regions[-to_drop]
##  }else{
##    regions
##  }
##}


##DEPRECATED: Now part of TCN function.  
###' Gets the counts of window GC content at the windows anchored on the positions given. 
###'
###' For each location referred to by the GRanges contained within reg_list, the GC content in two surrounding windows is calculated.  The first, the "forward window", start at the genomic location of interest + lstrip and is isize-lstrip-rstrip bases long.  The second, the "reverse window", ends at the genomic location of interest - lstrip and extends backwards (towards lower genomic coordinates) for isize-lstrip-rstrip bases.
###'
###' There are three different methods used to fetch this information depending on the size and length of the genomic ranges requested.  
###' When the number of genomic locations to interrogate is small and they cover only a small fraction of the genome, the sequence for each window is simply directly fetched.
###' When only a small fraction of genome is covered by a large number of independent regions, we fetch the sequence for an overlapping larger region and then extract the GC content at the locations of interest.
###' When a large fraction of the genome is required, the entire genome is fetched and GC content calculated in a sliding window.  The required windows are then extracted by sub setting this larger list.
###'
###' It is useful to specify multiple regions of interest at the same time if possible as this prevents the genome sequence having to fetched multiple times.
###'
###' After the GC has been calculated in each of the desired windows, optional aggregation is performed.  Returning either the GC count at each window, the frequency of windows with different counts for each separate GRange specified or the frequency of windows with different counts for each region in reg_list.
###'
###' @param reg_list A list (or GRangesList) of regions where we want to know the GC content in a window anchored at the contained locations.  For example, a list might contain a set of read locations in the first entry and a mappability mask of the genome in the second.
###' @param isize The insert size to use.
###' @param lstrip The number of bases to trim at the 5' end of the window.
###' @param rstrip The number of bases to trim at the 3' end of the window.
###' @param summarise One of 'Region', 'GRange', or 'BP'.  If 'BP', the GC count in each anchored window is returned.  If 'Region' or 'Grange', returns the frequency of window GC counts aggregated across each individual GRange entry (if 'GRange' specified) or each list entry in reg_list (if 'Region' specified).
###' @param genome The BSgenome object from which genome sequence is extracted.
###' @param min_cvg_frac In deciding how to fetch data, the number of entries in the covering GRanges objects is compared to the number in each reg_list.  If this fraction is larger than min_cvg_frac (i.e., we don't save many fetches by loading the covering regions rather then the targets directly), then the direct fetch method is used.
###' @param max_genome_frac If the requested region covers more than this fraction of the genome, then fetch the entire genome for simplicity and speed.
###' @param max_quick_regions Don't try and fetch things the quick way if we have more than this many distinct regions to fetch.
###' @return If summarise is 'Region', a list of matricies with length(reg_list) rows and columns giving the number of windows with specified GC count in that region.  Otherwise reg_list is returned with extra columns 'forw' and 'revr' added holding either the frequencies of GC windows (if summarise is 'GRange') or the individual GC counts in order (if summarise is 'BP').
##get_GC_counts = function(reg_list,isize,lstrip=3,rstrip=3,summarise=c('Region','GRange','BP'),genome=Hsapiens,min_cvg_frac=0.5,max_genome_frac=0.04,max_quick_regions=2000){
##  summarise = match.arg(summarise)
##  wlen = isize-rstrip-lstrip
##  if(summarise=='Region'){
##    #Create the output matrix
##    out_for = matrix(0,nrow=length(reg_list),ncol=wlen+1)
##    rownames(out_for)=names(reg_list)
##    colnames(out_for)= 0:wlen
##    out_rev = out_for
##  }else{
##    #Create place to store individual counts or summary vector if we need that
##    reg_list = lapply(reg_list,function(e){
##                        e$forw=vector('list',length(e))
##                        e$revr=vector('list',length(e))
##                        e})
##  }
##  #Get the chromosome lengths
##  clens = seqlengths(genome)
##  names(clens) = fix_chr(names(clens))
##  #Get the list of chromosomes we need to look at
##  chrs = unique(do.call(c,lapply(reg_list,function(e) unique(fix_chr(as.character(seqnames(e)))))))
##  #Order them sensibly
##  chrs = chrs[order(num_chr(chrs))]
##  #The remapping of things to UCSC genome chr format 
##  chrFormatMap = setNames(prefix_chr(chrs),chrs)
##  #Check that we haven't asked for bad regions
##  if(any(sapply(reg_list,function(e) any(start(e)<lstrip+wlen) | any(end(e) > clens[fix_chr(as.character(seqnames(e)))]-lstrip-wlen+1)))){
##    stop('Location requested without any valid window.')
##  }
##  #Make a coverage object so we know where we need data
##  cvg = reduce(GRanges(unlist(lapply(reg_list,function(e) fix_chr(as.character(seqnames(e)))),use.names=FALSE),
##                       IRanges(unlist(lapply(reg_list,start),use.names=FALSE),
##                               unlist(lapply(reg_list,end),use.names=FALSE))))
##  #Add flank that we'll need
##  start(cvg) = pmax(1,start(cvg)-lstrip-wlen+1)
##  end(cvg) = pmin(clens[fix_chr(as.character(seqnames(cvg)))],end(cvg)+lstrip+wlen-1)
##  #Reduce it again
##  cvg = reduce(cvg)
##  lcvg = length(cvg)
##  #Split by (fixed up) chromosome
##  cvg = split(cvg,as.character(seqnames(cvg)))
##  #Check if we can do any of them super quickly
##  todo = rep(TRUE,length(reg_list))
##  for(i in seq_along(reg_list)){
##    #What fraction of the genome is covered
##    genome_frac = sum(as.numeric(width(reg_list[[i]])))/sum(as.numeric(seqlengths(genome)))
##    #How long is the reduced set of regions we need to fetch compared to getting them all
##    cvg_frac = lcvg/length(reg_list[[i]])
##    #If we can do it quickly, do it quickly...
##    #Do the quick version if we're only having to fetch 20% or fewer regions if we do it the long way and the total fraction of the genome asked for is small so we don't run into memory issues
##    if(lcvg < max_quick_regions & cvg_frac > min_cvg_frac & genome_frac<max_genome_frac) {
##      tmp = get_GC_counts_for_region(reg_list[[i]],isize,lstrip,rstrip,summarise=summarise!='BP',genome=genome)
##      if(summarise=='Region'){
##        out_for[i,] = out_for[i,] + colSums(tmp$forw)
##        out_rev[i,] = out_rev[i,] + colSums(tmp$revr)
##      }else if(summarise=='GRange'){
##        reg_list[[i]]$forw = tmp$forw
##        reg_list[[i]]$revr = tmp$revr
##      }else{
##        reg_list[[i]]$forw = lapply(as.list(data.frame(t(tmp$forw))),`[[`,1)
##        reg_list[[i]]$revr = lapply(as.list(data.frame(t(tmp$revr))),`[[`,1)
##      }
##      #Now mark it as done
##      todo[i]=FALSE
##    }
##  }
##  #Process one chromosome at a time, unless we did them all quickly
##  if(any(todo)){
##    for(chr in chrs){
##      cat(sprintf('Processing chromosome %s\n',chr))
##      cat('Loading sequence.\n')
##      #Work out if we should get the whole chromosome, or just the relevant subset
##      chr_frac = sum(width(cvg[[chr]]))/clens[chr]
##      if(chr_frac > max_genome_frac){
##        genome_seq = getSeq(genome,prefix_chr(chr))
##      }else{
##        genome_seq = getSeq(genome,renameSeqlevels(cvg[[chr]],chrFormatMap))
##      }
##      #Mask anything that isn't pure ACGT in the window
##      drop_bases = paste0(grep('[^ACGT]',alphabet(genome_seq),value=TRUE),collapse='')
##      drop_bases_id = paste0(grep('[^ACGT]',alphabet(genome_seq),value=TRUE),collapse='|')
##      #Get the GC content in sliding windows from the sequence
##      cat('Binning GC.\n')
##      if(chr_frac > max_genome_frac) {
##        tmp = letterFrequencyInSlidingView(genome_seq,c('GC',drop_bases),view.width=wlen)
##        cat('Re-formatting.\n')
##        rm(genome_seq)
##        #This chomps memory for some reason.
##        gc()
##        tmp = Rle(ifelse(tmp[,drop_bases_id]==0,tmp[,'G|C'],NA))
##        gc()
##        cat('Counting at targets.\n')
##        #Get the counts at each position
##        for(j in seq_along(reg_list)){
##          if(!todo[j]){
##            next
##          }
##          #Which regions do we need to get
##          o = which(fix_chr(as.character(seqnames(reg_list[[j]])))==chr)
##          #Don't count if this is a reverse strand position
##          oof = o[as.character(strand(reg_list[[j]][o]))!='-']
##          oor = o[as.character(strand(reg_list[[j]][o]))!='+']
##          if(summarise=='Region'){
##            #Get the windows with the right offsets
##            out_for[j,] = out_for[j,] + ltable(tmp[shift(ranges(reg_list[[j]][oof]),lstrip)],levels=0:wlen)
##            out_rev[j,] = out_rev[j,] + ltable(tmp[shift(ranges(reg_list[[j]][oor]),-lstrip-wlen+1)],levels=0:wlen)
##          }else{
##            #Get windows at each target location
##            gcf = lapply(shift(ranges(reg_list[[j]][oof]),lstrip),function(e) tmp[e])
##            gcr = lapply(shift(ranges(reg_list[[j]][oor]),-lstrip-wlen+1),function(e) tmp[e])
##            if(summarise=='GRange'){
##              #Summarise at GRange level
##              reg_list[[j]]$forw[o] = lapply(gcf,ltable,levels=0:wlen)
##              reg_list[[j]]$revr[o] = lapply(gcr,ltable,levels=0:wlen)
##            }else{
##              reg_list[[j]]$forw[o] = lapply(gcf,as.numeric)
##              reg_list[[j]]$revr[o] = lapply(gcr,as.numeric)
##            }
##          }
##        }
##      }else{
##        tmp = lapply(genome_seq,letterFrequencyInSlidingView,letters=c('GC',drop_bases),view.width=wlen)
##        cat('Re-formatting.\n')
##        gc()
##        #Make it a RleList, subsetting is then easier
##        tmp = RleList(lapply(tmp,function(e) Rle(ifelse(e[,drop_bases_id]==0,e[,'G|C'],NA))))
##        names(tmp) = seq_along(tmp)
##        gc()
##        cat('Counting at targets.\n')
##        for(j in seq_along(reg_list)){
##          if(!todo[j]){
##            next
##          }
##          #Which regions do we need to get
##          o = which(fix_chr(as.character(seqnames(reg_list[[j]])))==chr)
##          #Get the target regions.  Need to be careful about chromosome names here.
##          fnomMap = setNames(fix_chr(seqlevels(reg_list[[j]])),seqlevels(reg_list[[j]]))
##          ot = findOverlaps(renameSeqlevels(reg_list[[j]][o],fnomMap),cvg[[chr]])
##          shits = subjectHits(ot)
##          if(any(sort(queryHits(ot))!=seq(length(o)))){
##            stop('This should not happen!')
##          }
##          #Work out the right offset for each region.
##          #-start+1 gets us the window at the positions in reg_list.  Then need an extra offset to get the window we actually want.
##          for_pos = shift(ranges(reg_list[[j]][o]),-start(cvg[[chr]][shits])+1+lstrip)
##          rev_pos = shift(ranges(reg_list[[j]][o]),-start(cvg[[chr]][shits])+1-lstrip-wlen+1)
##          #Order them sensibly to allow fast allocation
##          fro = order(as.numeric(shits),start(for_pos),end(for_pos))
##          rro = order(as.numeric(shits),start(rev_pos),end(rev_pos))
##          #Convert them to ranges lists for easy subsetting of GC counts
##          for_pos = GRanges(shits[fro],for_pos[fro])
##          rev_pos = GRanges(shits[rro],rev_pos[rro])
##          #Don't count if this is a reverse strand position
##          oof = as.character(strand(reg_list[[j]][o[fro]]))!='-'
##          oor = as.character(strand(reg_list[[j]][o[rro]]))!='+'
##          if(summarise=='Region'){
##            out_for[j,] = out_for[j,] + ltable(tmp[as(for_pos[oof],'RangesList')]@unlistData,levels=0:wlen)
##            out_rev[j,] = out_rev[j,] + ltable(tmp[as(rev_pos[oor],'RangesList')]@unlistData,levels=0:wlen)
##          }else{
##            #The fast and hopefully never wrong way
##            gcf = as.list(split(tmp[as(for_pos[oof],'RangesList')]@unlistData,Rle(seq_along(o[fro][oof]),width(reg_list[[j]][o[fro]][oof]))))
##            gcr = as.list(split(tmp[as(rev_pos[oor],'RangesList')]@unlistData,Rle(seq_along(o[rro][oor]),width(reg_list[[j]][o[rro]][oor]))))
##            if(summarise=='GRange'){
##              reg_list[[j]]$forw[o[fro][oof]] = lapply(gcf,ltable,levels=0:wlen)
##              reg_list[[j]]$revr[o[rro][oor]] = lapply(gcr,ltable,levels=0:wlen)
##            }else{
##              reg_list[[j]]$forw[o[fro][oof]] = lapply(gcf,as.numeric)
##              reg_list[[j]]$revr[o[rro][oor]] = lapply(gcr,as.numeric)
##            }
##          }
##        }
##      }
##    }
##  }
##  if(summarise=='Region'){
##    return(list(forw=out_for,revr=out_rev))
##  }else{
##    return(reg_list)
##  }
##}

##DEPRECATED: Now part of TCN function.  
###' Get's the GC content in windows anchored at the small number of genomic locations given
###'
###' Unlike get_GC_counts, this function is optimised to get the GC content at a small number of genomic locations.  It should probably never be directly called, instead use get_GC_counts which will call this function when appropriate.
###'
###' @param regions The GRanges object giving the positions at which GC should be calculated.
###' @param isize The insert size to use at each location (including lstrip and rstrip).
###' @param lstrip The number of bases skipped at the 5' end of the fragment.
###' @param rstrip The number of bases skipped at the 3' end of the fragment.
###' @param summarise If True, The total count of each GC value within each region is returned.  Otherwise, the raw GC values at each location are.
###' @param genome Either NULL or the output of bin_GC_genome.  If bin_GC_genome output, this is used to calculate GC at desired locations more quickly.
###' @return A list containing forward and reverse strand counts under as "forw" and "rev".  If summarise is False, each will contain a list of the same length as regions, named by the regions, containing the GC count based at each position in regions in turn.  If summarised, each entry is a matrix with columns giving the number of GCs and rows giving the different regions.
##get_GC_counts_for_region = function(regions,isize,lstrip=3,rstrip=3,summarise=TRUE,genome=Hsapiens) {
##  wlen = isize-rstrip-lstrip
##  #Build a new adjusted window object that gets the sequence we need
##  windows = GRanges(prefix_chr(seqnames(regions)),IRanges(pmax(1,start(regions)-(isize-rstrip-1)),pmin(seqlengths(genome)[prefix_chr(seqnames(regions))],end(regions)+(isize-rstrip-1))))
##  #Now get the sequence
##  tgt_seq = getSeq(genome,windows)
##  #Count GCs and bases that are bad
##  drop_bases = paste0(grep('[^ACGT]',alphabet(tgt_seq),value=TRUE),collapse='')
##  drop_bases_id = paste0(grep('[^ACGT]',alphabet(tgt_seq),value=TRUE),collapse='|')
##  tgt_gc = lapply(tgt_seq,letterFrequencyInSlidingView,c('GC',drop_bases),view.width=wlen)
##  tgt_gc = lapply(tgt_gc,function(e) ifelse(e[,drop_bases_id]==0,e[,'G|C'],NA))
##  #Work out what the left and right offset from the regions we care about was (after truncating and sequence fetch)
##  loff = start(regions)-start(windows)
##  roff = end(windows)-end(regions)
##  #Get all the forward strand GC values we can
##  #Need to strip off the bits we added, but also the bits that refer to regions outside our target
##  tgt_gc_for = lapply(seq_along(tgt_gc),function(i) c(tgt_gc[[i]][-(1:(loff[i]+lstrip))],rep(NA,(isize-rstrip-1)-roff[i])))
##  #The +1 is because to get the last n of vector of length N is (N-n+1):N
##  tgt_gc_rev = lapply(seq_along(tgt_gc),function(i) c(rep(NA,(isize-rstrip-1)-loff[i]),tgt_gc[[i]][-((length(tgt_gc[[i]])-roff[i]-lstrip+1):length(tgt_gc[[i]]))]))
##  names(tgt_gc_for) = as.character(regions)
##  names(tgt_gc_rev) = as.character(regions)
##  if(summarise){
##    #Count the members in each segment
##    seg_counts_for = t(sapply(tgt_gc_for,ltable,levels=0:wlen))
##    seg_counts_rev = t(sapply(tgt_gc_rev,ltable,levels=0:wlen))
##    return(list(forw=seg_counts_for,revr=seg_counts_rev))
##  }
##  #Just return the numbers
##  return(list(forw=tgt_gc_for,revr=tgt_gc_rev))
##}


##DEPRECATED: Now part of TCN function.  
###' Get's the coverage at a specified location.  Should only be used for small regions
###'
###' @param bfile The BAM file to load counts from.
###' @param tgts The GRanges at which to calculate the coverage.
###' @return The tgts object with a obsCov column containing the coverage at each location.
##calculateCoverage = function(bfile,tgts) {
##  require(GenomicFiles)
##  #Filter out bad reads
##  #Keep only the first mate of a read, discard anything even mildly suspicious
##  flags = scanBamFlag(isPaired = TRUE, isProperPair = TRUE, isUnmappedQuery = FALSE,
##         hasUnmappedMate = FALSE, isMinusStrand = NA, isMateMinusStrand = NA,
##         isFirstMateRead = TRUE, isSecondMateRead = FALSE, 
##         isSecondaryAlignment = FALSE, isNotPassingQualityControls = FALSE,
##         isDuplicate = NA)
##  bf = open(BamFile(bfile,yieldSize=10000000))
##  #Get the positions of any read in the target region
##  param = ScanBamParam(flags,what=c('rname','pos','mpos','qwidth','isize'))
##  yield = function(e) scanBam(e,param=param)
##  map = function(e) {
##    pos = unlist(lapply(e,`[[`,'pos'),use.names=FALSE)
##    mpos = unlist(lapply(e,`[[`,'mpos'),use.names=FALSE)
##    qwidth = unlist(lapply(e,`[[`,'qwidth'),use.names=FALSE)
##    isize = unlist(lapply(e,`[[`,'isize'),use.names=FALSE)
##    strand = pos < mpos
##    tmp = GRanges(unlist(lapply(e,`[[`,'rname'),use.names=FALSE),
##            IRanges(ifelse(strand,pos,mpos),width=abs(isize)),
##            strand = ifelse(strand,'+','-'),
##            readLen = qwidth,
##            insert = isize
##            )
##    subsetByOverlaps(tmp,tgts)
##  }
##  #Keep just the count in each region
##  reducer = function(x,y,...) {
##    c(x,y)
##    ##Is this the first time?
##    #if(inherits(x,'GenomicRanges')){
##    #  x = coverage(x)
##    #}
##    ##Get the common set of sequences
##    #chrs = sort(unique(c(seqlevels(y),names(x))))
##    ##Add in extra chromosomes to both elements
##    #seqlevels(y) = chrs
##    #for(chr in chrs[!(chrs %in% names(x))]){
##    #  x[[chr]] = Rle(integer())
##    #}
##    ##Calculate new coverage
##    #y = coverage(y)
##    ##Combine, ordering appropriately
##    #x[chrs]+y[chrs]
##  }
##  done = function(e) length(e[[1]][[1]]) == 0L
##  #Actually load a subset of the reads
##  cat('Calculating coverage from reads.\n')
##  cov = coverage(reduceByYield(bf,yield,map,reducer,done))
##  close(bf)
##  tgts$obsCov = as.list(cov[tgts])
##  return(tgts)
##}


##DEPRECATED: Now part of TCN function.  
###' Predicts the coverage at the specified location based on GC content.
###'
###' @param tgts The GRanges object containing the regions at which we want predicted coverage.
###' @param gcBias The gcBias object for this sample.
###' @return The tgts object with a predCov column containing the predicted coverage.
##predict_coverage_from_GC = function(tgts,gcBias){
##  #Adjust the start/stop so that we get the whole coverage
##  ltgts = tgts
##  start(ltgts) = pmax(1,start(ltgts)-gcBias$isize+1)
##  end(ltgts) = end(ltgts)+gcBias$isize-1
##  #Reduce it, so we don't get overlapping region
##  ltgts = reduce(ltgts)
##  #Get the predicted read number at each location
##  tmp = predict_counts_from_GC(ltgts,gcBias,'BP')
##  #Construct GRanges for each "read"
##  reads = GRanges(rep(seqnames(tmp),width(tmp)),
##                IRanges(start=rep(start(tmp),width(tmp)) + unlist(lapply(width(tmp),seq),use.names=FALSE)-1,
##                        width=gcBias$isize),
##                cov = unlist(tmp$forw,use.names=FALSE)
##                )
##  reads = c(reads,GRanges(rep(seqnames(tmp),width(tmp)),
##                IRanges(end = rep(start(tmp),width(tmp)) + unlist(lapply(width(tmp),seq),use.names=FALSE)-1,
##                        width=gcBias$isize),
##                cov = unlist(tmp$revr,use.names=FALSE)
##                ))
## 
##  #Get the coverage object
##  cov = coverage(reads,weight=reads$cov)
##  #Finally, extract the coverage at the regions we wanted in the first instance
##  tgts$predCov = as.list(cov[tgts])
##  return(tgts)
##}


##DEPRECATED: Now part of TCN function.  
###' Predict number of reads in ranges from GC content
###'
###' Uses the GC bias correction given to calculate how many reads to expect based on the GC content in windows anchored at locations in the target regions specified.
###'
###' @param tgts A GRanges object giving the locations at which we should predict counts.
###' @param gcBias The gcBias object for this sample.
###' @param summarise Should we return the number of predicted counts per base pair ('BP'), per entry in tgts ('GRange') or in total ('Region').
###' @return If summarise=='Region' then a list giving the predicted number of forward and reverse counts.  Otherwise returns the tgts object with columns forw and revr giving the predicted counts for that particular range (if summarise is 'GRange') or at each base pair in order (if summarise is 'BP').
##predict_counts_from_GC = function(tgts,gcBias,summarise=c('Region','GRange','BP')) {
##  summarise = match.arg(summarise)
##  tmp = get_GC_counts(list(tgts),gcBias$isize,gcBias$lstrip,gcBias$rstrip,summarise=summarise)
##  if(summarise=='Region'){
##    return(list(forw=sum(gcBias$forw['rate',]*tmp$forw[1,]),revr=sum(gcBias$revr['rate',]*tmp$revr[1,])))
##  }else if(summarise=='GRange'){
##    tgts$forw = apply(tmp[[1]]$forw,1,function(e) {sum(gcBias$forw['rate',] * e)})
##    tgts$revr = apply(tmp[[1]]$revr,1,function(e) {sum(gcBias$revr['rate',] * e)})
##  }else{
##    tgts$forw = lapply(tmp[[1]]$forw,function(e) {gcBias$forw['rate',][as.integer(e+1)]})
##    tgts$revr = lapply(tmp[[1]]$revr,function(e) {gcBias$revr['rate',][as.integer(e+1)]})
##  }
##  return(tgts)
##}


####
###predict_counts_from_GC = function(tgts,gcBias,returnAllCounts=FALSE) {
###  wlen = gcBias$isize-gcBias$rstrip-gcBias$lstrip
###  genome = gcBias$genome
###  #Make the rate vector so we can quickly predict values
###  frat = Rle(gcBias$forw['rate',])
###  rrat = Rle(gcBias$revr['rate',])
###  #Get the chromosome lengths
###  clens = seqlengths(genome)
###  names(clens) = fix_chr(names(clens))
###  #Get the list of chromosomes we need to look at
###  chrs = unique(fix_chr(as.character(seqnames(tgts))))
###  #Order them sensibly
###  chrs = chrs[order(num_chr(chrs))]
###  #The remapping of things to UCSC genome chr format 
###  chrFormatMap = setNames(prefix_chr(chrs),chrs)
###  #Check that we haven't asked for bad regions
###  if(any(start(tgts)<gcBias$lstrip+wlen) | any(end(tgts) > clens[fix_chr(as.character(seqnames(tgts)))]-gcBias$lstrip-wlen+1)){
###    stop('Location requested without any valid window.')
###  }
###  #Make a coverage object so we know where we need data
###  cvg = reduce(tgts)
###  cvg = GRanges(fix_chr(as.character(seqnames(cvg))),ranges(cvg))
###  #Add flank that we'll need
###  start(cvg) = pmax(1,start(cvg)-gcBias$lstrip-wlen+1)
###  end(cvg) = pmin(clens[fix_chr(as.character(seqnames(cvg)))],end(cvg)+gcBias$lstrip+wlen-1)
###  #Reduce it again
###  cvg = reduce(cvg)
###  #Split by (fixed up) chromosome
###  cvg = split(cvg,as.character(seqnames(cvg)))
###  #Process one chromosome at a time.  Store the predicted counts as extra columns in tgts
###  #Add columns to targets to store output
###  tgts$predCntsFor = rep(as.numeric(NA),length(tgts))
###  tgts$predCntsRev = rep(as.numeric(NA),length(tgts))
###  if(returnAllCounts){
###    counts_for = vector('list',length(tgts))
###    counts_rev = vector('list',length(tgts))
###  }
###  for(chr in chrs){
###    cat(sprintf('Processing chromosome %s\n',chr))
###    cat('Loading sequence.\n')
###    #Work out if we should get the whole chromosome, or just the relevant subset
###    chr_frac = sum(width(cvg[[chr]]))/clens[chr]
###    if(chr_frac > .5){
###      genome_seq = getSeq(genome,prefix_chr(chr))
###    }else{
###      genome_seq = getSeq(genome,renameSeqlevels(cvg[[chr]],chrFormatMap))
###    }
###    #This gives us the entire genomes GC content in sliding windows, in Rle format.
###    drop_bases = paste0(grep('[^ACGT]',alphabet(genome_seq),value=TRUE),collapse='')
###    drop_bases_id = paste0(grep('[^ACGT]',alphabet(genome_seq),value=TRUE),collapse='|')
###    #Mask anything that isn't pure ACGT in the window
###    cat('Binning GC.\n')
###    if(chr_frac > 0.5) {
###      tmp = letterFrequencyInSlidingView(genome_seq,c('GC',drop_bases),view.width=wlen)
###      cat('Re-formatting.\n')
###      #This chomps memory for some reason.
###      gc()
###      #Predict no counts if the window contains non ACGT
###      tmp = Rle(ifelse(tmp[,drop_bases_id]==0,tmp[,'G|C'],0))
###      gc()
###      cat('Counting at targets.\n')
###      #Get the GC at each position
###      o = which(fix_chr(as.character(seqnames(tgts)))==chr)
###      #Predict forward reads from GC
###      oo = o[as.character(strand(tgts[o]))!='-']
###      tgts[oo]$predCntsFor = aggregate(frat[tmp[shift(ranges(tgts[oo]),gcBias$lstrip)]+1],IRanges(start=cumsum(c(0,width(tgts[oo])[-length(oo)]))+1,width=width(tgts[oo])),sum)
###      #tgts[oo]$predCntsFor = sapply(split(frat[tmp[shift(ranges(tgts[oo]),gcBias$lstrip)]+1],Rle(seq_along(oo),width(tgts[oo]))),sum)
###      oo = o[as.character(strand(tgts[o]))!='+']
###      tgts[oo]$predCntsRev = aggregate(rrat[tmp[shift(ranges(tgts[oo]),-gcBias$lstrip-wlen+1)]+1],IRanges(start=cumsum(c(0,width(tgts[oo])[-length(oo)]))+1,width=width(tgts[oo])),sum)
###      #tgts[oo]$predCntsRev = sapply(split(rrat[tmp[shift(ranges(tgts[oo]),-gcBias$lstrip-wlen+1)]+1],Rle(seq_along(oo),width(tgts[oo]))),sum)
###    }else{
###      tmp = lapply(genome_seq,letterFrequencyInSlidingView,letters=c('GC',drop_bases),view.width=wlen)
###      cat('Re-formatting.\n')
###      gc()
###      #Make it a RleList, subsetting is then easier
###      tmp = RleList(lapply(tmp,function(e) Rle(ifelse(e[,drop_bases_id]==0,e[,'G|C'],0))))
###      names(tmp) = seq_along(tmp)
###      gc()
###      cat('Counting at targets.\n')
###      #Which regions do we need to get
###      o = which(fix_chr(as.character(seqnames(tgts)))==chr)
###      #Get the target regions.  Need to be careful about chromosome names here.
###      fnomMap = setNames(fix_chr(seqlevels(tgts)),seqlevels(tgts))
###      ot = findOverlaps(renameSeqlevels(tgts[o],fnomMap),cvg[[chr]])
###      shits = subjectHits(ot)
###      if(any(sort(queryHits(ot))!=seq(length(o)))){
###        stop('This should not happen!')
###      }
###      #Work out the right offset for each region.
###      #-start+1 gets us the window at the positions in tgts  Then need an extra offset to get the window we actually want.
###      for_pos = shift(ranges(tgts[o]),-start(cvg[[chr]][shits])+1+gcBias$lstrip)
###      rev_pos = shift(ranges(tgts[o]),-start(cvg[[chr]][shits])+1-gcBias$lstrip-wlen+1)
###      #Order them sensibly
###      fro = order(as.numeric(shits),start(for_pos),end(for_pos))
###      rro = order(as.numeric(shits),start(rev_pos),end(rev_pos))
###      #Convert them to ranges lists for easy subsetting of GC counts
###      for_pos = GRanges(shits[fro],for_pos[fro])
###      rev_pos = GRanges(shits[rro],rev_pos[rro])
###      #Don't count if this is a reverse strand position
###      oo = as.character(strand(tgts[o][fro]))!='-'
###      #Predict forward counts from GC
###      #Do this in one shot.  Might be prone to error (although I've tested it without finding an issue).  Should be faster than the loop below...
###      #tgts[o[fro][oo]]$predCntsFor = sapply(split(frat[tmp[as(for_pos[oo],'RangesList')]@unlistData+1],Rle(seq_along(which(oo)),width(tgts[o[fro][oo]]))),sum)
###      #The slow and steady way...
###      for(i in seq_along(which(oo))){
###        tgts[o[fro][oo][i]]$predCntsFor = sum(frat[tmp[[as.character(seqnames(for_pos[oo][i]))]][ranges(for_pos[oo][i])]+1])
###        #Store the counts if we need to
###        if(returnAllCounts){
###          counts_for[[o[fro][oo][i]]] = tmp[[as.character(seqnames(for_pos[oo][i]))]][ranges(for_pos[oo][i])]
###        }
###      }
###      ###A quick sanity check that we're storing things in the right places
###      ##if(any(width(tgts[o[fro][oo]])!=width(for_pos[oo]))){
###      ##  stop("I've made a huge mistake.")
###      ##}
###      #Repeat for other strand
###      oo = as.character(strand(tgts[o][rro]))!='+'
###      #tgts[o[rro][oo]]$predCntsRev = sapply(split(frat[tmp[as(rev_pos[oo],'RangesList')]@unlistData+1],Rle(seq_along(which(oo)),width(tgts[o[rro][oo]]))),sum)
###      #The slow and steady way...
###      for(i in seq_along(which(oo))){
###        tgts[o[rro][oo][i]]$predCntsRev = sum(frat[tmp[[as.character(seqnames(rev_pos[oo][i]))]][ranges(rev_pos[oo][i])]+1])
###        #Store the counts if we need to
###        if(returnAllCounts){
###          counts_rev[[o[rro][oo][i]]] = tmp[[as.character(seqnames(rev_pos[oo][i]))]][ranges(rev_pos[oo][i])]
###        }
###      }
###      ###A quick sanity check that we're storing things in the right places
###      ##if(any(width(tgts[o[rro][oo]])!=width(rev_pos[oo]))){
###      ##  stop("I've made a huge mistake.")
###      ##}
###    }
###  }
###  if(returnAllCounts){
###    return(list(counts_for=counts_for,counts_rev=counts_rev))
###  }
###  return(tgts)
###}


