# {{{ plan LS 2023

# if the gds representation is to fulfil the ambition of acting like a MethyLumiSet (or methylset)  object, methods should
# provide the same illusion that the data is a matrix.   Currently there isn't a dim() method for gds or gdsn classes.  there
# is a `[` method for gdsn, which is the main engine for extracting data.  The `[` method for gds just hands off to gdsn 
# (ie it provides the convenience of not having to extract eg the betas before subseting but still just subsets a node).

# for usability a key point will be the distinction between extracting data and irreversibly changing the file representation.  
# users will intuitively be reluctant to make irreversible changes, but this has to be balanced by the limited number of copies 
# of large data sets that it's practical to keep.  

# The workflow model should probably be that in the dataset assembly/qc phase you work with a read/write handle to he gds file 
# and do irreversible filtering and combining ops, and in the analysis phase you work with a readonly handle.  This would argue
# for methods whose behaviour depends on the handle, (and should have clear feedback, possibly asking for confirmation for 
# irreversible ops).

# additional possibilities:  functions could check if their return is being assigned and if so return a handle to a new gds file 
# and/or create a new file if a path is given (ie default path=NULL)

# one concept that we have discussed, but I don't think we have implemented, is infrastructure facilitating filters.  An extended gds 
# object could have one or more vectors representing desired subsets of rows and columns for particular purposes

# }}}

#{{{ chainsaw
#' Chainsaw -- modify gds file by subsetting all nodes
#'
#' Currently the `[` function for the gds.class objects used by bigmelon only subsets a 
#' single node. This function does more like what you would normally expect a subsetting function 
#' to do, it returns a subset of the entire object.  It may in future be a replacement for `[.gds.class`.
#'
#' This function is intended for use in the preprocessing and QC phase of a DNA methylation workflow.
#' For efficiency, bigmelon stores data in a file, and the gds.class object is a file handle. 
#' True to its name, chainsaw chops the underlying file, this is a side affect of the function and 
#' is not affected by assignment of the return value.
#'
#' @param gfile A gds.class object.
#' @param i Specifies rows (ie probes) in the desired subset, similar to behaviour of `[`
#' @param j Specifies columns (ie sampless) in the desired subset, similar to behaviour of `[`
#' @param v If true, spew many messages.
#' @param cleanup If true, run a cleanup function that can substantially reduce the file size.

#' @return  a gds.class object.  This is a handle to the same file that the gfile argument 
#' points to.  It's not generally useful to have two handles to the same file, but it may make code 
#' more readable.  In interactive use, if not assigned, the returned object is usefully pretty-printed.

chainsaw <- function( gfile, i='', j='', v=FALSE, cleanup=TRUE ){

   history.submitted <- as.character(Sys.time())

   # is this a gds file handle?
   stopifnot(inherits(gfile, "gds.class"))
   # is the file where it's supposed to be?
   stopifnot(file.exists(gfile$file))
   # is it open RW?
   stopifnot(!gfile$readonly)



   # any useful gds file in a bigmelon workflow will have a methylated and/or betas node 
   # take the dims from that and subset all similarly dimmed nodes accordingly using `[.gdsn.class`
   # <<<< note the different approach in the S4 row and column names functions >>>>
   # paths node points to "fData/Probe_ID" "pData/barcode" as row and colnames.


   #### TOUTDOUX ###############

   # handle i == ''  check
   # handle j == ''  check
   # handle fData    check
   # handle pData    check
   # use names       check
   # add history     check

   # there is a choice here of either going through the expected nodes by name (a),
   # or walking the tree and handling each node depending on what it is (b). 
   # b is attractive because this is not a fully sealed object, and users might reasonably
   # add nodes for transformed or normalised versions of the data. 

   # current vers is a compromise that subsets top level matrices -- 
   # and deals with fData, pData and history explicitly
   # nodes with 0 columns, 0 rows or names starting with 'tmp' are skipped.
   # (there are sometimes tmp_* nodes introduced by gdsfmt, haven't found details)
   
   dimes <- sapply(ls.gdsn(gfile), function(x) dim(index.gdsn(gfile,x)))
   if (v) message(length(dimes), " nodes")
   mates <- sapply(dimes,length) ==2
   if (v) message(sum(mates), " matrices")
   d <- dimes$methylated
   if (is.null(d)) d <- dimes$betas
   if (is.null(d)) stop ('no methylated or betas node\n')
   
   if (length(i) ==1 && i=='') i <- 1:d[1]
   if (length(j) ==1 && j=='') j <- 1:d[2]

   # char rownames
   if (is.character(i)){
     i <- match( i, rownames(gfile))
     if(!any(!is.na(i))) stop("no rownames match, giving up")
     if(any(is.na(i))) message("note some rownames don't match, filling with NAs")
  }   
   
   # char coalnames
   if (is.character(j)){
     j <- match( j, colnames(gfile))
     if(!any(!is.na(j))) stop("no samplenames match, giving up")
     if(any(is.na(j))) message("note some samplenames don't match, filling with NAs")
  }   

   vict <- dimes[mates] # 
   for (nod in names(vict)) {

      if( any(vict[[nod]] ==0) || grepl('^tmp', nod) ) next
      eye <- 1:vict[[nod]][1]
      jay <- 1:vict[[nod]][2]
      if (vict[[nod]][1] == d[1]  ) {
         eye <-i
         if (v) message('subsetting rows of ', nod)
      }
      if (vict[[nod]][2] == d[2]  ) {
         jay <-j
         if (v) message('subsetting columns of ', nod)
      }
   assign.gdsn(index.gdsn(gfile,nod), seldim=list(eye, jay))
   }  
   
   #### handle fData here, subset i
   for (col in ls.gdsn(index.gdsn(gfile,'fData/'))){
      nod <- index.gdsn(gfile,paste0('fData/',col))
      if (dim(nod) == d[1]){
         if (v) message('subsetting fData ', col)
         assign.gdsn(nod, seldim=list(i))
      }
   }
   
   #### handle pData here, subset j
   for (row in ls.gdsn(index.gdsn(gfile,'pData/'))){
      nod <- index.gdsn(gfile,paste0('pData/',row))
      if (dim(nod) == d[2]){
         if (v) message('subsetting pData ',row)
         assign.gdsn(nod, seldim=list(j))
      }
   }
   
   #### handle history here, add new dim
   history.finished <- as.character(Sys.time())
   history.command <- as.character(paste( 'subset ', length(i), ' probes and ', length(j),' samples.' ))
        h <- data.frame(submitted = history.submitted,
                        finished = history.finished,
                        command = history.command,
                        stringsAsFactors = FALSE)
        for(j in colnames(h)){
            h_coln <- j
            h_index_str <- paste("history/", j, sep = "")
            h_child_n <- index.gdsn(gfile, h_index_str, silent = TRUE)
            append.gdsn(h_child_n, val = h[ ,h_coln])
        }

   message( gfile$filename , " has been modified ")
   #### return gfile (will print new structure)
   if(cleanup) cleanup.gds(gfile$filename)
   gfile
}

#}}}

#{{{ dim.gds.class
#' dim.gds.class  S3 method returning dimensions of data represented by a gds file handle.

dim.gds.class <- function (gfile, v=FALSE){
   # note v=TRUE works only if called directly, generic dim doesn't allow 2 arrghs
   # is this a gds file handle
   # is the file where it's supposed to be

   # any useful gds file in a bigmelon workflow will have a methylated and/or betas node 
   dimes <- sapply(ls.gdsn(gfile), function(x) dim(index.gdsn(gfile,x)))
   if (v) message(length(dimes), " nodes")
   mates <- sapply(dimes,length) ==2
   if (v) message(sum(mates), " matrices")
   d <- dimes$methylated
   if (is.null(d)) d <- dimes$betas
   if (is.null(d)) stop ('no methylated or betas node\n')
   if (v) {
      message (d[1], ' probes, ', d[2], ' samples')
      invisible(d)
   }
   
   d

     
}   

#}}}

#{{{ dim.gdsn.class
# dim.gdsn.class S3 method returning dimensions of gds nodes if a matrix

dim.gdsn.class <- function (obj){

   # is obj a gds node
   stopifnot(inherits(obj, "gdsn.class"))
   # is the node a matrix (don't need to test, returns "", NULL, or length

   objdesp.gdsn(obj)$dim

}

#}}}

