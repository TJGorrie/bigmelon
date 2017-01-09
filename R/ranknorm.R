quickquan <- function(gds, node, onetwo = NULL, rank = FALSE, new.node = NULL, perc = 1, ...){
    # getquantiles, the REALLY fast way to quickly normalize
    if(rank & is.null(new.node)) stop('Generating Ranks, Please specify a new.node name')
    if(perc > 1 | perc < 0) stop('perc must be bounded between 0 and 1')
    datnod <- node
    if(length(datnod) == 1) datnod <- index.gdsn(gds, as.character(datnod))
    dim <- objdesp.gdsn(datnod)$dim
    nc <- 1:dim[2]
    nc <- sample(nc, length(nc)*perc, replace=FALSE)
    roll <- rep(0, dim[1])
    if(is.null(onetwo)){
        # Initialising rolling sum
        nobs <- rep(dim[1], length(nc))
        i <- (0:(dim[1]-1))/(dim[1]-1)
        # Sorting and 'growing' rolling sum.
        if(rank){
            n.t <- add.gdsn(gds, as.character(new.node), storage = 'float64', valdim=c(dim[1],0), val = NULL, replace = TRUE)
            n.a <- add.gdsn(gds, paste0('isna', new.node), storage = 'int8', valdim=c(dim[1],0), val = NULL, replace = TRUE, visible = FALSE)
        }
        for(x in 1:dim[2]){
            val <- readex.gdsn(datnod, sel = list(NULL, x))
            if(x %in% nc){ # If Perc?
            S <- rep(NA, dim[1])
            si <- sort(val, method = "quick", index.return = TRUE)
            nobsj <- length(si$x)
            # Used later for NA's
            if(nobsj < dim[1]){
                nobs[x] <- nobsj
                S <- approx((0:(nobsj-1))/(nobsj-1), si$x, i, ties = "ordered")$y
            } else {
                S <- si$x
            }
            roll <- roll + S
            } # If Perc?
            if(rank){
                ranks <- rank(val)
                #ranks[onetwo=='I'] <- rank(val[onetwo=='I'])
                #ranks[onetwo=='II'] <- rank(val[onetwo=='II'])
                append.gdsn(n.t, ranks)
                append.gdsn(n.a, as.numeric(is.na(val)))
            }
        }  
        # Calculating 'rowMeans'
        rm <- roll/length(nc)
        if(rank){
            put.attr.gdsn(n.t, 'ranked', val = TRUE)
            put.attr.gdsn(n.t, 'is.na', val = paste0('isna',new.node))
            inter <- i
            #inter[onetwo=='I'] <- i.i
            #inter[onetwo=='II'] <- ii.i
            put.attr.gdsn(n.t, 'inter', val = inter)
            put.attr.gdsn(n.t, 'quantiles', val = rm)
            #put.attr.gdsn(n.t, 'onetwo', val = onetwo)
        }
#
#
#
    } else { # IF DESIGN IS SUPPLIED
#
#
#
        roll <- rep(0, dim[1])
        inum    <- sum(onetwo=='I')
        i.nobs  <- rep(inum, length(nc))
        i.i     <- (0:(inum-1))/(inum-1)
        # Type II
        iinum   <- sum(onetwo=='II')
        ii.nobs <- rep(iinum, length(nc))
        ii.i    <- (0:(iinum-1))/(iinum-1)
        if(rank){
            n.t <- add.gdsn(gds, as.character(new.node), storage = 'float64', valdim=c(dim[1],0), val = NULL, replace = TRUE)
            n.a <- add.gdsn(gds, paste0('isna', new.node), storage = 'int8', valdim=c(dim[1],0), val = NULL, replace = TRUE, visible = FALSE)
        }
        # Sorting + Rolling Sum 
        for(x in 1:dim[2]){
            val   <- readex.gdsn(datnod, sel=list(NULL, x)) 
            if(x %in% nc){ # If perc
            # Type I
            i.S      <- rep(NA, inum)
            i.si     <- sort(val[onetwo=='I' ], method = "quick", index.return = TRUE)
            i.nobsj  <- length(i.si$x)
            if(i.nobsj < inum){
                i.nobs[x] <- i.nobsj
                i.S <- approx((0:(i.nobsj-1))/(i.nobsj-1), i.si$x, i.i, ties = "ordered")$y
            } else {
                i.S <- i.si$x
            }
            roll[onetwo=='I']  <- roll[onetwo=='I' ] + i.S

            # Type II
            ii.S     <- rep(NA, iinum)
            ii.si    <- sort(val[onetwo=='II'], method = "quick", index.return = TRUE)
            ii.nobsj <- length(ii.si$x)

            if(ii.nobsj < iinum){
                ii.nobs[x] <- ii.nobsj
                ii.S <- approx((0:(ii.nobsj-1))/(ii.nobsj-1), ii.si$x, ii.i, ties = "ordered")$y
            } else {
                ii.S <- ii.si$x
            }
            roll[onetwo=='II'] <- roll[onetwo=='II'] + ii.S
            } # End perc
            if(rank){
                ranks <- rep(NA, dim[1])
                ranks[onetwo=='I'] <- rank(val[onetwo=='I'])
                ranks[onetwo=='II'] <- rank(val[onetwo=='II'])
                append.gdsn(n.t, ranks)
                append.gdsn(n.a, as.numeric(is.na(val)))
            }
        }
        # rowmeans (Generate quantiles)
        rm <- roll/length(nc)
        if(rank){
            put.attr.gdsn(n.t, 'ranked', val = TRUE)
            put.attr.gdsn(n.t, 'is.na', val = paste0('isna',new.node))
            inter <- rep(NA, dim[1])
            inter[onetwo=='I'] <- i.i
            inter[onetwo=='II'] <- ii.i
            put.attr.gdsn(n.t, 'inter', val = inter)
            put.attr.gdsn(n.t, 'quantiles', val = rm)
            put.attr.gdsn(n.t, 'onetwo', val = onetwo)
        }
    }
    # If not storing ranks return quantiles.
    if(!rank) return(rm) else return(0)
}

dasenrank <- function(gds, mns, uns, onetwo, roco, calcbeta = NULL, ...){ # {{{
    # Assuming that mns and uns are 1 element character strings!not gdsn.nodes
    if(length(mns) == 1) mns <- index.gdsn(gds, mns)
    if(length(uns) == 1) uns <- index.gdsn(gds, uns)
    if(length(onetwo) == 1)  onetwo <- read.gdsn(index.gdsn(gds, onetwo))
    if(class(onetwo) == 'gdsn.class') onetwo <- read.gdsn(onetwo)
    f <- createfn.gds("temp.gds", allow.duplicate = TRUE)
    dim <- objdesp.gdsn(mns)$dim
    # NORMALIZING
    dfsfit.gdsn(f, targetnode = mns, roco = roco, newnode = "mnsc",
                onetwo = onetwo)
    dfsfit.gdsn(f, targetnode = uns, roco = NULL, newnode = "unsc",
                onetwo = onetwo)
    # Get Rank + Quantiles
    quickquan(gds, index.gdsn(f,'mnsc'), onetwo = onetwo, rank=TRUE, new.node = 'mnsrank', perc=1) 
    quickquan(gds, index.gdsn(f,'unsc'), onetwo = onetwo, rank=TRUE, new.node = 'unsrank', perc=1)
    #getrankandquantiles(gds, index.gdsn(f, 'mnsc'), onetwo, 'mnsrank')
    #getrankandquantiles(gds, index.gdsn(f, 'unsc'), onetwo, 'unsrank')
    # COMPLETE Return nothing
    if(is.null(calcbeta)){
        message('Run \'computebetas\' to calculate betas!')
    } else { 
        computebetas(gds, calcbeta, 'mnsrank', 'unsrank', fudge = 100)
    }
    closefn.gds(f)
    unlink('temp.gds', force = TRUE)
} # }}}

computebetas <- function(gds, new.node, mns, uns, fudge = 100, ...){ # {{{
    if(length(mns) == 1) mns <- index.gdsn(gds, mns)
    if(length(uns) == 1) uns <- index.gdsn(gds, uns)    
    dim <- objdesp.gdsn(mns)$dim
    n.t <- add.gdsn(gds, new.node, storage = 'float64', valdim=c(dim[1],0), val = NULL, replace = TRUE)
    for(x in 1:dim[2]){
        meth <- mns[, x, name = FALSE]
        unmeth <- uns[, x, name = FALSE]
        beta <- meth/(meth + unmeth + fudge)
        append.gdsn(n.t, beta)
    }
} # }}}

'[.gdsn.class' <- function(x, i, j, name = TRUE, drop = TRUE){ # {{{
  # Method of subsetting gdsn.class objects w.o reading in entire object.
  # Wrapper for readex.gdsn
  # arg: "name": Will point towards "fData/Probe_ID" and "pData/barcode" for row and col names. (by default)
  #              name = F will enable for faster indexing.
  # TODO: [] method? a.k.a vector method.
  # TODO: Warnings, if any.
  base <- getfolder.gdsn(x)
  dim <- objdesp.gdsn(x)$dim
  
  # x[ , ]
  if(missing(i) & missing(j)){ # {{{
    mat <- read.gdsn(x)
    if(name){
      rownames(mat) <- read.gdsn(index.gdsn(base,
                                 read.gdsn(index.gdsn(base, "paths"))[1]))

      colnames(mat) <- read.gdsn(index.gdsn(base,
                                 read.gdsn(index.gdsn(base, "paths"))[2]))
    }
  } # }}}

  # x[ , j]
  if(missing(i) & !missing(j)){ # {{{
    j1 <- j
    if(is.character(j1)) j <- match(j1, read.gdsn(index.gdsn(base,
                                                  read.gdsn(index.gdsn(base, "paths"))[2]))) # ok
    if(is.logical(j1))   j <- (1:objdesp.gdsn(x)$dim[2])[j1] # ok

    ncol <- j
    mat <- as.matrix(readex.gdsn(x, sel = list(NULL, ncol)))
    if(length(j)==1) mat <- as.matrix(mat)
    if(name){
      rownames(mat) <- read.gdsn(index.gdsn(base,
                                 read.gdsn(index.gdsn(base, "paths"))[1]))
      colnames(mat) <- readex.gdsn(index.gdsn(base,
                                              read.gdsn(index.gdsn(base, "paths"))[2]), sel = ncol)
    }
  } # }}}

  # x[ i, ]
  if(!missing(i) & missing(j)){  # {{{
    i1 <- i
    if(is.character(i1)) i <- match(i1, read.gdsn(index.gdsn(base,
                                                  read.gdsn(index.gdsn(base, "paths"))[1]))) # ok
    if(is.logical(i1))   i <- (1:objdesp.gdsn(x)$dim[2])[i1] # ok
    nrow <- i
    if(length(i)==1){ # Calling a single row makes naming difficult this is a work-around.
      mat <- as.matrix(t(readex.gdsn(x, sel = list(nrow, NULL))))
    } else {
      mat <- as.matrix(readex.gdsn(x, sel = list(nrow, NULL)))
    }

    if(name){
      rownames(mat) <- readex.gdsn(index.gdsn(base,
                                              read.gdsn(index.gdsn(base, "paths"))[1]), sel = nrow)
      colnames(mat) <- read.gdsn(index.gdsn(base,
                                            read.gdsn(index.gdsn(base, "paths"))[2]))
    }

  } # }}}

  # x[ i, j] 
  if(!missing(i) & !missing(j)){# {{{

    i1 <- i
    if(is.character(i1)) i <- match(i1, read.gdsn(index.gdsn(base,
                                                  read.gdsn(index.gdsn(base, "paths"))[1]))) # ok
    if(is.logical(i1))   i <- (1:objdesp.gdsn(x)$dim[2])[i1] # ok

    j1 <- j
    if(is.character(j1)) j <- match(j1, read.gdsn(index.gdsn(base,
                                                             read.gdsn(index.gdsn(base, "paths"))[2]))) # ok
    if(is.logical(j1))   j <- (1:objdesp.gdsn(x)$dim[2])[j1] # ok

    nrow <- i
    ncol <- j
    mat <- readex.gdsn(x, sel = list(nrow, ncol))

    # i=1, j=1
    if( length(i) == 1 & length(j)  == 1) mat <- matrix(mat)
    # i>1, j=1
    if(!length(i) == 1 & length(j)  == 1) mat <- matrix(mat)
    # i=1, j>1
    if( length(i) == 1 & !length(j) == 1) mat <- t(matrix(mat)) 

    if(name){
      rownames(mat) <- readex.gdsn(index.gdsn(base,
                                   read.gdsn(index.gdsn(base, "paths"))[1]), sel = nrow)
      colnames(mat) <- readex.gdsn(index.gdsn(base,
                                              read.gdsn(index.gdsn(base, "paths"))[2]), sel = ncol)   
    }
  } # }}}

  ranked <- get.attr.gdsn(x)[['ranked']]
  if(is.null(ranked)) ranked <- FALSE
  if(ranked){
      quantiles <- get.attr.gdsn(x)[['quantiles']]
      inter <- get.attr.gdsn(x)[['inter']]
      ot <- get.attr.gdsn(x)[['onetwo']]
      design <- FALSE
      if(!is.null(ot)) design <- TRUE 
      isna <- index.gdsn(base, get.attr.gdsn(x)[['is.na']])[i, j, name = TRUE, drop = FALSE]
      for(z in 1:ncol(isna)){
          acol <- isna[,z]
          re <- rep(NA, length(acol))
          if(design){
          re[ot=='I'&&!acol] <- approx(inter[ot=='I'], quantiles[ot=='I'], (mat[ot=='I'&&!acol,z]-1)/(sum(ot=='I')-1), ties = "ordered")$y
          re[ot=='II'&&!acol] <- approx(inter[ot=='II'], quantiles[ot=='II'], (mat[ot=='II'&&!acol,z]-1)/(sum(ot=='II')-1), ties = "ordered")$y
          } else {
          re[!acol] <- approx(inter, quantiles, (mat[!acol,z]-1)/(length(inter)-1), ties='ordered')$y
          }
          mat[,z] <- re
      }
  }  
  return(mat)
}


