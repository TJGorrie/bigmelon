quickquan <- function(gds, node, onetwo = NULL, rank = FALSE, new.node = NULL, perc = 1, ...){
    # When rank=FALSE quantiles are generated quite quickly.
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
            n.t <- add.gdsn(gds, as.character(new.node), storage = 'float64',
                            valdim=c(dim[1],0), val = NULL, replace = TRUE)
            n.a <- add.gdsn(gds, paste0('isna', new.node), storage = 'int8', 
                            valdim=c(dim[1],0), val = NULL, replace = TRUE, 
                            visible = FALSE)
        }
        for(x in 1:dim[2]){
            # If Rank = T read data 
            if(rank) val <- readex.gdsn(datnod, sel = list(NULL, x))
            if(x %in% nc){ # If Perc?
                # If Rank=F but x is in nc, read data! This will save time
                # by not reading data when unneeded!
                if(!rank) val <- readex.gdsn(datnod, sel = list(NULL, x))
                S <- rep(NA, dim[1])
                si <- sort(val, method = "quick", index.return = TRUE)
                nobsj <- length(si$x)
                # Used later for NA's
                if(nobsj < dim[1]){
                    nobs[x] <- nobsj
                    S <- approx((0:(nobsj-1))/(nobsj-1),
                                si$x, i, ties = "ordered")$y
                } else {
                    S <- si$x
                }
                roll <- roll + S
            } # If Perc?
            if(rank){
                ranks <- rank(val)
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
            put.attr.gdsn(n.t, 'inter', val = inter)
            put.attr.gdsn(n.t, 'quantiles', val = rm)
        }
    } else { # IF DESIGN IS SUPPLIED
        roll <- rep(0, dim[1])
        inum    <- sum(onetwo=='I')
        i.nobs  <- rep(inum, length(nc))
        i.i     <- (0:(inum-1))/(inum-1)
        # Type II
        iinum   <- sum(onetwo=='II')
        ii.nobs <- rep(iinum, length(nc))
        ii.i    <- (0:(iinum-1))/(iinum-1)
        if(rank){
            n.t <- add.gdsn(gds, as.character(new.node), storage = 'float64', 
                            valdim = c(dim[1],0), val = NULL, replace = TRUE)
            n.a <- add.gdsn(gds, paste0('isna', new.node), storage = 'int8', 
                            valdim = c(dim[1],0), val = NULL, replace = TRUE, 
                            visible = FALSE)
        }
        # Sorting + Rolling Sum 
        for(x in 1:dim[2]){
            if(rank) val   <- readex.gdsn(datnod, sel=list(NULL, x)) 
            if(x %in% nc){ # If perc
                if(!rank) val   <- readex.gdsn(datnod, sel=list(NULL, x)) 
                # Type I
                i.S      <- rep(NA, inum)
                i.si     <- sort(val[onetwo=='I' ], method = "quick", index.return = TRUE)
                i.nobsj  <- length(i.si$x)
                if(i.nobsj < inum){
                    i.nobs[x] <- i.nobsj
                    i.S <- approx((0:(i.nobsj-1))/(i.nobsj-1), 
                                    i.si$x, i.i, ties = "ordered")$y
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
                    ii.S <- approx((0:(ii.nobsj-1))/(ii.nobsj-1), 
                                    ii.si$x, ii.i, ties = "ordered")$y
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
    # If not storing ranks return generated quantiles.
    if(!rank) return(rm) else return(0)
}

dasenrank <- function(gds, mns, uns, onetwo, roco, calcbeta = NULL, ...){# {{{
    # Assuming that mns and uns are 1 element strings, not gdsn.class
    if(length(mns) == 1) mns <- index.gdsn(gds, mns)
    if(length(uns) == 1) uns <- index.gdsn(gds, uns)
    if(length(onetwo) == 1)  onetwo <- read.gdsn(index.gdsn(gds, onetwo))
    if(class(onetwo) == 'gdsn.class') onetwo <- read.gdsn(onetwo)
    f <- createfn.gds('temp.gds', allow.duplicate = TRUE)
    dim <- objdesp.gdsn(mns)$dim
    # NORMALIZING
    dfsfit.gdsn(f, targetnode = mns, roco = roco, newnode = "mnsc",
                onetwo = onetwo)
    dfsfit.gdsn(f, targetnode = uns, roco = NULL, newnode = "unsc",
                onetwo = onetwo)
    # Get Rank + Quantiles
    quickquan(gds, index.gdsn(f, 'mnsc'), onetwo = onetwo, 
                rank = TRUE, new.node = 'mnsrank', ...)
    quickquan(gds, index.gdsn(f, 'unsc'), onetwo = onetwo, 
                rank = TRUE, new.node = 'unsrank', ...)
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
    n.t <- add.gdsn(gds, new.node, storage = 'float64', valdim=c(dim[1],0), 
                    val = NULL, replace = TRUE)
    for(x in 1:dim[2]){
        # This is slow, depending on the number of samples. 
        meth <- mns[, x, name = FALSE]
        unmeth <- uns[, x, name = FALSE]
        beta <- meth/(meth + unmeth + fudge)
        append.gdsn(n.t, beta)
    }
} # }}}
