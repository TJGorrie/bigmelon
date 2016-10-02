qn.gdsn <- function(gds,
                    target,
                    newnode
                    ){ # {{{ 
    if(!length(newnode) == 1) stop("newnode does not have a length of 1!")
    datnod <- target
    dim <- objdesp.gdsn(datnod)$dim
    # Initialising rolling sum
    roll <- rep(0, dim[1])
    nobs <- rep(dim[1], dim[2])
    i <- (0:(dim[1]-1))/(dim[1]-1)
    # Sorting and 'growing' rolling sum.
    for(x in 1:dim[2]){
        val <- readex.gdsn(datnod, sel = list(NULL, x))
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
    }  
    # Calculating 'rowMeans'
    rm <- roll/dim[2]
    # Creating new gdsnode with 'newnode' name
    n.t <- add.gdsn(gds, newnode, storage = "float64", 
                    valdim = c(dim[1], 0), val = NULL, replace = TRUE)
    # Reranking + Replacing Values
    for(x in 1:dim[2]){
        val <- readex.gdsn(datnod, sel = list(NULL, x)) 
        r <- rank(val)
        # If NAs exist - Preserving NAs
        if(nobs[x] < dim[1]) {
            isna <- is.na(val)
            val[!isna] <- approx(i, rm, (r[!isna]-1)/(nobs[x]-1),
                                ties = "ordered")$y
        } else {
            val <- approx(i, rm, (r-1)/(dim[1]-1), ties = "ordered")$y
        }
        # Writing new values to node
        append.gdsn(n.t, val)
    }
} # }}}

