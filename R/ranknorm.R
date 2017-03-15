.sortvector <- function(vec){
    dim <- length(vec)
    i <- (0:(dim-1))/(dim-1)
    S <- rep(NA, dim)
    si <- sort(vec, method = 'quick', index.return = TRUE)
    nobsj <- length(si$x)
    if(nobsj < dim[1]){
        S <- approx((0:(nobsj-1))/(nobsj-1),
                    si$x, i, ties = 'ordered')$y
    } else {
        S <- si$x
    }
    return(S)
}
# Quite slow, think about breaking it down.
getquantilesandranks <- function(gds, node, onetwo, rank.node = NULL, perc = 1){
    if(!is.null(rank.node)&!is.character(rank.node)) stop('rank.node needs to be a character string is supplied.')
    if(!(perc <= 1 & perc > 0)) stop('perc must be greater than 0 and less than or equal to 1.')
    datnod <- node
    # If node is single element vector - index it.
    if(length(datnod) == 1) datnod <- index.gdsn(gds, as.character(datnod))
    dim <- objdesp.gdsn(datnod)$dim
    if(!length(onetwo) == dim[1]) stop('Length of onetwo and nrow of gdsn node do not match!')
    nc <- rep(TRUE, dim[2])
    if(perc < 1) nc <- which(nc)%in%sample(which(nc), round(length(nc)*perc), replace = FALSE)
    f <- createfn.gds("quickquan.gds", allow.duplicate = TRUE)
    quants <- rep(NA, dim[1])
    # Split + Sort Array by probe type
    nodeI <- add.gdsn(node = f, name = 'nodeI', storage = 'float64',
                      valdim = c( sum(onetwo=='I'), 0), val = NULL, replace = TRUE)
    nodeII <- add.gdsn(node = f, name = 'nodeII', storage = 'float64',
                        valdim = c(sum(onetwo == 'II'), 0), val = NULL, replace = TRUE)
    apply.gdsn(node = datnod,
                target.node = nodeI,
                as.is = 'gdsnode',
                margin = 2,
                selection = list(onetwo == 'I', nc),
                var.index = 'relative',
                FUN = function(index, x){ .sortvector(x) } )

    # This can be parallelized!
    quants[onetwo == 'I'] <- apply.gdsn(node = nodeI,
                                        as.is = 'double',
                                        var.index = 'none',
                                        margin = 1,
                                        FUN = mean)
    apply.gdsn(node = datnod,
                target.node = nodeII,
                as.is = 'gdsnode',
                margin = 2,
                selection = list(onetwo == 'II', nc),
                var.index = 'relative',
                FUN = function(index, x){ .sortvector(x) } )

    # This can be parallelized!!!
    quants[onetwo == 'II'] <- apply.gdsn(node = nodeII,
                                        as.is = 'double',
                                        var.index = 'none',
                                        margin = 1,
                                        FUN = mean)
    if(!is.null(rank.node)){
        rn <- add.gdsn(node = gds, name = rank.node, storage = 'float64',
                        valdim = c(length(onetwo), 0), val = NULL,
                        replace = TRUE)
        rnna <- add.gdsn(node = gds, name = paste0('isna', rank.node),
                        storage = 'int8', valdim = c(length(onetwo), 0),
                        val = NULL, replace = TRUE, visible = FALSE)
        apply.gdsn(node = datnod,
                    target.node = list(rn, rnna),
                    as.is = 'gdsnode',
                    margin = 2,
                    var.index = 'none',
                    FUN = function(x, onetwo){
                        ranks <- rep(NA, length(x))
                        ranks[onetwo=='I'] <- rank(x[onetwo=='I'])
                        ranks[onetwo=='II']<- rank(x[onetwo=='II'])
                        out <- list(ranks, as.numeric(is.na(x)))
                        return(out)
                    }, onetwo = onetwo
                    )
    }
    inter <- rep(0, dim[1])
    inter[onetwo == 'I'] <- (0:(sum(onetwo=='I')-1))/(sum(onetwo=='I')-1)
    inter[onetwo == 'II']<- (0:(sum(onetwo=='II')-1))/(sum(onetwo=='II')-1)
    if(!is.null(rank.node)){
        put.attr.gdsn(rn, 'ranked', val = TRUE)
        put.attr.gdsn(rn, 'is.na', val = paste0('isna', rank.node))
        put.attr.gdsn(rn, 'inter', val = inter)
        put.attr.gdsn(rn, 'quantiles', val = quants)
        put.attr.gdsn(rn, 'onetwo', val = onetwo)
    return(0)
    }

    closefn.gds(f)
    unlink('quickquan.gds', force = TRUE)
    output <- list(quantiles = quants, inter = inter, onetwo = onetwo)
    return(output)
}

dasenrank <- function(gds, mns, uns, onetwo, roco, calcbeta = NULL, perc = 1){# {{{
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
    getquantilesandranks(gds, index.gdsn(f, 'mnsc'), onetwo = onetwo,
                rank.node = 'mnsrank', perc = perc, ...)
    getquantilesandranks(gds, index.gdsn(f, 'unsc'), onetwo = onetwo,
                rank.node = 'unsrank', perc = perc, ...)
    # COMPLETE Return nothing
    if(is.null(calcbeta)){
        message('Run \'computebeta.gds(gds, calcbeta, \'mnsrank\', \'unsrank\', fudge = 100)\' to calculate betas!')
    } else {
        message('Calculating Betas...')
        computebeta.gds(gds, calcbeta, 'mnsrank', 'unsrank', fudge = 100)
    }
    closefn.gds(f)
    unlink('temp.gds', force = TRUE)
} # }}}

computebeta.gds <- function(gds, new.node, mns, uns, fudge = 100){ # {{{
    if(length(mns) == 1) mns <- index.gdsn(gds, mns)
    if(length(uns) == 1) uns <- index.gdsn(gds, uns)
    dim <- objdesp.gdsn(mns)$dim
    n.t <- add.gdsn(gds, new.node, storage = 'float64',
                valdim=c(dim[1],0), val = NULL, replace = TRUE)
    for(x in 1:dim[2]){
        # This may be slow, depending on the number of samples - but memory efficient.
        meth <- mns[, x, name = FALSE]
        unmeth <- uns[, x, name = FALSE]
        beta <- meth/(meth + unmeth + fudge)
        append.gdsn(n.t, beta)
    }
} # }}}
