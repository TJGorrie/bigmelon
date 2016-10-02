dfsfit.gdsn <- function(gds,
                        targetnode,
                        newnode,
                        roco,
                        onetwo
                        ){ # {{{
    # Converting supplied 'nodes' into environment
    if(length(onetwo) == 1) onetwo <- read.gdsn(index.gdsn(gds, onetwo))
    if(class(onetwo) == 'gdsn.class') onetwo <- read.gdsn(onetwo)
    if(length(roco) == 1) roco <- read.gdsn(index.gdsn(gds, roco))
    if(class(roco) == 'gdsn.class')   roco <- read.gdsn(roco)
    datnod <- targetnode
    dim <- objdesp.gdsn(datnod)$dim
    # Replicating: dfs2.R 'apply(mn, 2, dfs2, onetwo)'
    mdf <- apply.gdsn(datnod, margin = 2, FUN = function(val, onetwo){
        one <- density(val[onetwo == 'I'], na.rm = TRUE, n = 2^15,
                        from = 0, to = 5000)
        two <- density(val[onetwo == 'II'], na.rm = TRUE, n = 2^15,
                        from = 0, to = 5000)
        one$x[which.max(one$y)] - two$x[which.max(two$y)] 
        }, as.is = "double", onetwo = onetwo)
    # Replicating: dfsfit.R 
    if(!is.null(roco)){
        scol <- as.numeric(substr(roco,6,6))
        srow <- as.numeric(substr(roco,3,3))
        fit  <- try(  lm(mdf ~ srow + scol ), silent=TRUE) 
        if(!inherits(fit, "try-error")){ 
            mdf <- fit$fitted.values
        } else { 
            message ('Sentrix position model failed, skipping') 
        }
    }
    # Creating newnode:
    n.t <- add.gdsn(gds, newnode, storage = "float64", 
                    valdim = c(dim[1],0), val = NULL, replace=TRUE) 
    for(x in 1:dim[2]){
        val <- readex.gdsn(datnod, sel = list(NULL, x)) 
        val[onetwo=='I'] <- val[onetwo=='I'] - rep(mdf[x], sum(onetwo=='I'))
        # Commit to New Node.
        append.gdsn(n.t, val)
    }
} # }}}
