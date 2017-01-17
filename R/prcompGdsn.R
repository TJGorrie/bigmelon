# gds.class "method" for prcomp
# Notable changes: 
# [1] Scaling and centering is NOT performed by default as
#     prcomp is usually performed on normalized data.
#     Notable Time save as samples are expected to be normalized already.
# [2] Number of learned principal components is optional,
#     5 are default, however this is assuming that majority of variance can be
#     found in the first few PCs hence why the output is slimmed down.
#     No Notable time save (May affect timesave when retx=T) 
#     Due to how the principal components are learned 
#     it is not a memory efficient operation
#     Therefore the forthcoming approach is chosen:
# [3] Perform pca on a subset of data (optionally), by default we opt for 1% 
#     As it is felt that considering the wide berth of probes 1% is enough to
#     learn the rough indication of the principal components. 
#     Notable Time save (obviously), limits memory usage as svd is performed
#     on a small subset. 
# 
prcomp.gds.class <- function(x, node.name, center = FALSE, scale. = FALSE, 
                            rank. = NULL, retx=FALSE, tol = NULL,
                            perc = 0.01, npcs = NULL, 
                            parallel = NULL, method = c('quick', 'sorted'),
                            verbose = FALSE, ...){
    method <- match.arg(method)
    if(method=="quick") parallel <- NULL
    # Parallel only amenable to calculating variance for sorted method.
    if(!is.null(parallel)){
        if(!inherits(parallel, 'cluster')){
            if(requireNamespace('parallel') & as.integer(parallel) > 1){
                if(verbose) message('Making Cluster of ', parallel, ' cores.')
                parallel <- parallel::makeCluster(parallel)
            } else {
                message('Cannot make cluster, continuing on single core.')
                parallel <- NULL
            }
        }
    }

    # x = gfile
    # node.name = node to pca
    # method = quick or sorted
    # parallel = if Multicore to be used for applicable steps!
    x0 <- x # Original gfile 
    x <- inode <- index.gdsn(x, node.name)
    dim <- objdesp.gdsn(x)$dim

    f <- createfn.gds("temp.gds", allow.duplicate=TRUE)
    n.t <- add.gdsn(f, val = NULL, valdim = c(dim[1], 0),
                    replace = TRUE, name = "scale", storage = "float64")
    cen <- NULL
    sc <- NULL
    if(center|scale.){ # If normalized Center | Scale = FALSE
        if(verbose) message('Scaling data... ')
#        if(!is.null(parallel)){
#           censc <- clusterApply.gdsn(cl = parallel,
#                              gds.fn = x0$filename,
#                              node.name = node.name,
#                              margin = 2,
#                              as.is="double",
#                              FUN=function(val, center, scale.,n.t){
#                              y <- scale(val, center=center, scale=scale.)
#                              append.gdsn(n.t, y)
#                              if(!is.null(attr(y, "scaled:center"))){
#                              out <- attr(y, "scaled:center")
#                              }
#                              if(!is.null(attr(y, "scaled:scale"))){
#                              names(out) <- attr(y, "scaled:scale")
#                              }
#                              }, center=center, scale.=scale., n.t = n.t)
#        } else {
#           censc <- apply.gdsn(x, margin = 2, as.is='double', FUN=function(val, center, scale., n.t){
#                               y <- scale(val, center=center, scale=scale.)
#                              append.gdsn(n.t, y)
#                              if(!is.null(attr(y, "scaled:center"))){
#                              out <- attr(y, "scaled:center")
#                              }
#                              if(!is.null(attr(y, "scaled:scale"))){
#                              names(out) <- attr(y, "scaled:scale")
#                              }
#                              }, center=center, scale.=scale., n.t = n.t)
#        }
         # apply.gdsn could work but getting cen and sc is difficult.
       for(i in 1:dim[2]){
            val <- readex.gdsn(x, sel = list(NULL, i))
            y <- scale(val, center = center, scale = scale.)
            if(!is.null(attr(y, "scaled:center"))){
                cen[i] <- attr(y, "scaled:center")
            }
            if(!is.null(attr(y, "scaled:scale"))){
                sc[i] <- attr(y, "scaled:scale")
            }
            append.gdsn(n.t, y)
        }
        x <- n.t 
    }
    n <- dim[1]
    p <- dim[2]
    k <- if(!is.null(rank.)) {
        stopifnot(length(rank.) == 1, is.finite(rank.),
                    as.integer(rank.) > 0)
        min(as.integer(rank.), n, p, npcs)
            ## Note that La.svd() *still* needs
            # a (n x p) and a (p x p) auxiliary
        } else {
        min(npcs, n, p)
        }
    # Current method would be to run a % of probes and learn a fewer pcs...
    
    sel <- rep(FALSE, times = n)
    if(method == "sorted"){
        if(verbose) message('Identifying largest variance in rows...')
        if(!is.null(parallel)){
            if(center|scale.){
              closefn.gds(f)
              f <- openfn.gds('temp.gds', allow.duplicate = TRUE)
              x <- n.t <- index.gdsn(f, 'scale') 
            }
            variances <- clusterApply.gdsn(cl = parallel,
                                           gds.fn = if(center | scale.) 'temp.gds' else x0$filename,
                                           node.name = if(center | scale.) 'scale' else node.name,
                                           margin = 1,
                                           as.is = 'double',
                                           FUN = function(val){
                                           isna <- sum(is.na(val))>0
                                           fn <- quantile(val, na.rm = TRUE)
                                           if(isna) out <- NA else out <- fn[4] - fn[2]
                                           out
                                           } )
        } else { 
            variances <- apply.gdsn(node = x, margin = 1, as.is = 'double', FUN = function(val){
                                           isna <- sum(is.na(val)) > 0
                                           fn <- quantile(val, na.rm = TRUE)
                                           if(isna) out <- NA else out <- fn[4] - fn[2]
                                           out
                                           } ) }
    sel[order(variances, decreasing = TRUE)[1:round(n*perc)]] <- TRUE
    } else if(method=="quick"){
        if(verbose) message('Stochastically selecting rows ... ')
        sel[sample(1:n, round(n*perc), replace = FALSE)] <- TRUE
    }
    # Not Memory efficient though!
    mat <- na.omit(readex.gdsn(node = x, sel = list(sel, NULL)))
    if(verbose) message('Learning first ', k, ' Principal components.')
    s <- svd(mat, nu = 0, nv = k)
    j <- seq_len(k)
    #    s$d <- s$d / sqrt(max(1, n - 1))
    # Incase NA's are removed, we need to account for the change selection.
    s$d <- s$d / sqrt(max(1, round(nrow(mat)*perc) - 1))
    if(!is.null(tol)){
        ## we get rank at least one even for a 0 matrix.
        rank <- sum(s$d > (s$d[1L]*tol))
        if(rank < k){
            j <- seq_len(k <- rank)
            s$v <- s$v[,j , drop = FALSE]
        }
    }
    #    dimnames(s$v) <- list(colnames(x), paste0("PC", j)) # Broken
    r <- list(sdev = s$d, rotation = s$v,
                center = if(is.null(cen)) FALSE else cen,
                scale = if(is.null(sc)) FALSE else sc)
    # There isn't a clever way to perform matrix multiplications yet.
    if (retx) {
    # Original Statement:   r$x <- x %*% s$v
    # Slow (Memory Intensive) Method:
    r$x <- read.gdsn(x) %*% s$v
    }
    closefn.gds(f)
    unlink("temp.gds")
    if(!is.null(parallel)){
        if(verbose) message('Stopping Cluster')
        parallel::stopCluster(parallel)
    }
    class(r) <- "prcomp"
    r
}
