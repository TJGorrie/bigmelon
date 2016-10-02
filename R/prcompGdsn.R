# prcomp.gdsn
# gdsn.class "method" for prcomp
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
prcomp.gdsn.class <- function(x, center = FALSE, scale. = FALSE, tol = NULL,
                                rank. = NULL, retx=FALSE,
                                perc = 0.01, npcs = NULL, ...){
    # x = node 
    dim <- objdesp.gdsn(x)$dim
    f <- createfn.gds("temp.gds", allow.duplicate=TRUE)
    n.t <- add.gdsn(f, val = NULL, valdim = c(dim[1], 0),
                    replace = TRUE, name = "scale", storage = "float64")
    x0 <- x
    cen <- NULL
    sc <- NULL
    if(center|scale.){ # If normalized Center | Scale = FALSE
        # Would apply.gdsn suffice here?
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
    nrow <- rep(FALSE, times = n)
    nrow[sample(1:n, round(length(1:n)*perc), replace = FALSE)] <- TRUE
    # Not Memory efficient though! 
    s <- svd(na.omit(readex.gdsn(node = x,
            sel = list(nrow, NULL))), nu = 0, nv = k)
    j <- seq_len(k)
    #    s$d <- s$d / sqrt(max(1, n - 1))
    s$d <- s$d / sqrt(max(1, round(length(1:n)*perc) - 1))
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
    unlink("temp2.gds")
    class(r) <- "prcomp"
    r
}
