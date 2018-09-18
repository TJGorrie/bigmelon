.getEstimate2 <- function(mat, design, coef, B = NULL, permutations = NULL, full = FALSE){
  # mat = node (i.e betas)
  p <- objdesp.gdsn(mat)$dim[2]
  n <- objdesp.gdsn(mat)$dim[1]
  v <- design[, coef]
  A <- design[, -coef, drop = FALSE]
  qa <- qr(A)
  S <- diag(nrow(A)) - tcrossprod(qr.Q(qa)) # ncol * ncol matrix, "small"
  vv <- if(is.null(B)){
    matrix(v, ncol = 1)
  } else {
    if (is.null(permutations)) {
      replicate(B, sample(v))
    } else {
      apply(permutations, 2, function(i) v[i])
    }
  }
  sv <- S %*% vv
  vsv <- diag(crossprod(vv, sv))
  #### Code Changes here!
  if(full){
    df.residual <- p - qa$rank - 1
    if(is.null(B)){
      o <- apply.gdsn(node = mat, margin = 1, as.is = 'list',
      #target.node = list(n.b, n.r),
          FUN = function(x, S, vv, vsv, sv, df.residual){
            sy <- x %*% S
            b <- (x %*% crossprod(S, vv))/vsv
            tcross <- tcrossprod(b, sv)
            sigma <- sum((sy - tcross)^2)/df.residual
            list('B'=b, 'sigma'= sigma)
          }, 
          S = S, vv = vv, vsv = vsv, sv = sv, 
          df.residual = df.residual
      )
    } else {
      o <- apply.gdsn(node=mat, margin=1, as.is = 'list', 
      #target.node = list(n.b, n.r),
         FUN = function(x, S, vv, vsv, sv, B, df.residual){
           tmp <- sy <- x %*% S
           sigma <- b <- (x %*% crossprod(S, vv))/vsv
           for(j in seq_len(B)){
             tmp <- tcrossprod(b[,j], sv[,j])
             sigma[j] <- sum((sy-tmp)^2)
           }
           sigma <- sqrt(sigma/df.residual)
           list('B'= b, 'sigma'=sigma)
         }, S = S, vv = vv, vsv = vsv, sv = sv, 
         df.residual = df.residual, B = B
      )
    }
    coef <- if(is.null(B)) sapply(o, '[[', 'B') else t(sapply(o, '[[', 'B'))
    sigma <- if(is.null(B)) sapply(o, '[[', 'sigma') else t(sapply(o, '[[', 'sigma'))
    out <- list(coef = coef, # n * B big
        sigma = sigma, # n * B big
        stdev.unscaled = sqrt(1/vsv), 
        df.residual = df.residual)
    if(is.null(B)) out$stdev <- as.numeric(out$stdev)
  } else {
    out <- apply.gdsn(node=mat, margin = 1, as.is = 'list',
       FUN = function(x, S, vv, vsv){
          b <- (x %*% crossprod(S, vv))/vsv
        },
        S = S, vv = vv, vsv = vsv
    )
    out <- do.call(rbind, out)
  }
  return(out)
}


bumphunterEngine.gdsn <-function(mat, design, chr = NULL, pos,
   cluster = NULL, coef = 2,
   cutoff = NULL, pickCutoff = FALSE,
   pickCutoffQ = 0.99,
   maxGap = 500,
   nullMethod = c("permutation","bootstrap"),
   smooth = FALSE,
   smoothFunction = locfitByCluster,
   useWeights = FALSE,
   B = ncol(permutations),
   permutations = NULL,
   verbose = TRUE, ...){
      nullMethod  <- match.arg(nullMethod)
      n <- objdesp.gdsn(mat)$dim[1]
      p <- objdesp.gdsn(mat)$dim[2]
      if (is.null(B))  B = 0
      if (!is.matrix(permutations) & !is.null(permutations)) stop(
         "permutations must be NULL or a matrix."
      )
    if (!is.null(permutations)) {
    if (nrow(design) != nrow(permutations)) stop(
     "Number of rows of 'design' must match number of rows of 'permutations'."
    )
    if (B != ncol(permutations)) {
      warning("Ignoring the supplied B. Using 'ncol(permutations)' instead.")
      B = ncol(permutations)
    }
  }
  if (ncol(design) > 2 & B > 0 & nullMethod == "permutation"){
    message("[bumphunterEngine] The use of the permutation test is not recommended with multiple covariates, (ncol(design)>2). Consider changing 'nullMethod' changed to 'bootstrap' instead. See vignette for more information.")
    warning("The use of the permutation test is not recommended with multiple covariates, (ncol(design)>2). Consider changing 'nullMethod' changed to 'bootstrap' instead. See vignette for more information.")
  }
  if (p != nrow(design))
    stop("Number of columns of 'mat' must  match number of rows of 'design'")
  if (B < 0)
    stop("'B' has to be an integer greater or equal to zero")
  if (!(is.null(cutoff) || length(cutoff) %in% 1:2))
    stop("'cutoff' has to be either NULL or a vector of length 1 or 2")
  if (length(cutoff) == 2)
    cutoff <- sort(cutoff)
  if (is.null(cutoff) && !pickCutoff)
    stop("Must pick a cutoff, either using 'cutoff' or 'pickCutoff'")
  if (!is.null(cutoff) && pickCutoff) {
    pickCutoff <- FALSE
    warning("'cutoff' was supplied so ignoring 'pickCutoff=TRUE'")
  }
  if (pickCutoff && (length(pickCutoffQ) != 1 || pickCutoff <  0 || pickCutoffQ > 1))
    stop("Using `pickCutoff = TRUE' requires that 'pickCutoffQ' is a single number between 0 and 1")
  if (pickCutoff && B < 1)
    stop("Must do at least one permution to pick a cutoff")
  ####
  if (!getDoParRegistered())
    registerDoSEQ()
  workers <- getDoParWorkers()
  backend <- getDoParName()
  version <- getDoParVersion()
  ####
  subverbose <- max(as.integer(verbose) - 1L, 0)
  if (verbose) {
    if (workers == 1) {
      mes <- "[bumphunterEngine] Using a single core (backend: %s, version: %s)."
      message(sprintf(mes, backend, version))
    }
    else {
      mes <- "[bumphunterEngine] Parallelizing using %s workers/cores (backend: %s, version: %s)."
      message(sprintf(mes, workers, backend, version))
    }
  }
  if (is.null(chr))
    chr <- rep("Unspecified", length(pos))
  if (is.factor(chr))
    chr <- as.character(chr)
  if (is.null(cluster))
    cluster <- bumphunter::clusterMaker(chr, pos, maxGap = maxGap)
  if (verbose)
    message("[bumphunterEngine] Computing coefficients.")
  if (useWeights & smooth) { 
    tmp <- .getEstimate2(mat = mat, design = design,
                         coef = coef, full = TRUE)
    rawBeta <- tmp$coef
    weights <- 1/tmp$sigma
    rm(tmp)
  } else {
    rawBeta <- .getEstimate2(mat = mat, design = design, coef = coef,
                             full = FALSE)
    weights <- NULL
  }
  if (smooth) {
    if (verbose)
      message("[bumphunterEngine] Smoothing coefficients.")
    beta <- smoother(y = rawBeta, x = pos, cluster = cluster,
                     weights = weights, smoothFunction = smoothFunction,
                     verbose = subverbose, ...)
    Index <- which(beta$smoothed)
    beta <- beta$fitted
  } else {
    beta <- rawBeta
    Index <- seq(along = beta)
  }
  if (B > 0) {
    if (nullMethod == "permutation"){
      if (verbose)
        message("[bumphunterEngine] Performing ", B, " permutations.")
      if (useWeights && smooth) {
        tmp <- .getEstimate2(mat, design, coef, B,
            permutations, full = TRUE
        )
        permRawBeta <- tmp$coef
        weights <- tmp$sigma
        rm(tmp)
      }
      else {
        permRawBeta <- .getEstimate2(mat, design, coef, B, 
           permutations, full = FALSE
        )
        weights <- NULL
      }
      NullBeta<-as.matrix(permRawBeta)	
    }
    if (nullMethod == "bootstrap"){
      message("[bumphunterEngine] Performing ", B, " bootstraps.")
      qr.X <- qr(design)
      ##rescale residuals
      h <- diag(tcrossprod(qr.Q( qr(design))))
      ##create the null model to which we add bootstrap resids
      design0 <- design[,-coef,drop=FALSE]
      qr.X0 <- qr(design0)
      boots <- createfn.gds('bs.gds', allow.duplicate = TRUE)
      res <- add.gdsn(node = boots, name='resids', val = NULL, 
         storage = 'float64', valdim = c(p,0)
      )
      null <- add.gdsn(node = boots, name='null', val = NULL, 
         storage = 'float64', valdim = c(p,0)
      )
      apply.gdsn(node = mat, margin = 1, as.is = 'gdsnode', 
         target.node = list(x=res, y=null), 
         FUN = function(x, s1, s2, n1){
            res <- t(s1 %*% x)/s2
                null <- t(n1 %*% x)
                list(x=res, y=null)
         },
         s1 = t(diag(nrow(design)) - tcrossprod(qr.Q(qr.X))),
         s2 = sqrt(1-h), n1 = tcrossprod(qr.Q(qr.X0)) )

      ##Now do the boostraps
      chunksize <- ceiling(B/workers) 
      bootIndexes<-replicate(B, sample(1:p,replace=TRUE),simplify=TRUE)
      
      # Tyler is here :(
      tmp <- lapply(seq_len(ncol(bootIndexes)),
          FUN = function(x, resids, null, s1,s2,s3,s4,s5,coef, useWeights){
            outList <- apply.gdsn(list(x=resids, y=null), margin=c(2,2), as.is='list',
              FUN = function(X, j, s1, s2, s3, s4, s5, useWeights, coef){
                # create null model
                matstar <- X$y + X$x[j]
                # compute estimate
                nullbetas <- backsolve(s1, crossprod(s2, matstar))[coef]
                if(useWeights) {
                  # compute sigma
                  sigma <- sqrt(sum((s4%*%matstar)^2)/s5)
                  outList <- list(coef = nullbetas, sigma = sigma)
                } else {
                  outList <- nullbetas
                }
                return(outList)
              }, j = bootIndexes[,x], 
              s1 = s1, 
              s2 = s2, 
              s3 = s3, 
              s4 = s4,
              s5 = s5, 
              useWeights = useWeights, 
              coef = coef
            )
            if(useWeights) return(list(
               coef = sapply(outList, '[[', 'coef'), 
               sigma = sapply(outList, '[[', 'sigma'))
            )
            else return(unlist(outList))
          }, 
          resids = index.gdsn(boots, 'resids'), 
          null = index.gdsn(boots, 'null'), 
          s1 = qr.R(qr.X), 
          s2 = qr.Q(qr.X), 
          useWeights = useWeights, 
          coef = coef, 
          s3 = tcrossprod(qr.Q(qr.X)), 
          s4 = t(diag(nrow(design))-tcrossprod(qr.Q(qr.X))), 
          s5 =  (nrow(design) - qr.X$rank))
      if (useWeights && smooth) { # Here...
        bootRawBeta <- do.call(Map, c(cbind, tmp))$coef 
        # or sapply(tmp, '[[' ,'coef')
        weights <- do.call(Map, c(cbind, tmp))$sigma
      } else {
        bootRawBeta <- do.call(cbind, tmp)
        weights <- NULL
      }
      NullBeta<-bootRawBeta
      rm(tmp)
      rm(bootRawBeta)
      closefn.gds(boots)
      unlink(boots[[1]])
    }
    
    if (verbose)
      message("[bumphunterEngine] Computing marginal ",nullMethod," p-values.")
      sumGreaterOrEqual <- rowSums(
         bumphunter:::greaterOrEqual(abs(NullBeta),
         abs(as.vector(rawBeta)))
      )
      pvs <- (sumGreaterOrEqual + 1L)/(B + 1L)
     if (smooth) {
       if (verbose)
         message("[bumphunterEngine] Smoothing ",nullMethod," coefficients.")
         permBeta <- smoother(y = NullBeta, x = pos, cluster = cluster,
         weights = weights, smoothFunction = smoothFunction,
         verbose = subverbose, ...)$fitted
    } else permBeta <- NullBeta
    if (is.null(cutoff))
      cutoff <- quantile(abs(permBeta), pickCutoffQ, na.rm = TRUE)
    if (verbose)
      message(sprintf("[bumphunterEngine] cutoff: %s",
      round(cutoff, 3)))
  }
  if (verbose)
    message("[bumphunterEngine] Finding regions.")
  tab <- regionFinder(x = beta, chr = chr, pos = pos, cluster = cluster,
      cutoff = cutoff, ind = Index, verbose = FALSE)
  if (nrow(tab) == 0) {
    if (verbose)
      message("[bumphunterEngine] No bumps found!")
    return(list(table = NA, coef = rawBeta, fitted = beta,
      pvaluesMarginal = NA))
  } else {
    if (verbose)
      message(sprintf("[bumphunterEngine] Found %s bumps.",
         nrow(tab))
      )
  }
  if (B < 1) {
    return(list(table = tab, coef = rawBeta, fitted = beta,
      pvaluesMarginal = NA))
  }
  if (verbose)
    message("[bumphunterEngine] Computing regions for each ",nullMethod,".")
  chunksize <- ceiling(B/workers)
  subMat <- NULL
  nulltabs <- foreach(subMat = iter(permBeta, by = "col", chunksize = chunksize),
    .combine = "c", .packages = "bumphunter") %dorng% {
      apply(subMat, 2, regionFinder, chr = chr, pos = pos,
            cluster = cluster, cutoff = cutoff, ind = Index,
            verbose = FALSE)
    }
  attributes(nulltabs)[["rng"]] <- NULL
  if (verbose)
    message("[bumphunterEngine] Estimating p-values and FWER.")
  L <- V <- A <- as.list(rep(0, B))
  for (i in 1:B) {
    nulltab <- nulltabs[[i]]
    if (nrow(nulltab) > 0) {
      L[[i]] <- nulltab$L
      V[[i]] <- nulltab$value
      A[[i]] <- nulltab$area
    }
  }
  computation.tots <- function(tab, V, L) {
    Lvalue <- cbind(tab$L, abs(tab$value))
    chunksize <- ceiling(nrow(Lvalue)/workers)
    subL <- NULL
    tots <- foreach(subL = iter(Lvalue, by = "row", chunksize = chunksize),
        .combine = "cbind", .packages = "bumphunter") %dorng%
        {
          apply(subL, 1, function(x) {
            res <- sapply(seq(along = V), function(i) {
              sum(bumphunter:::greaterOrEqual(L[[i]], x[1]) &
                    bumphunter:::greaterOrEqual(abs(V[[i]]),
                                                x[2]))
            })
            c(mean(res > 0), sum(res))
          })
        }
    attributes(tots)[["rng"]] <- NULL
    rate1 <- tots[1, ]
    pvalues1 <- tots[2, ]/sum(sapply(nulltabs, nrow))
    return(list(rate1 = rate1, pvalues1 = pvalues1))
  }
  ##    ptime1 <- proc.time()
  comp <- computation.tots(tab = tab, V = V, L = L)
  rate1 <- comp$rate1
  pvalues1 <- comp$pvalues1
  ##    ptime2 <- proc.time()
  computation.tots2 <- function(tab, A) {
    Avalue <- matrix(tab$area, ncol = 1)
    chunksize <- ceiling(nrow(Avalue)/workers)
    subA <- NULL
    tots2 <- t(foreach(subA = iter(Avalue, by = "row", chunksize = chunksize),
           .combine = "cbind", .packages = "bumphunter") %dorng%
           {
             sapply(subA, function(x) {
               return(sapply(seq(along = A), function(i) {
                 sum(bumphunter:::greaterOrEqual(A[[i]], x[1]))
               }))
             })
           })
    if (is.vector(tots2)) {
      tots2 <- matrix(tots2, nrow = 1)
    }
    rate2 <- rowMeans(tots2 > 0)
    pvalues2 <- rowSums(tots2)/sum(sapply(nulltabs, nrow))
    return(list(rate2 = rate2, pvalues2 = pvalues2))
  }
  ##ptime1 <- proc.time()
  comp <- computation.tots2(tab = tab, A = A)
  rate2 <- comp$rate2
  pvalues2 <- comp$pvalues2
  ##ptime2 <- proc.time()
  tab$p.value <- pvalues1
  tab$fwer <- rate1
  tab$p.valueArea <- pvalues2
  tab$fwerArea <- rate2
  tab <- tab[order(tab$fwer, -tab$area), ]
  algorithm <- list(version = packageDescription("bumphunter")$Version,
      coef = coef, cutoff = cutoff, pickCutoff = pickCutoff,
      pickCutoffQ = pickCutoffQ, smooth = smooth, maxGap = maxGap, nullMethod=nullMethod,
      B = B, permutations = permutations, useWeights = useWeights,
      smoothFunction = deparse(substitute(smoothFunction)))
  ret <- list(table = tab, coef = rawBeta, fitted = beta, pvaluesMarginal = pvs,
      null = list(value = V, length = L), algorithm = algorithm)
  class(ret) <- "bumps"
  return(ret)
}
