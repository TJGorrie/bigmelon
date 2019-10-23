# Get quantiles of (pre/un)normalized node.
getquantiles <- function(gds, node, perc, onetwo, rank = FALSE, new.node = NULL){
    x <- index.gdsn(gds, node)
    ranked <- get.attr.gdsn(x)[['ranked']]
    if(!is.null(ranked)){
        out <- list(quantiles = get.attr.gdsn(x)[['quantiles']],
                    inter = get.attr.gdsn(x)[['inter']],
                    onetwo = get.attr.gdsn(x)[['onetwo']])
        return(out)
    } else {
       out <- getquantilesandranks(gds = gds, node = node, onetwo = onetwo, perc = perc, rank.node = NULL)
       return(out)
    }
}

# EstimateCellCounts.gds
estimateCellCounts.gds <- function(
    gds,
    gdPlatform = c("450k", "EPIC", "27k"),
    mn = NULL,
    un = NULL,
    bn = NULL,
    perc = 1,
    compositeCellType = "Blood",
#    processMethod = "auto",
    probeSelect = "auto",
    cellTypes = c("CD8T","CD4T","NK","Bcell","Mono","Gran"),
    referencePlatform = c("IlluminaHumanMethylation450k",
        "IlluminaHumanMethylationEPIC",
        "IlluminaHumanMethylation27k"),
    returnAll = FALSE,
    meanPlot = FALSE,
    verbose = TRUE,
    ...) {
    # My shameless /scopy/s implmentation of minfi::estimateCellCounts
    # For those who do not want use minfi::read.metharray
    referencePlatform <- match.arg(referencePlatform)
    rgPlatform <- gdPlatform <- match.arg(gdPlatform)
    if(rgPlatform == 'EPIC') rgPlatform <- '450k'
    # No method in bigmelon to derive annotation from object, therefore specify.
    # Sanity Checking from minfi...
    if((compositeCellType == "CordBlood") && (!"nRBC" %in% cellTypes)){
        message("[estimateCellCounts] Consider including 'nRBC' in argument 'cellTypes' for cord blood estimation.\n")
    }
    referencePkg <- sprintf("FlowSorted.%s.%s", compositeCellType, rgPlatform)
    if(!require(referencePkg, character.only = TRUE)){
        stop(sprintf("Could not find reference data package for compositeCellType '%s' and referencePlatform '%s' (inferred package name is '%s')",
                     compositeCellType, rgPlatform, referencePkg))
    }
    data(list = referencePkg)
    referenceRGset <- get(referencePkg)
    if(! "CellType" %in% names(colData(referenceRGset)))
        stop(sprintf("the reference sorted dataset (in this case '%s') needs to have a phenoData column called 'CellType'"),
             names(referencePkg))
    if(sum(colnames(gds) %in% colnames(referenceRGset)) > 0)
        stop("the sample/column names in the user set must not be in the reference data ")
    if(!all(cellTypes %in% referenceRGset$CellType))
        stop(sprintf("all elements of argument 'cellTypes' needs to be part of the reference phenoData columns 'CellType' (containg the following elements: '%s')",
                     paste(unique(referenceRGset$cellType), collapse = "', '")))
    if(length(unique(cellTypes)) < 2)
        stop("At least 2 cell types must be provided.")
    # Here is where things get different.
    # First assumption: Data is prenormalized - we will assume dasen is used.
    # Bigmelon has a wrapper that enables the quantification of quantiles for both M and U without having to rerank and re-sort data.
    # Therefore, bigmelon will check if quantiles have been computed already.
    if(is.null(mn) & !'mnsrank'%in% ls.gdsn(gds)) mn <- 'methylated'
    if(is.null(mn) & 'mnsrank' %in% ls.gdsn(gds)) mn <- 'mnsrank'
    if(is.null(un) & !'unsrank'%in% ls.gdsn(gds)) un <- 'unmethylated'
    if(is.null(un) & 'unsrank' %in% ls.gdsn(gds)) un <- 'unsrank'
    if(is.null(bn)) bn <- 'betas'
    referencePd <- colData(referenceRGset)
    referenceMset <- preprocessRaw(referenceRGset)
    if(gdPlatform == 'EPIC'){
    message('No Reference Set of EPIC, converting array to 450k...')
    message('This method is NOT memory efficient!')
    M <- gds[rownames(referenceMset), , node = mn]
    U <- gds[rownames(referenceMset), , node = un]
    rownames(M) <- rownames(U) <- rownames(referenceMset)
    ot <- getProbeType(referenceMset)
    sMI <- wateRmelon:::.normalizeQuantiles2(M[ot=='I',])
    sMII <- wateRmelon:::.normalizeQuantiles2(M[ot=='II',])
    sUI <- wateRmelon:::.normalizeQuantiles2(U[ot=='I',])
    sUII <- wateRmelon:::.normalizeQuantiles2(U[ot=='II',])
    mquan <- list(quantiles = rep(0, nrow(M)),
                    inter = rep(0, nrow(M)),
                    onetwo = ot
                )
    mquan[['quantiles']][ot == 'I'] <- sMI[[1]]
    mquan[['inter']][ot == 'I'] <- sMI[[2]]
    mquan[['quantiles']][ot == 'II'] <- sMII[[1]]
    mquan[['inter']][ot == 'II'] <- sMII[[2]]
    uquan <- list(quantiles = rep(0, nrow(M)),
                    inter = rep(0, nrow(M)),
                    onetwo = ot)
    uquan[['quantiles']][ot == 'I'] <- sUI[[1]]
    uquan[['inter']][ot == 'I'] <- sUI[[2]]
    uquan[['quantiles']][ot == 'II'] <- sUII[[1]]
    uquan[['inter']][ot == 'II'] <- sUII[[2]]
    mquan[['rn']] <- uquan[['rn']] <- rownames(M)

    } else {
    ot <- fot(gds)
    mquan <- getquantiles(gds = gds, node = mn, onetwo = ot, perc = perc)
    # rownames not working(?)
    mquan[['rn']] <- read.gdsn(index.gdsn(gds,read.gdsn(index.gdsn(gds, "paths"))[1]))
    uquan <- getquantiles(gds = gds, node = un, onetwo = ot, perc = perc)
    # rownames not working(?)
    uquan[['rn']] <- read.gdsn(index.gdsn(gds,read.gdsn(index.gdsn(gds, "paths"))[1]))
    # We can skip combining data-sets.
    # Instead preprocessRaw reference to replace data with quantiles.
    }
    nmet <- wateRmelon:::.impose(getMeth(referenceMset), mquan)
    nume <- wateRmelon:::.impose(getUnmeth(referenceMset), uquan)
    rm(referenceRGset)

    # Everything else continues as normal.
    referenceMset <- minfi::MethylSet(Meth=na.omit(nmet), Unmeth=na.omit(nume), colData=referencePd, annotation(referenceMset))

    if(verbose) message("[estimateCellCounts] Picking probes for composition estimation.\n")
    compData <- minfi:::pickCompProbes(referenceMset, cellTypes = cellTypes, compositeCellType = compositeCellType, probeSelect = probeSelect)
    coefs <- compData$coefEsts
    rm(referenceMset)

    if(verbose) message("[estimateCellCounts] Estimating composition.\n")
    coefdat <- gds[rownames(coefs),,node=bn]
    rownames(coefdat) <- rownames(coefs)
    counts <- minfi:::projectCellType(coefdat, coefs)
    # counts <- minfi:::projectCellType(gds[rownames(coefs), , node = bn ], coefs)
    if (meanPlot) {
        smeans <- compData$sampleMeans
        smeans <- smeans[order(names(smeans))]
        sampleMeans <- c(colMeans(gds[rownames(coefs),, node = bn]), smeans)
        sampleColors <- c(rep(1, ncol(coefdat)), 1 + as.numeric(factor(names(smeans))))
        plot(sampleMeans, pch = 21, bg = sampleColors)
        legend("bottomleft", c("blood", levels(factor(names(smeans)))),
               col = 1:7, pch = 15)
    }
    if(returnAll) {
        list(counts = counts, compTable = compData$compTable)
    } else {
        counts
    }
}
