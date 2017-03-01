# Get quantiles of (pre/un)normalized node.
getq <- function(gds, node, perc, onetwo, rank = FALSE, new.node = NULL){
    x <- index.gdsn(gds, node)
    ranked <- get.attr.gdsn(x)[['ranked']]
    if(is.null(ranked)) ranked <- FALSE
    if(ranked){
        out <- list(quantiles = get.attr.gdsn(x)[['quantiles']],
                    inter = get.attr.gdsn(x)[['inter']],
                    ot = get.attr.gdsn(x)[['onetwo']])
        return(out)
    } else {
       quickquan(gds = gds, node = node, onetwo = onetwo, perc = perc, rank=rank, new.node = new.node)
    }
}

impose <- function(matrix, quan){
    ot <- quan[['onetwo']]
    quantiles <- quan[['quantiles']]
    names(ot) <- quan[['rn']]
    inter <- quan[['inter']]
    blank <- matrix(NA, length(ot), ncol(matrix))
    rownames(blank) <- names(ot)
    blank[rownames(matrix), ] <- matrix
    for(z in 1:ncol(matrix)){
        isna <- is.na(blank[,z])
        r <- rep(0, length(ot))
        r[ot=='I'] <- rank(blank[ot=='I', z])
        r[ot=='II'] <- rank(blank[ot=='II', z])
        blank[ot == 'I' & (!isna), z] <-  approx(inter[ot == 'I'],
            quantiles[ot == 'I'],
            (r[ot=='I'&(!isna)] - 1) /(sum(ot == 'I')-1), ties = "ordered")$y

        blank[ot == 'II' & (!isna),z] <- approx(inter[ot == 'II'],
            quantiles[ot == 'II'],
            (r[ot=='II'&(!isna)] - 1)/(sum(ot == 'II')-1), ties = "ordered")$y
    }
    b2 <- na.omit(blank)
    b2<- b2[rownames(matrix),]
    return(b2)
}

# EstimateCellCounts.gdsn
estimateCellCounts.gdsn <- function(
    gds,
    gdPlatform = c("450k", "EPIC", "27k"),
    mn = NULL,
    un = NULL,
    perc = 0.25,
    compositeCellType = "Blood",
    processMethod = "auto",
    probeSelect = "auto",
    cellTypes = c("CD8T","CD4T","NK","Bcell","Mono","Gran"),
    referencePlatform = c("IlluminaHumanMethylation450k",
        "IlluminaHumanMethylationEPIC",
        "IlluminaHumanMethylation27k"),
    returnAll = FALSE,
    meanPlot = FALSE,
    verbose=TRUE,
    ...) {
    # My shameless /scopy/s implmentation of minfi::estimateCellCounts
    # For those who do not want use minfi::read.metharray
    referencePlatform <- match.arg(referencePlatform)
    rgPlatform <- match.arg(gdPlatform) # No method in bigmelon to derive annotation from object, therefore specify.
    if(!sub("IlluminaHumanMethylation", "", referencePlatform) == rgPlatform){
        stop("Reference and gdPlatform must match.")
    }
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
    # MORE SANITY
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
    if(is.null(mn)&!'mnsrank'%in%ls.gdsn(gds)) mn <- 'methylated'
    if(is.null(mn)&'mnsrank'%in%ls.gdsn(gds)) mn <- 'mnsrank'
    if(is.null(un)&!'unsrank'%in%ls.gdsn(gds)) un <- 'unmethylated'
    if(is.null(un)&'unsrank'%in%ls.gdsn(gds)) un <- 'unsrank'
    ot <- fot(gds)
    mquan <- getq(gds = gds, node = mn, onetwo = ot, perc = perc)
    # rownames not working(?)
    mquan[['rn']] <- read.gdsn(index.gdsn(gds,read.gdsn(index.gdsn(gds, "paths"))[1]))
    uquan <- getq(gds = gds, node = un, onetwo = ot, perc = perc)
    # rownames not working(?)
    uquan[['rn']] <- read.gdsn(index.gdsn(gds,read.gdsn(index.gdsn(gds, "paths"))[1]))
    # We can skip combining data-sets.
    # Instead preprocessRaw reference to replace data with quantiles.
    referencePd <- colData(referenceRGset)
    referenceMset <- preprocessRaw(referenceRGset)
    nmet <- impose(getMeth(referenceMset), mquan)
    nume <- impose(getUnmeth(referenceMset), uquan)
    rm(referenceRGset)

    referenceMset <- minfi::MethylSet(Meth=na.omit(nmet), Unmeth=na.omit(nume), colData=referencePd, annotation(referenceMset))

    if(verbose) message("[estimateCellCounts] Picking probes for composition estimation.\n")
    compData <- minfi:::pickCompProbes(referenceMset, cellTypes = cellTypes, compositeCellType = compositeCellType, probeSelect = probeSelect)
    coefs <- compData$coefEsts
    rm(referenceMset)

    if(verbose) message("[estimateCellCounts] Estimating composition.\n")
    counts <- minfi:::projectCellType(betas(gds)[rownames(coefs), ], coefs)

    if (meanPlot) {
        smeans <- compData$sampleMeans
        smeans <- smeans[order(names(smeans))]
        sampleMeans <- c(colMeans(betas(gds)[rownames(coefs), ]), smeans)

        sampleColors <- c(rep(1, ncol(mSet)), 1 + as.numeric(factor(names(smeans))))
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
