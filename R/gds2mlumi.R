gds2mlumi <- function(gds, i, j){
    history.submitted = as.character(Sys.time())
    x <- gds
    if("NBeads"%in%ls.gdsn(x)){
        aDat <- assayDataNew(
                    betas = x[i, j, node = "betas",
                                name = TRUE, drop = FALSE],
                    pvals = x[i, j, node = "pvals",
                                name = TRUE, drop = FALSE],
                    NBeads = x[i, j, node = "NBeads", 
                                name = TRUE, drop = FALSE],
                    methylated = x[i, j, node = "methylated", 
                                    name = TRUE, drop = FALSE],
                    unmethylated = x[i, j, node = "unmethylated", 
                                        name = TRUE, drop = FALSE])
    } else {
        aDat <- assayDataNew(
                    betas = x[i, j, node = "betas",
                                name = TRUE, drop = FALSE],
                    pvals = x[i, j, node = "pvals",
                                name = TRUE, drop = FALSE],
                    methylated = x[i, j, node = "methylated",
                                    name = TRUE, drop = FALSE],
                    unmethylated = x[i, j, node = "unmethylated",
                                        name = TRUE, drop = FALSE])
    }
    # Creating MethyLumiSet
    x.lumi = new("MethyLumiSet", assayData=aDat)
    pdat <- pData(x)
    rownames(pdat) <- colnames(x)
    pData(x.lumi) <- pdat[j, , drop = FALSE]
    fdat <- fData(x)
    rownames(fdat) <- rownames(x)
    fData(x.lumi) <- fdat[i, , drop = FALSE]
    if(length(grep("QC", ls.gdsn(x), ignore.case = TRUE))>1){
        qcm <- QCmethylated(x)
        qcu <- QCunmethylated(x)
        colnames(qcm) <- colnames(qcu) <- colnames(x)
        rownames(qcm) <- rownames(qcu) <- QCrownames(x)
        qc <- new("MethyLumiQC",
                    assayData = assayDataNew(
                        methylated = qcm[ , j, drop = FALSE],
                        unmethylated = qcu[ , j, drop = FALSE])
                    )
        x.lumi@QC <- qc
    }
    # Forgotton things: 
    #  x.lumi@protocolData <- protocolData(NChannelSet)
    #  x.lumi@annotation <- annotation(NChannelSet)
    #  x.lumi@QC@annotation <- annotation(NChannelSet)
    # fvarLabels(x.lumi) <- possibleLabels[1:ncol(fdat)]
    # fvarMetadata(x.lumi)[,1] <- possibleMetadata[1:ncol(fdat)]
    # pval.detect(x.lumi) <- pval # default value
    history.finished <- as.character(Sys.time())
    history.command <- "Converted to methylumi with gds2mlumi (bigmelon)"
    x.lumi@history <- rbind(getHistory(x), 
                            data.frame( submitted = history.submitted, 
                                        finished = history.finished,
                                        command = history.command))
#    rownames(x.lumi) <- rownames(x)[i]
    return(x.lumi)
}

gds2mset <- function(gds, i, j, anno = NULL){ 
    x <- gds
    if(!is.null(anno)){ 
        if(!anno %in% c("27k", "450k", "epic")){
        stop("anno needs to be: \'27k\', \'450k\', \'epic\'")
        }
    }
    M <- x[i = i, j = j,   "methylated", name = TRUE, drop = FALSE]
    U <- x[i = i, j = j, "unmethylated", name = TRUE, drop = FALSE]
    pd <- pData(x)[j, , drop = FALSE]
    rownames(pd) <- colnames(x)[j]
    #    pd <- annotatedDataFrameFrom(object = as.matrix(pd), byrow = TRUE)
    if(!is.null(anno)){
        if(anno == "27k"){
            anno <- c("IlluminaHumanMethylation27k", "ilmn12.hg19")
        } else if(anno == "450k"){ 
            anno <- c("IlluminaHumanMethylation450k", "ilmn12.hg19")
        } else if(anno == "epic"){ 
            anno <- c("IlluminaHumanMethylationEPIC", "ilmn12.hg19")
        } else if(anno == "unknown"){
            anno <- c("Unknown", "Unknown")
        }
    }
    # Guess Array Type - will not get correct array if performed on subset.
    if(is.null(anno)){
        nr <- nrow(fData(x))
        # Will guess array type based on number of rows, will fail on subsets!
        if(nr > 50000 & nr < 500000){
            anno <- c("IlluminaHumanMethylation450k", "ilmn12.hg19")
        } else if(nr >= 500000){ 
            anno <- c("IlluminaHumanMethylationEPIC", "ilmn12.hg19")
        } else if(nr <=50000){ 
            anno <- c("IlluminaHumanMethylation27k", "ilmn12.hg19")
        }
    }
    names(anno) <- c("array", "annotation")
    out <- MethylSet(Meth = M, Unmeth = U, colData = pd, annotation = anno)
    out@preprocessMethod <- c(
        rg.norm = "Converted from gdsfmt to MethylSet (bigmelon)",
        minfi = as.character(packageVersion("minfi")),
        manifest = NA #packageVersion(getManifest(anno))
        )
    out
}
