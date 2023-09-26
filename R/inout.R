# newgds -- generate a gdsfile stub {{{
newgds <- function(file){
    message('creating gdsfile...' )
    bmln <- createfn.gds(file)
    n <- add.gdsn(bmln, name = "description")
    put.attr.gdsn(n, "FileFormat", "DNA_Methylation")
    bmln
} # AS - OK }}}

# findgdsobj -- function to check if file is already linked to by an existing gds object {{{
# in the current workspace (global env)
findgdsobj <- function(gds){
    fhandle <- gds
    for(obj in ls(envir = .GlobalEnv)){
        if(inherits(get(obj), 'gds.class') &&
            get(obj)$filename == normalizePath(gds)){
            temp <- deparse(obj)
            temp <- gsub("[[:blank:]\"]", "", temp)
            warning(paste(get(obj)$filename,
                    "is already open, appending to the linked gds object"))
            return(fhandle <- get(obj))
        }
    }
    fhandle
} # AS - OK }}}

# handle -- take a gds.class object, the name of a gds.class object or a gds filename {{{
# and return a writable gds.class object
handle <- function(gds){
    newfile <- as.logical(0)
    # gds can be a filepath or gds.class object or an object name
    # gds.class object given -> use as handle
    if(inherits(gds, 'gds.class' ) && file_test('-f', gds$filename)){
        handle <- gds
    } else if(inherits(gds, 'character')){ # Filepath given
        #file saved in global env
        if(exists(gds) && inherits(g <- get(gds), 'gds.class')){
            g <- get(gds)
            handle <- openfn.gds(g$filename, readonly=FALSE, allow.duplicate=TRUE, allow.fork=TRUE)
        } else if(file_test('-f', gds)){ # On Disk
            #file exists, try to open
            # if this fails, check if it is already linked by a gds object
            handle <- tryCatch({
                openfn.gds(gds, readonly=FALSE, allow.duplicate=TRUE, allow.fork=TRUE)
                }, error = function(e){
                    return(findgdsobj(gds))
                }, finally = {
                })
        } else if(!(file_test('-f', gds))){ # Need to create file
            #handle <- newgds(gds)
            handle <- gds
            newfile <- as.logical(1)
        }
    #gds is not a valid gds.class object or filepath
    } else {
        stop(paste(substitute(gds),
            "is not a valid gds.class object or filepath"))
    }
    list(handle,newfile)
} #}}}

# app2gds -- append one or more arrays of data to gdsfile {{{
# Not functioning for minfi objects!
app2gds <- function(m, bmln){
    history.submitted <- as.character(Sys.time())
    # check that m is a methylumi set
    if(!(inherits(m, 'MethyLumiSet'))){
        stop(paste(deparse(substitute(m)), "is not a MethyLumiSet object"))
    }

    # call handle to deal with file checks
    rehandle <- handle(bmln)
    bmln <- rehandle[[1]]
    newfile <- rehandle[[2]]
    if(!(newfile)){
        if(! length(rownames(bmln)) == dim(m)[1]){
            stop(paste(deparse(substitute(m)), "has a different length to bmln"))
        }
        # More rigid checks of rownames need to be down prior to importing
        message(paste('appending to', bmln$filename))
        if(exists("betas", assayData(m))){
            message('betas... ' )
            append.gdsn(index.gdsn(bmln, "betas"), val = betas(m))
        }
        if(exists("pvals", assayData(m))){
            message('pvals... ' )
            append.gdsn(index.gdsn(bmln, "pvals"), val = pvals(m))
        }
        if(exists("methylated", assayData(m))){
            message('methylated... ' )
            append.gdsn(index.gdsn(bmln, "methylated"), val = methylated(m))
        }
        if(exists("unmethylated", assayData(m))){
            message('unmethylated... ' )
            append.gdsn(index.gdsn(bmln, "unmethylated"),
                        val = unmethylated(m))
        }
        if(exists("NBeads", assayData(m))){
            message('NBeads... ' )
            append.gdsn(index.gdsn(bmln, "NBeads"),
                        val = assayDataElement(m, "NBeads"))
        }
        # fData should be the same for the appended data, no need to add
        message('fData... ' )
        # deal with blank col names in fData
        nb_col <- which(colnames(fData(m)) == "")
        colnames(fData(m))[nb_col] <- nb_col
        message('pData... ' )
        mpdat <- data.frame(lapply(pData(m), as.character),
                            stringsAsFactors = FALSE)
        for(i in colnames(mpdat)){
            # iterate through colnames in m
            pData_coln <- i
            # select corresponding child node of pData in the gds file
            pData_index_str <- paste("pData/", i, sep = "")
            pData_child_n <- index.gdsn(bmln, pData_index_str, silent = TRUE)
            # append
            append.gdsn(pData_child_n, val = mpdat[ , pData_coln])
        }
        if(!is.null(m@QC)){ # If QC is null, not appended.
            message('qcdata... ' )
            if(class(index.gdsn(bmln,
                                'QCmethylated', silent = TRUE)) != 'NULL'){
                append.gdsn(index.gdsn(bmln, "QCmethylated"),
                            val = methylated(QCdata(m)))
            }
            if(class(index.gdsn(bmln,
                                'QCunmethylated', silent = TRUE)) != 'NULL'){
                append.gdsn(index.gdsn(bmln, "QCunmethylated"),
                            val = unmethylated(QCdata(m)))
            }
        }
        # add append operation to history
        numsamp <- nrow(pData(m))
        history.finished <- as.character(Sys.time())
        history.command <- as.character(paste(numsamp,
                                            "arrays appended (bigmelon)"))
        h <- data.frame(submitted = history.submitted,
                        finished = history.finished,
                        command = history.command,
                        stringsAsFactors = FALSE)
        for(j in colnames(h)){
            h_coln <- j
            h_index_str <- paste("history/", j, sep = "")
            h_child_n <- index.gdsn(bmln, h_index_str, silent = TRUE)
            append.gdsn(h_child_n, val = h[ ,h_coln])
        }
        # include history from arrays that have been added
        mhist <- data.frame(lapply(m@history, as.character),
                            stringsAsFactors = FALSE)
    for (k in colnames(mhist)){
        hist_coln <- k
        hist_index_str <- paste("history/",k,sep="")
        hist_child_n <- index.gdsn(bmln,hist_index_str,silent=TRUE)
            append.gdsn(hist_child_n, val=mhist[,hist_coln])
    }
    } else {
    # call es2gds to create new file
    bmln <- es2gds(m,bmln)
    }
    # Close and open in write mode?
    closefn.gds(bmln)
    bmln <- openfn.gds(bmln[[1]], readonly=FALSE, allow.duplicate=TRUE, allow.fork=TRUE)

    return(bmln)
} #}}}

# iadd2 -- add data from all idat files that are stored within a single directory to a gds file {{{
iadd2 <- function(path, gds, chunksize = NULL, force=TRUE,...){
    rown <- TRUE
    if(force){
        thefile <- try(openfn.gds(gds, allow.duplicate=TRUE), silent=T)
        if(!inherits(thefile, 'try-error')){
            # TODO this needs fixing see iadd function
            rown <- read.gdsn(index.gdsn(thefile, 'fData/Probe_ID'))
            closefn.gds(thefile)
        }
    }
    gdsfile <- gds
    barcodes <- bfp(path)
    if(is.null(chunksize) & length(barcodes) > 500){
    chunksize <- 500
    message('More than 500 barcodes identified. Switching to batch mode!')
    }
    if(!is.null(chunksize)){
        if(!is.numeric(chunksize)) stop("chunksize must be numeric!")
        chunks <- seq(1, length(barcodes), chunksize)
        for(i in chunks){
            sets <- barcodes[i:(i+(chunksize-1))]
            # TODO this needs fixing see iadd function
            ml <- methylumIDATepic(barcodes = sets[!is.na(sets)],
                                    n = TRUE, oob = FALSE, idatPath = path, force=force)[rown,]
            gdsfile <- app2gds(ml, gdsfile)
        }
    } else {
        ml <- methylumIDATepic(barcodes, n = TRUE,
                                oob = FALSE, idatPath = path, force=force)[rown,]
        gdsfile <- app2gds(ml, gdsfile)
    }
    gdsfile
} #}}}

# iadd -- add data from multiple, specified, idat files providing to a specified gds file. {{{
iadd <- function (bar, gds, n = TRUE, force=TRUE, target_cpgs = NULL, ...){
    rown <- TRUE
    if(is.null(target_cpgs)){
        if(force){
            thefile <- try(openfn.gds(gds, allow.duplicate=TRUE), silent=T)
            if(!inherits(thefile, 'try-error')){
                rown <- read.gdsn(index.gdsn(thefile, 'fData/Probe_ID'))
                closefn.gds(thefile)
            }
        }
    } else {
        rown <- target_cpgs
    }
    ifile <- basename(bar)

    ### this assumes only 2 underscores, fails when names have _ in prefix, eg from GEO
    ### need to count from end instead or just chop off [Red|Grn].idat
    #pieces <- strsplit(ifile, "[_.]")
    #slide <- sapply(pieces, '[', 1)
    #pos <- sapply(pieces, '[', 2)
    #bar <- paste(slide, pos, sep = "_")
    

    bar <- sub('_Red.idat','', sub('_Grn.idat','',ifile))


    # This breaks when target_cpgs are out of index?
    mlu <- methylumIDATepic(barcodes = bar, n=n, force=force, ...)
    if(force){
        # TODO ensure that rows get forced to match original gds rownames
        # This is too unstable...
        recon <- names(assayData(mlu))
        results <- list()
        for(tabl in recon){
             results[[tabl]] <- t(t(assayDataElement(mlu, tabl)[,1][rown]))
        }
        assayData(mlu) <- results
        fData(mlu) <- fData(mlu)[rown,]
    }
    output <- app2gds(mlu, gds)
    return(output)
}

finalreport2gds <- function(finalreport, gds, ...){
    mset <- methylumiR(finalreport, ...)
    bmln <- es2gds(mset, gds)
    return(bmln)
} #}}}

# idats2gds -- idats2gds will add data from a set of barcodes into the same gds {{{
idats2gds <- function(barcodes, gds, n=TRUE, force=FALSE, ...){

    if(force){
        # Get a list of Greens.
        grns <- sprintf('%s_Grn.idat', barcodes)
        message('Determining IDAT lengths and ChipType')
        idx <- lapply(grns,function(x){
            f <- illuminaio::readIDAT(x)
            r <- rownames(f$Quants)
            return(list(chip=f$ChipType, rn=r))
        })
        chiptype <- unique(sapply(idx, '[[', 'chip'))
        if(length(chiptype) > 1) stop('idats2gds cannot process arrays from two separate platforms')
        manifest <- switch(chiptype,
            'BeadChip 8x5' = 'IlluminaHumanMethylationEPIC',
            'BeadChip 12x8' = 'IlluminaHumanMethylation450k',
            'BeadChip 12x1' = 'IlluminaHumanMethylation27k',
        )
        if (ChipType == "BeadChip 8x5" && length(idx$rn) > 1.1e+06) {
            manifest = "IlluminaHumanMethylationEpicv2"
        }
        manifest <- minfi::getManifest(manifest)
        cpgs_a <- cpgs_b <- c(getProbeInfo(manifest, type = "I")$Name, getProbeInfo(manifest, type = "II")$Name)
        names(cpgs_a) <- c(getProbeInfo(manifest, type = "I")$AddressA, getProbeInfo(manifest, type = "II")$AddressA)
        names(cpgs_b) <- c(getProbeInfo(manifest, type = "I")$AddressB, getProbeInfo(manifest, type = "II")$AddressB)
        message('Determining CpG intersection')
        tot_cpgs <- lapply(idx, function(y, cpgs_a, cpgs_b){
            unique(na.omit(c(cpgs_a[y[[2]]], cpgs_b[y[[2]]])))
        }, cpgs_a=cpgs_a, cpgs_b=cpgs_b)

        target_cpgs <- sort(Reduce(intersect, tot_cpgs))
        names(target_cpgs) <- target_cpgs
        message(sprintf('Expected Number of CpGs: %s (Number may be lower!)', length(target_cpgs)))
    }

    for(bar in barcodes){
        message(sprintf('Processing: %s...', bar))
        output <- iadd(bar, gds = gds, n = n, force = force, target_cpgs = target_cpgs, ...)
    }

    return(output)
} # }}}
