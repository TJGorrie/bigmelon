# generate a gdsfile stub
newgds <- function(file){
    message('creating gdsfile...' )
    bmln <- createfn.gds(file)
    n <- add.gdsn(bmln, name = "description")
    put.attr.gdsn(n, "FileFormat", "DNA_Methylation")
    bmln
} # AS - OK

# function to check if file is already linked to by an existing gds object
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
} # AS - OK
# take a gds.class object, the name of a gds.class object or a gds filename
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
            handle <- openfn.gds(g$filename, readonly=FALSE)
        } else if(file_test('-f', gds)){ # On Disk
            #file exists, try to open
            # if this fails, check if it is already linked by a gds object
            handle <- tryCatch({
                openfn.gds(gds, readonly=FALSE)
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
}

# append one or more arrays of data to gdsfile
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
        if(! length(rownames(bmln)) == length(rownames(m))){
            stop(paste(deparse(substitute(m)), "has a different length to bmln"))
        }
        m <- m[rownames(bmln),]

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
    bmln
}

iadd2 <- function(path, gds, chunksize = NULL, force=TRUE,...){
    rown <- TRUE
    if(force){
        thefile <- try(openfn.gds(gds, allow.duplicate=TRUE), silent=T)
        if(!inherits('try-error', thefile)) rown <- rownames(thefile)
        closefn.gds(thefile)
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
            ml <- methylumIDATepic(barcodes = sets[!is.na(sets)],
                                    n = TRUE, oob = FALSE, idatPath = path)
            gdsfile <- app2gds(ml, gdsfile)
        }
    } else {
        ml <- methylumIDATepic(barcodes, n = TRUE,
                                oob = FALSE, idatPath = path)
        gdsfile <- app2gds(ml, gdsfile)
    }
    gdsfile
}

iadd <- function (bar, gds, n = TRUE, force=TRUE, ...){
    # TODO add check to ensure the arguments match to wacky dimension problems.
    rown <- TRUE
    if(force){
        thefile <- try(openfn.gds(gds, allow.duplicate=TRUE), silent=T)
        if(!inherits('try-error', thefile)) rown <- rownames(thefile)
        closefn.gds(thefile)
    }
    ifile <- basename(bar)
    pieces <- strsplit(ifile, "[_.]")
    slide <- sapply(pieces, '[', 1)
    pos <- sapply(pieces, '[', 2)
    bar <- paste(slide, pos, sep = "_")
    mlu <- methylumIDATepic(barcodes = bar, n=n, ...)[rown,]
    app2gds(mlu, gds)
}

finalreport2gds <- function(finalreport, gds, ...){
    mset <- methylumiR(finalreport, ...)
    bmln <- es2gds(mset, gds)
    return(bmln)
}
