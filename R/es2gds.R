es2gds <- function(m, file, qc = TRUE){
    # Not memory efficient as we are working from
    # pre-loaded Expression Set derivatives.
    history.submitted <- as.character(Sys.time())
    if(!inherits(m, c('MethylSet', 'MethyLumiSet', 'RGChannelSet'))){
        stop(paste(deparse(substitute(m)),
            "is not a MethyLumiSet or MethylSet or RGChannelSet object"))
    }
    if(!(file_test('-f', file))){
        message(paste(file, "doesn't exist, creating new file"))
        bmln <- newgds(file)
    } else if ( (file_test('-f', file))){
        stop('File already exists')
    }
    # methylumi set method
    if(inherits(m, 'MethyLumiSet')){
        # check for existence of file
        # todo: check sample names for overlap
        if(exists("betas", assayData(m))){
            message('betas... ' )
            add.gdsn(bmln, "betas", val = betas(m))
        }
        if(exists("pvals", assayData(m))){
            message('pvals... ' )
            add.gdsn(bmln, "pvals", val = pvals(m))
        }
        if(exists("methylated", assayData(m))){
            message('methylated... ' )
            add.gdsn(bmln, "methylated", val = methylated(m))
        }
        if(exists("unmethylated", assayData(m))){
            message('unmethylated... ' )
            add.gdsn(bmln, "unmethylated", val = unmethylated(m))
        }
        if(exists("NBeads", assayData(m))){
            message('NBeads... ' )
            add.gdsn(bmln, "NBeads", val = assayDataElement(m, "NBeads"))
        }
        message('fData... ' )
        # deal with empty column names in fData
        nb_col <- which(colnames(fData(m)) == "")
        colnames(fData(m))[nb_col] <- nb_col
        add.gdsn(bmln, "fData", 
        val = data.frame(lapply(fData(m), as.character),
                stringsAsFactors = FALSE), check = FALSE)
        message('pData... ' )
        add.gdsn(bmln, "pData", 
        val = data.frame(lapply(pData(m), as.character),
                stringsAsFactors = FALSE))
        message('qcdata... ' )
        if(qc){
            if(!is.null(m@QC)){
                add.gdsn(bmln, "QCmethylated",
                        val = a <- methylated(QCdata(m)))
                add.gdsn(bmln, "QCunmethylated",
                        val = unmethylated(QCdata(m)))
                add.gdsn(bmln, "QCrownames", val = rownames(a))
            } else {
                message('   ... skipped')
            }
        }
        message('finished creating gdsfile' )
        add.gdsn(bmln, "history", 
        val = data.frame(lapply(m@history, as.character),
                        stringsAsFactors = FALSE))
        history.command <- "MethylumiSet converted to gds (bigmelon)"
        history.finished <- as.character(Sys.time())
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
    } 
    #BEGIN Minfi Methods
    # BEGIN MethylSet method
    if(inherits(m, 'MethylSet')){
        # Betas
        message('betas... ' )
        add.gdsn(bmln, "betas", val = getBeta(m))
        # M
        message('methylated... ' )
        add.gdsn(bmln, "methylated", val = getMeth(m))
        # U
        message('unmethylated... ' )
        add.gdsn(bmln, "unmethylated", val = getUnmeth(m))
        # fData
        message('fData... ')
        fd <- data.frame(lapply(getAnnotation(m), as.character),
                stringsAsFactors = FALSE)
        add.gdsn(bmln, "fData", val = fd, check = FALSE)
        message('pData... ')
        # Assuming data-less pData(m)
        pd <- data.frame(lapply(cbind(rownames(pData(m)), pData(m)),
                        as.character), stringsAsFactors = FALSE)
        colnames(pd) <- c("barcode", colnames(pData(m)))
        add.gdsn(bmln, "pData", val = pd)
        # QC
        if(qc){
            message('qcdata... ')
            message('    ... skipped')
        }
        message('finished creating gdsfile' )
        add.gdsn(bmln, "history", 
        val = data.frame(lapply(data.frame(submitted = "",
                                            finished = "",
                                            command = ""),
                            as.character), stringsAsFactors = FALSE))
        history.command <- "MethylSet converted to gds (bigmelon)"
        history.finished <- as.character(Sys.time())
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
    # END MethylSet Method
    # BEGIN RGChannelSet method
    }
    if(inherits(m, 'RGChannelSet')){
        m2 <- preprocessRaw(m)
        message('betas... ')
        add.gdsn(bmln, "betas", val = getBeta(m2))
        message('methylated... ')
        add.gdsn(bmln, "methylated", val = getMeth(m2))
        message('unmethylated... ')
        add.gdsn(bmln, "unmethylated", val = getUnmeth(m2))
        if(inherits(m, 'RGChannelSetExtended')){
            message('pvals... ')
            add.gdsn(bmln, "pvals", val = detectionP(m))
            message('NBeads... ')
            add.gdsn(bmln, "NBeads", val = as.matrix(beadcount(m)))
        }
        message('fData... ')
        fd <- data.frame(lapply(getAnnotation(m), as.character),
                            stringsAsFactors = FALSE)
        add.gdsn(bmln, "fData", val = fd,check = FALSE)
        message('pData... ')
        pd <- data.frame(lapply(cbind(rownames(pData(m)), pData(m)),
                            as.character), stringsAsFactors = FALSE)
        colnames(pd) <- c("barcode", colnames(pData(m)))
        add.gdsn(bmln, "pData", val = pd)
        message('qcdata... ')
        if(qc){
            ctrls <- getProbeInfo(m, type = "Control")
            ctrls <- ctrls[ctrls$Address %in% rownames(m), ]
            add.gdsn(bmln, "QCmethylated",
                        val = a <- getGreen(m)[ctrls$Address, ])
            add.gdsn(bmln, "QCunmethylated",
                        val = getRed(m)[ctrls$Address, ])
            add.gdsn(bmln, "QCrownames", val = rownames(a))
        } else { 
            message('     ...skipped')
        }
        message('finished creating gdsfile' )
        add.gdsn(bmln, "history", 
                    val = data.frame(lapply(data.frame(submitted = "",
                                                    finished = "",
                                                    command = ""),
                        as.character), stringsAsFactors = FALSE)
                )
        history.command <- "RGChannelSet converted to gds (bigmelon)"
        history.finished <- as.character(Sys.time())
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
    # END RGChannelSet(Extended) Method
    }
    # Point to rowcol
    # Search non-specifically for potential row and colnames
    rown <- list()
    # Change search terms - in any event!
    for(i in c("targetID", "name", "Probe_ID")){
        rown[[i]] <- grep(i, ls.gdsn(index.gdsn(bmln, "fData")),
                            ignore.case = TRUE)
    }
    rowna <- ls.gdsn(index.gdsn(bmln, "fData"))[unlist(rown)[1]]
    coln <- list()
    # Change search terms - in any event!
    for(i in c("sampleID", "barcode", "name")){ 
        coln[[i]] <- grep(i, ls.gdsn(index.gdsn(bmln, "pData")),
                            ignore.case = TRUE)
    }
    colna <- ls.gdsn(index.gdsn(bmln, "pData"))[unlist(coln)[1]]
    message('Directing \'rownames\' to ', paste0("fData/", rowna),
            ' by default, change with redirect.gds if incorrect.')
    message('Directing \'colnames\' to ', paste0("pData/", colna),
            ' by default, change with redirect.gds if incorrect.')
    add.gdsn(bmln, "paths", val = c(paste0("fData/", rowna),
                                    paste0("pData/", colna)))
    bmln
}
