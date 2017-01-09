redirect.gds<- function(gds, rownames, colnames){
    # Wrapper to change the gdsn node "paths" whererow and colnames are stored
    a <- try(index.gdsn(gds, rownames), silent = TRUE)
    if(inherits(a, "try-error")){
        stop(paste(rownames, "does not exist in gds object!"))
    }  
    b <- try(index.gdsn(gds, colnames), silent = TRUE)
    if(inherits(b, "try-error")){
        stop(paste(colnames, "does not exist in gds object!"))  
    }
    add.gdsn(gds, name = "paths", replace = TRUE, val = c(rownames, colnames))
    cat("Changing 'rownames' path to:", rownames, ". \n")
    cat("Changing 'colnames' path to:", colnames, ". \n") 
}

backup.gdsn <- function(gds = NULL, node){
    # Quick function to quickly copy specific nodes into a new-folder 
    # called "backup". not recommended for standard workflow. see copyto.gdsn
    if(is.null(gds)) gds <- getfolder.gdsn(node)
    if(!("backup" %in% ls.gdsn(gds))) addfolder.gdsn(gds, "backup")
    copyto.gdsn(index.gdsn(gds, "backup"), node)
}

'[.gds.class' <- function(x, i, j, node, name = TRUE, drop = TRUE){
    # Method of subsetting gds.class objects while specifying node.
    # Create gdsn.class for selected node
    dat <- index.gdsn(x, node)
    # Pass to '[.gdsn.class'
    dat[i = i, j = j, name = name, drop = drop]
} # }}}

#'[.gdsn.class' <- function(x, i, j, name = TRUE, drop = TRUE){ # {{{
    # Method of subsetting gdsn.class objects w.o reading in entire object.
    # Wrapper for readex.gdsn
    # arg: "name": Will point towards "fData/Probe_ID" and
    # "pData/barcode" for row and col names. (by default)
    #              name = F will enable for faster indexing.
    # TODO: [] method? a.k.a vector method.
    # TODO: Warnings, if any.
#    base <- getfolder.gdsn(x)
    # x[ , ]
#    if(missing(i) & missing(j)){ # {{{
#        mat <- read.gdsn(x)
#        if(name){
#            rownames(mat) <- read.gdsn(index.gdsn(base,
#                                read.gdsn(index.gdsn(base, "paths"))[1]))
#            colnames(mat) <- read.gdsn(index.gdsn(base,
#                                read.gdsn(index.gdsn(base, "paths"))[2]))
#        }
#    } # }}}
    # x[ , j]
#    if(missing(i) & !missing(j)){ # {{{
#        j1 <- j
#        if(is.character(j1)){
#            j <- match(j1, read.gdsn(index.gdsn(base,
#                            read.gdsn(index.gdsn(base, "paths"))[2]))) # ok
#        }
#        if(is.logical(j1)) j <- (1:objdesp.gdsn(x)$dim[2])[j1] # ok
#        ncol <- j
#        mat <- as.matrix(readex.gdsn(x, sel = list(NULL, ncol)))
#        if(length(j) == 1) mat <- as.matrix(mat)
#        if(name){
#            rownames(mat) <- read.gdsn(index.gdsn(base,
#                                read.gdsn(index.gdsn(base, "paths"))[1]))
#            colnames(mat) <- readex.gdsn(index.gdsn(base,
#                                read.gdsn(index.gdsn(base, "paths"))[2]),
#                                sel = ncol)
#        }
#    } # }}}
#    # x[ i, ]
#    if(!missing(i) & missing(j)){  # {{{
#        i1 <- i
#        if(is.character(i1)){
#            i <- match(i1, read.gdsn(index.gdsn(base,
#                            read.gdsn(index.gdsn(base, "paths"))[1]))) # ok
#        }
#        if(is.logical(i1)) i <- (1:objdesp.gdsn(x)$dim[2])[i1] # ok
#        nrow <- i
#        if(length(i) == 1){ # Calling a single row makes naming difficult.
#            mat <- as.matrix(t(readex.gdsn(x, sel = list(nrow, NULL))))
#       } else {
#            mat <- as.matrix(readex.gdsn(x, sel = list(nrow, NULL)))
#        }
#        if(name){
#            rownames(mat) <- readex.gdsn(index.gdsn(base,
#                                read.gdsn(index.gdsn(base, "paths"))[1]),
#                                            sel = nrow)
#            colnames(mat) <- read.gdsn(index.gdsn(base,
#                                read.gdsn(index.gdsn(base, "paths"))[2]))
#        }
#    } # }}}
    # x[ i, j] 
#    if(!missing(i) & !missing(j)){ # {{{
#        i1 <- i
#        if(is.character(i1)){
#            i <- match(i1, read.gdsn(index.gdsn(base,
#                            read.gdsn(index.gdsn(base, "paths"))[1]))) # ok
#        }
#        if(is.logical(i1))   i <- (1:objdesp.gdsn(x)$dim[2])[i1] # ok
#        j1 <- j
#        if(is.character(j1)){
#            j <- match(j1, read.gdsn(index.gdsn(base,
#                            read.gdsn(index.gdsn(base, "paths"))[2]))) # ok
#        }
#        if(is.logical(j1)) j <- (1:objdesp.gdsn(x)$dim[2])[j1] # ok
#        nrow <- i
#        ncol <- j
#        mat <- readex.gdsn(x, sel = list(nrow, ncol))
#        # i=1, j=1
#        if( length(i) == 1 & length(j)  == 1) mat <- matrix(mat)
#        # i>1, j=1
#        if(!length(i) == 1 & length(j)  == 1) mat <- matrix(mat)
#        # i=1, j>1
#        if( length(i) == 1 & !length(j) == 1) mat <- t(matrix(mat)) 
#        if(name){
#            rownames(mat) <- readex.gdsn(index.gdsn(base,
#                                read.gdsn(index.gdsn(base, "paths"))[1]),
#                                        sel = nrow)
#            colnames(mat) <- readex.gdsn(index.gdsn(base,
#                                read.gdsn(index.gdsn(base, "paths"))[2]),
#                                        sel = ncol)   
#        }
#    } # }}}
#    return(mat[ , , drop = drop])
#} # }}}

#colnames <- function (x, do.NULL = TRUE, prefix = "col") 
setMethod(
    f = "colnames",
    signature(x = "gds.class"),
    definition = function(x, do.NULL = TRUE, prefix = NULL){
        read.gdsn(index.gdsn(x,read.gdsn(index.gdsn(x, "paths"))[2]))
    }
) # OK

setMethod(
    f = "colnames",
    signature(x = "gdsn.class"),
    definition = function(x, do.NULL = TRUE, prefix = NULL){ 
        read.gdsn(index.gdsn(getfolder.gdsn(x),
            read.gdsn(index.gdsn(getfolder.gdsn(x), "paths"))[2]))
    }
)

#rownames <- function (x, do.NULL = TRUE, prefix = "col") 
setMethod(
    f = "rownames",
    signature(x = "gds.class"),
    definition = function(x, do.NULL = TRUE, prefix = NULL){
        read.gdsn(index.gdsn(x,read.gdsn(index.gdsn(x, "paths"))[1]))
    }
) # OK

setMethod(
    f = "rownames",
    signature(x = "gdsn.class"),
    definition = function(x, do.NULL = TRUE, prefix = NULL){
        read.gdsn(index.gdsn(getfolder.gdsn(x),
            read.gdsn(index.gdsn(getfolder.gdsn(x), "paths"))[1]))
    }
)

# Standard eset functions to grab matrices, including indexing:
# Calling gdsn node but will allow direct subsetting with '['
# Alternative work around would be setGeneric("betas") function(object, ...)
## setGeneric("betas", function(object, ...) standardGeneric("betas"))
#  object[i, j, node = "betas", name = TRUE, drop = FALSE]
# which overwrites other "betas" methods. 
# This also applies to other eset methods.

setMethod(
    f = "betas",
    signature(object = "gds.class"),
    definition = function(object){  
        index.gdsn(object, 'betas')
    }
) # OK

setMethod(
    f = "methylated",
    signature(object = "gds.class"),
    definition = function(object){
        index.gdsn(object, 'methylated')
    }
) # OK

setMethod(
    f = "unmethylated",
    signature(object = "gds.class"),
    definition = function(object){
        index.gdsn(object, 'unmethylated')
    }
) # OK

setMethod(
    f = "pvals",
    signature(object = "gds.class"),
    definition = function(object){
        index.gdsn(object, 'pvals')
    }
) # OK

# These are small and can be read completely into memory.
setMethod(
    f = "fData",
    signature(object = "gds.class"),
    definition = function(object){
        fd <- index.gdsn(object, 'fData')[ , ,name = FALSE, drop = FALSE]
        rownames(fd) <- rownames(object)
        fd
    }
) # OK

setMethod(
    f = "pData",
    signature(object = "gds.class"),
    definition = function(object){
        index.gdsn(object, 'pData')[ , , name = FALSE, drop = FALSE]
    }
) # OK

setMethod(
    f = "getHistory",
    signature(object = "gds.class"),
    definition = function(object){
        read.gdsn(index.gdsn(object, 'history'))
    }
) # AS - OK 

setGeneric("QCmethylated", function(object){standardGeneric("QCmethylated")})
setMethod(
    f = "QCmethylated",
    signature(object = "gds.class"),
    definition = function(object){
        out <- read.gdsn(index.gdsn(object, 'QCmethylated'))
        rownames(out) <- QCrownames(object)
        out
    }
) # AS - OK 

setGeneric("QCunmethylated",function(object){
    standardGeneric("QCunmethylated")
})
setMethod(
    f = "QCunmethylated",
    signature(object = "gds.class"),
    definition = function(object){
        out <- read.gdsn(index.gdsn(object, 'QCunmethylated'))
        rownames(out) <- QCrownames(object)
        out
    }
) # AS - OK 

setGeneric("QCrownames", function(object){standardGeneric("QCrownames")})
setMethod(
    f = "QCrownames",
    signature(object = "gds.class"),
    definition = function(object){
        read.gdsn(index.gdsn(object, 'QCrownames'))
    }
) # AS - OK  

setMethod(
    f = "betaqn",
    signature(bn = "gds.class"),
    definition = function(bn){
    if(bn$readonly) stop("gds object is in Read-Only mode, please reload!")
        # Create Temp Nod
        f <- createfn.gds("temp.gds", allow.duplicate = TRUE)
        history.submitted <- as.character(Sys.time())
        # Normalise betas into Temporary object
        qn.gdsn(gds = f,
                target = index.gdsn(bn, "betas"),
                newnode = "betas")
        # Create new node in original - replacing old node.
        n.t <- add.gdsn(bn,
                        name = "betas",
                        valdim = c(objdesp.gdsn(
                                    index.gdsn(f, "betas"))$dim[1], 0),
                        val = NULL,
                        storage = "float64",
                        replace = TRUE)
        # Append new betas to original gds col by col.
        for(i in 1:dim[2]){
            val <- readex.gdsn(index.gdsn(f, "betas"), sel = list(NULL, i))
            append.gdsn(node = n.t, val = val)
        }
        history.finished <- as.character(Sys.time())
        history.command <- "Normalized with betaqn method (bigmelon)"
        h <- data.frame(submitted = history.submitted,
                        finished = history.finished,
                        command = history.command,
                        stringsAsFactors = FALSE)
        for(j in colnames(h)){
            h_coln <- j
            h_index_str <- paste("history/", j, sep = "")
            h_child_n <- index.gdsn(bn, h_index_str, silent = TRUE)
            append.gdsn(h_child_n, val = h[,h_coln])  
        }
        # Deleting Temp File.
        closefn.gds(f)
        unlink("temp.gds", force = TRUE) 
    }
)

setMethod(
    f = "naten",
    signature(mn = "gds.class"),
    definition = function(mn, fudge = 100, ret2 = FALSE, node = "betas"){
    if(mn$readonly) stop("gds object is in Read-Only mode, please reload!")
        object <- mn
        history.submitted <- as.character(Sys.time())
        # Normalize using gds method - temp node created within!
        naten.gds(gds = object,
                    node,
                    mns = index.gdsn(object, "methylated"),
                    uns = index.gdsn(object, "unmethylated"),
                    fudge,
                    ret2
                    )
        history.finished <- as.character(Sys.time())
        history.command <- "Normalized with naten method (bigmelon)"
        h <- data.frame(submitted = history.submitted, 
                        finished = history.finished, 
                        command = history.command,
                        stringsAsFactors = FALSE)
        for(j in colnames(h)){
            h_coln <- j
            h_index_str <- paste("history/",j,sep="")
            h_child_n <- index.gdsn(object,h_index_str,silent=TRUE)
            append.gdsn(h_child_n, val=h[,h_coln])
        }
    }
) # TGS

setMethod(
    f = "nanet",
    signature(mn = "gds.class"),
    definition = function(mn, fudge = 100, ret2 = FALSE, node = "betas"){
        history.submitted <- as.character(Sys.time())
        if(mn$readonly) stop("gds object in Read-Only mode, please reload!")
        object <- mn
        nanet.gds(gds = object,
                    node,
                    mns = index.gdsn(object, "methylated"),
                    uns = index.gdsn(object, "unmethylated"),
                    fudge,
                    ret2
                    )
        history.finished <- as.character(Sys.time())
        history.command <- "Normalized with nanet method (wateRmelon)"
        h <- data.frame(submitted = history.submitted, 
                        finished = history.finished, 
                        command = history.command,
                        stringsAsFactors = FALSE
                        )
        for(j in colnames(h)){
            h_coln <- j
            h_index_str <- paste("history/",j,sep="")
            h_child_n <- index.gdsn(object,h_index_str,silent=TRUE)
            append.gdsn(h_child_n, val=h[,h_coln])
        }
    }
)

setMethod(
    f = "nanes",
    signature(mns = "gds.class"),
    definition = function(mns, fudge=100, ret2=FALSE, node = "betas", ...){
        history.submitted <- as.character(Sys.time())
        if(mns$readonly) stop("gds object in Read-Only mode")   
        object <- mns
        nanes.gds(gds = object,
                    node,
                    mns = index.gdsn(object, "methylated"),
                    uns = index.gdsn(object, "unmethylated"),
                    onetwo = fData(object)[, grep('DESIGN',
                                                colnames(fData(object)),
                                                ignore.case = TRUE)[1]],
                    fudge,
                    ret2
                    )
        history.finished <- as.character(Sys.time())
        history.command <- "Normalized with nanes method (wateRmelon)"
        h <- data.frame(submitted = history.submitted, 
                        finished = history.finished, 
                        command = history.command,
                        stringsAsFactors = FALSE
                        )
        for(j in colnames(h)){
            h_coln <- j
            h_index_str <- paste("history/",j,sep="")
            h_child_n <- index.gdsn(object,h_index_str,silent=TRUE)
            append.gdsn(h_child_n, val=h[,h_coln])
        }
    }
) # TGS

setMethod(
    f = "danes",
    signature(mn = "gds.class"),
    definition = function(mn, fudge = 100, ret2 = FALSE, node = "betas",...){
        history.submitted <- as.character(Sys.time())  
        if(mn$readonly) stop("gds object in Read-Only mode, please reload!")
        object <- mn
        danes.gds(gds = object,
                    node, 
                    mns = index.gdsn(object, "methylated"),
                    uns = index.gdsn(object, "unmethylated"),
                    onetwo = fData(object)[, grep('DESIGN',
                                                colnames(fData(object)),
                                                ignore.case = TRUE)[1]],
                    fudge,
                    ret2
                    )
        history.finished <- as.character(Sys.time())
        history.command <- "Normalized with danes method (wateRmelon)"
        h <- data.frame(submitted = history.submitted, 
                        finished = history.finished, 
                        command = history.command,
                        stringsAsFactors = FALSE
                        )
        for(j in colnames(h)){
            h_coln <- j
            h_index_str <- paste("history/",j,sep="")
            h_child_n <- index.gdsn(object,h_index_str,silent=TRUE)
            append.gdsn(h_child_n, val=h[,h_coln])
        }  
    }
) 

setMethod(
    f = "danet",
    signature(mn = "gds.class"),
    definition = function(mn, fudge = 100, ret2 = FALSE, node = "betas", ...){
        history.submitted <- as.character(Sys.time())   
        if(mn$readonly) stop("gds object in Read-Only mode, please reload!")
        object <- mn 
        danet.gds(  gds = object,
                    node,
                    mns = index.gdsn(object, "methylated"),
                    uns = index.gdsn(object, "unmethylated"),
                    onetwo = fData(object)[,grep('DESIGN',
                                                colnames(fData(object)),
                                                ignore.case = TRUE)[1]],
                    roco = unlist(data.frame(
                            strsplit(colnames(object), '_'),
                            stringsAsFactors = FALSE)[2,]),
                    fudge,
                    ret2
                    )
        history.finished <- as.character(Sys.time())
        history.command <- "Normalized with danet method (wateRmelon)"
        h <- data.frame(submitted = history.submitted, 
                        finished = history.finished, 
                        command = history.command,
                        stringsAsFactors=FALSE
                        )
        for(j in colnames(h)){
            h_coln <- j
            h_index_str <- paste("history/",j,sep="")
            h_child_n <- index.gdsn(object,h_index_str,silent=TRUE)
            append.gdsn(h_child_n, val=h[,h_coln])
        }
    }
) # TGS

setMethod(
    f = "daten1",
    signature(mn = "gds.class"),
    definition = function(mn, fudge = 100, ret2 = FALSE, node = "betas", ...){
        history.submitted <- as.character(Sys.time()) 
        if(mn$readonly) stop("gds object in Read-Only mode, please reload!")
        object <- mn  
        daten1.gds(gds = object,
                    node,
                    mns = index.gdsn(object, "methylated"),
                    uns = index.gdsn(object, "unmethylated"),
                    onetwo = fData(object)[,grep('DESIGN',
                                                colnames(fData(object)),
                                                ignore.case = TRUE)[1]],
                    roco = unlist(data.frame(
                                strsplit(colnames(object), '_'),
                                stringsAsFactors = FALSE)[2,]),
                    fudge,
                    ret2
                    )
        history.finished <- as.character(Sys.time())
        history.command <- "Normalized with daten1 method (wateRmelon)"
        h <- data.frame(submitted = history.submitted, 
                        finished = history.finished, 
                        command = history.command,
                        stringsAsFactors = FALSE
                        )
        for(j in colnames(h)){
            h_coln <- j
            h_index_str <- paste("history/",j,sep="")
            h_child_n <- index.gdsn(object,h_index_str,silent=TRUE)
            append.gdsn(h_child_n, val=h[,h_coln])
        }
    }
)

setMethod(
    f = "daten2",
    signature(mn = "gds.class"),
    definition = function(mn, fudge = 100,ret2 = FALSE,node = "betas", ...){
        history.submitted <- as.character(Sys.time())
        if(mn$readonly) stop("gds object in Read-Only mode, please reload!")
        object <- mn 
        daten2.gds(gds = object,
                    node,
                    mns = index.gdsn(object, "methylated"),
                    uns = index.gdsn(object, "unmethylated"),
                    onetwo =  fData(object)[,grep('DESIGN',
                                            colnames(fData(object)), 
                                            ignore.case = TRUE)[1]],
                    roco = unlist(data.frame(
                                strsplit(colnames(object), '_'),
                                stringsAsFactors = FALSE)[2,]),
                    fudge,
                    ret2
                    )
        history.finished <- as.character(Sys.time())
        history.command <- "Normalized with daten2 method (wateRmelon)"
        h <- data.frame(submitted = history.submitted, 
                        finished = history.finished, 
                        command = history.command,
                        stringsAsFactors = FALSE
                        )
        for(j in colnames(h)){
            h_coln <- j
            h_index_str <- paste("history/",j,sep="")
            h_child_n <- index.gdsn(object,h_index_str,silent=TRUE)
            append.gdsn(h_child_n, val=h[,h_coln])
        }
    }
)

setMethod(
    f = "nasen",
    signature(mns = "gds.class"),
    definition = function(mns, ret2 = FALSE, fudge = 100, node = "betas"){
        history.submitted <- as.character(Sys.time())
        if(mns$readonly) stop("gds object in Read-Only mode, please reload!")
        object <- mns
        nasen.gds(gds = object,
                    node,
                    mns = index.gdsn(object, "methylated"),
                    uns = index.gdsn(object, "unmethylated"),
                    onetwo =  fData(object)[,grep('DESIGN',
                                            colnames(fData(object)), 
                                            ignore.case = TRUE)[1]],
                    fudge,
                    ret2
                    )
        history.finished <- as.character(Sys.time())
        history.command <- "Normalized with nasen method (wateRmelon)"
        h <- data.frame(submitted = history.submitted, 
                        finished = history.finished, 
                        command = history.command,
                        stringsAsFactors = FALSE
                        )
        for(j in colnames(h)){
            h_coln <- j
            h_index_str <- paste("history/",j,sep="")
            h_child_n <- index.gdsn(object,h_index_str,silent=TRUE)
            append.gdsn(h_child_n, val=h[,h_coln])
        }
    }
)

setMethod(
    f = "dasen",
    signature(mns = "gds.class"),
    definition = function(mns, fudge = 100, ret2 = FALSE, node ="betas", ...){
        history.submitted <- as.character(Sys.time())   
        if(mns$readonly) stop("gds object in Read-Only mode, please reload!")
        object <- mns
        dasen.gds(object,
                    node,
                    mns = index.gdsn(object, "methylated"),
                    uns = index.gdsn(object, "unmethylated"),
                    onetwo =  fData(object)[,grep('DESIGN',
                                                colnames(fData(object)), 
                                                ignore.case = TRUE)[1]],
                    roco = unlist(data.frame(
                                strsplit(colnames(object), '_'),
                                stringsAsFactors = FALSE)[2,]),
                    fudge,
                    ret2
                    )
        history.finished <- as.character(Sys.time())
        history.command <- "Normalized with dasen method (wateRmelon)"
        h <- data.frame(submitted = history.submitted, 
                        finished = history.finished, 
                        command = history.command,
                        stringsAsFactors=FALSE
                        )
        for(j in colnames(h)){
            h_coln <- j
            h_index_str <- paste("history/",j,sep="")
            h_child_n <- index.gdsn(object,h_index_str,silent=TRUE)
            append.gdsn(h_child_n, val=h[,h_coln])
        }
    }
)

setMethod(
    f = "danen",
    signature(mns = "gds.class"),
    definition = function(mns, fudge = 100, ret2 = FALSE, node = "betas",...){
        history.submitted <- as.character(Sys.time())
        if(mns$readonly) stop("gds object in Read-Only mode, please reload!")
        object <- mns
        danen.gds(gds = object,
                    node,
                    mns = index.gdsn(object, "methylated"),
                    uns = index.gdsn(object, "unmethylated"),
                    onetwo =  fData(object)[,grep('DESIGN',
                                            colnames(fData(object)),
                                            ignore.case = TRUE)[1]],
                    roco = unlist(data.frame(
                                strsplit(colnames(object), '_'),
                                stringsAsFactors=FALSE)[2,]),
                    fudge,
                    ret2
                    )
        history.finished <- as.character(Sys.time())
        history.command <- "Normalized with danen method (wateRmelon)"
        h <- data.frame(submitted = history.submitted, 
                        finished = history.finished, 
                        command = history.command,
                        stringsAsFactors=FALSE
                        )
        for(j in colnames(h)){
            h_coln <- j
            h_index_str <- paste("history/",j,sep="")
            h_child_n <- index.gdsn(object,h_index_str,silent=TRUE)
            append.gdsn(h_child_n, val=h[,h_coln])
        }  
    }
)

setMethod(
    f = "exprs",
    signature(object = "gds.class"),
    definition = function(object){
        exp <- data.frame(betas(object)[,, name = TRUE], 
                            row.names = rownames(object),
                            check.rows = FALSE,
                            check.names = FALSE,
                            stringsAsFactors = FALSE)
        colnames(exp) <- colnames(object)
        exp
    }
) # AS - OK  # TGS -

setMethod(
    f = "pfilter",
    signature(mn = "gds.class"),
    definition = function(mn, perCount = NULL,
                            pnthresh = NULL, perc = NULL, pthresh = NULL ){
        if(mn$readonly) stop("gds object in Read-Only mode, please reload!")
        history.submitted <- as.character(Sys.time())
        object <- mn
        if("NBeads" %in% ls.gdsn(mn)){ 
            nb <- index.gdsn(object, "NBeads")
        } else { 
            cat("NBeads missing, using betas instead... \n")
            nb   <- index.gdsn(object, "betas")
        }
        bc    <- nb
        mn    <- NULL
        bn    <- NULL
        un    <- NULL
        pn    <- pvals(object)
        da    <- NULL
        l    <- pfilter.gds(mn = mn, un = un, bn = bn, da = da, 
                                pn = pn, bc = bc, perCount, pnthresh, perc,
                                pthresh) 
        history.finished <- as.character(Sys.time())
        history.command <- "pfilter applied (bigmelon)"
        h <- data.frame(submitted = history.submitted,
                        finished = history.finished,
                        command = history.command,
                        stringsAsFactors=FALSE
                        )
        for(j in colnames(h)){
            h_coln <- j
            h_index_str <- paste("history/",j,sep="")
            h_child_n <- index.gdsn(object,h_index_str,silent=TRUE)
            append.gdsn(h_child_n, val=h[,h_coln])
        }
        lpro <- l$probes
        lsam <- l$samples
        subSet(object,which(lpro),which(lsam))
    }  
)   # AS - OK # TGS - OK

#subsetting : i = features, j = samples
# TODO: Test Logical Indexing and Character Subsetting also(?)
setGeneric("subSet",function(x,i,j,...,drop=FALSE){standardGeneric("subSet")})
setMethod(
    f = "subSet",
    signature(x = "gds.class"),
    definition = function(x, i, j, ..., drop = FALSE)  {
        if(x$readonly) stop("gds object is in Read-Only mode, please reload!")
        history.submitted <- as.character(Sys.time())
        nodules <- ls.gdsn(x) # Important!
        f <- createfn.gds("temp.gds", allow.duplicate = TRUE)
        if("betas" %in% nodules){
            # Copy node to temp node colbycol 
            trait <- tolower(objdesp.gdsn(betas(x))$trait)
            n.t <- add.gdsn(f, name = "betas",
                            valdim = c(objdesp.gdsn(betas(x))$dim[1],0),
                            val = NULL, storage = trait)
            for(a in 1:objdesp.gdsn(betas(x))$dim[2]){
                append.gdsn(node = n.t, val = readex.gdsn(betas(x),
                            sel=list(NULL, a)))
            }
            # Replace Old Node with Subset colbycol
            n.n <- add.gdsn(x, name = "betas", valdim = c(length(i), 0),
                            val = NULL, replace = TRUE, storage = trait)
            for(z in j){
                append.gdsn(node = n.n, val = readex.gdsn(betas(f),
                            sel = list(i, z)))
            }
        }

        if("pvals" %in% nodules){
            trait <- tolower(objdesp.gdsn(pvals(x))$trait) 
            n.t <- add.gdsn(f, name = "pvals",
                            valdim = c(objdesp.gdsn(pvals(x))$dim[1],0),
                            val = NULL, storage = trait)
            for(a in 1:objdesp.gdsn(pvals(x))$dim[2]){
                append.gdsn(node = n.t, val = readex.gdsn(pvals(x),
                            sel = list(NULL, a)))
            }
            # Replace Old Node with Subset colbycol
            n.n <- add.gdsn(x, name = "pvals", valdim = c(length(i), 0),
                            val = NULL, replace = TRUE, storage = trait)
            for(z in j){
                append.gdsn(node = n.n, val = readex.gdsn(pvals(f),
                            sel = list(i, z)))
            }  
        }

        if("methylated" %in% nodules){
            trait <- tolower(objdesp.gdsn(methylated(x))$trait)
            n.t <- add.gdsn(f, name="methylated", 
                            valdim=c(objdesp.gdsn(methylated(x))$dim[1],0),
                            val = NULL, storage = trait)
            for(a in 1:objdesp.gdsn(methylated(x))$dim[2]){
                append.gdsn(node = n.t, val = readex.gdsn(methylated(x),
                            sel=list(NULL,a)))
            }
            # Replace Old Node with Subset colbycol
            n.n <- add.gdsn(x, name="methylated", valdim=c(length(i), 0),
                            val=NULL, replace=TRUE, storage = trait)
            for(z in j){
                append.gdsn(node = n.n, val = readex.gdsn(methylated(f),
                            sel = list(i, z)))
            }  
        }

        if("unmethylated" %in% nodules){
            trait <- tolower(objdesp.gdsn(unmethylated(x))$trait)
            n.t <- add.gdsn(f, name="unmethylated", 
                            valdim=c(objdesp.gdsn(unmethylated(x))$dim[1],0),
                            val = NULL, storage = trait)
            for(a in 1:objdesp.gdsn(unmethylated(x))$dim[2]){
                append.gdsn(node = n.t, val = readex.gdsn(unmethylated(x),
                            sel=list(NULL,a)))
            }
            # Replace Old Node with Subset colbycol
            n.n <- add.gdsn(x, name="unmethylated", valdim=c(length(i), 0),
                            val = NULL, replace = TRUE, storage = trait)
            for(z in j){
            append.gdsn(node = n.n, val = readex.gdsn(unmethylated(f),
                        sel=list(i, z)))
            } 
        }

        if("NBeads" %in% nodules){
            nb <- index.gdsn(x, "NBeads")
            trait <- tolower(objdesp.gdsn(nb)$trait)
            n.t <- add.gdsn(f, name="NBeads",
                            valdim=c(objdesp.gdsn(nb)$dim[1],0),
                            val = NULL, storage = trait)
            for(a in 1:objdesp.gdsn(nb)$dim[2]){
                append.gdsn(node = n.t,
                            val = readex.gdsn(nb, sel=list(NULL,a)))
            }
            # Replace Old Node with Subset colbycol
            nb2 <- index.gdsn(f, "NBeads")
            n.n <- add.gdsn(x, name="NBeads", valdim=c(length(i), 0),
                            val=NULL, replace=TRUE, storage = trait)
            for(z in j){
                append.gdsn(node = n.n,
                            val = readex.gdsn(nb2, sel=list(i, z)))
            } 
        }
        # These are small enough to warrant not worrying about memory use.
        if("fData" %in% nodules){  
            fdatasubs <- fData(x)[i,,drop=FALSE]
            add.gdsn(x, name="fData", valdim=dim(fdatasubs),
                    val=fdatasubs, replace=TRUE)  
        }

        if("pData" %in% nodules){  
            pdatasubs <- pData(x)[j,,drop=FALSE]
            add.gdsn(x, name="pData", valdim=dim(pdatasubs), 
                    val=pdatasubs, replace=TRUE)  
        }

        if("QCmethylated" %in% nodules){  
            qcmethsubs <- QCmethylated(x)[,j]
            add.gdsn(x, name="QCmethylated", valdim=dim(qcmethsubs), 
                    val= qcmethsubs, replace=TRUE)  
        }

        if("QCunmethylated" %in% nodules){  
            qcumethsubs <- QCunmethylated(x)[,j]
            add.gdsn(x, name="QCunmethylated", valdim=dim(qcumethsubs), 
                    val = qcumethsubs, replace=TRUE)
        }

        history.finished <- as.character(Sys.time())
        dim <- objdesp.gdsn(betas(x))$dim
        history.command <- paste("Subset of", dim[1],
                                "rows and", dim[2], "samples")
        h <- data.frame(submitted = history.submitted,
                        finished = history.finished,
                        command = history.command,
                        stringsAsFactors=FALSE
                        )
        for(j in colnames(h)){
            h_coln <- j
            h_index_str <- paste("history/",j,sep="")
            h_child_n <- index.gdsn(x,h_index_str,silent=TRUE)
            append.gdsn(h_child_n, val=h[,h_coln])
        }
        closefn.gds(f)
        unlink("subsettemp.gds")
    }
) # AS - OK # TGS OKish

#prcomp prcompGdsn contains S3 method!
#setMethod(
#  f = "prcomp",
#  signature(x = "gdsn.class"),
#  definition = function(x, retx = FALSE, center = FALSE, scale. = FALSE,
#                        perc = 0.01, npcs=5, ...){
#    prcomp.gdsn(x, retx, center, scale., perc, npcs)
#  }
#)

# outlyx - instead of utilizing prcomp.gdsn (since pcout isn't bigmelon 
# friendly - we opt for generic use of prcomp)
setMethod(
    f = "outlyx",
    signature(x = "gdsn.class"),
    definition = function(x, iqr = TRUE, iqrP = 2, pc=1, mv = TRUE,
                            mvP = 0.15, plot = TRUE, perc = 0.01, ...){
        dimx <- objdesp.gdsn(x)$dim
        outlyx(x[sample(1:dimx[1], dimx[1]*perc, replace = FALSE), ],
                iqr, iqrP, pc, mv, mvP, plot)
    }
)

setMethod( # Automatically select betas
    f = "outlyx",
    signature(x = "gds.class"),
    definition = function(x, iqr = TRUE, iqrP = 2, pc = 1, mv = TRUE, 
                            mvP = 0.15, plot = TRUE, perc = 0.01, ...){
        x <- betas(x)
        dimx <- objdesp.gdsn(x)$dim
        outlyx(x[sample(1:dimx[1], dimx[1]*perc, replace=FALSE), ],
                iqr, iqrP, pc, mv, mvP, plot)
    }
)

# agep # Extract coeff names. #Acts odd on marmalaid? 
setMethod(
    f = "agep",
    signature(betas = "gds.class"),
    definition = function(betas, coeff = NULL, verbose = FALSE){
        betas <- index.gdsn(betas, "betas")
        if(is.null(coeff)){ 
            data(coef)
            coeff <- coef 
        }
        betas <- betas[names(coeff)[-1], , name = TRUE, drop = FALSE]
        rownames(betas) <- names(coeff)[-1]
        # (Violently) Produces "small beta matrix".
        agep(betas, coeff, verbose) 
    }
)

setMethod(
    f = "agep",
    signature(betas = "gdsn.class"),
    definition = function(betas, coeff = NULL, verbose = FALSE){
        if(is.null(coeff)){
            data(coef)
            coeff <- coef
        }
        betas <- betas[names(coeff)[-1], , name = TRUE, drop = FALSE]
        rownames(betas) <- names(coeff)[-1]
        agep(betas, coeff, verbose)
    }
)

setGeneric(name= "qual")
# qual (?) # Do col by col computation of metrics. Collapse output.
setMethod(
    f= "qual",
    signature(norm="gdsn.class", raw="gdsn.class"),
    definition = function(norm, raw){
        dimnorm <- objdesp.gdsn(norm)$dim
        dimraw  <- objdesp.gdsn(raw)$dim
        if(!all(dimnorm == dimraw)) stop("Nodes are not the same dimensions")
        res <- t(sapply(1:dimnorm[2], function(x){
            dif <- norm[,x,name=F] - raw[,x,name=F]
            rmsd <- sqrt(mean(dif^2, na.rm = TRUE))
            sdd  <- sd(dif, na.rm = TRUE)
            sadd <- sd(abs(dif), na.rm = TRUE)
            srms <- rmsd/sdd
            out <- c(rmsd, sdd, sadd, srms)
            out
        } ) )
        rownames(res) <- colnames(norm)
        colnames(res) <- c("rmsd", "sdd", "sadd", "srms")
        res
    }
)

# bscon see bscon_methy.R
setMethod(
    f= "bscon",
    signature(x = "gds.class"),
    definition = function(x){
        nodules <- ls.gdsn(x)
        if("QCmethylated"%in%nodules){
            green.Channel <- QCmethylated(x)
        } else {
            stop("Green channel QC data could not be found")
        }
        if("QCunmethylated"%in%nodules){
            red.Channel <- QCunmethylated(x)
        } else {
            stop("Red channel QC data could not be found")
        }
        QCrows <- QCrownames(x)
        bisulfite.green <- green.Channel[grep("^B.*C.*I", QCrows),] 
        bisulfite.red <- red.Channel[grep("^B.*C.*I", QCrows),]
        bsI.green <- bisulfite.green[1:12,]
        bsI.red <- bisulfite.red[1:12,]
        bsII.green <- green.Channel[grep( '^B.*C.*II', QCrows),]
        bsII.red <- red.Channel[grep( '^B.*C.*II', QCrows),]
        BSI.betas <- rbind( bsI.green[1:3,], bsI.red[7:9,])/((rbind(
                            bsI.green[1:3,], bsI.red[7:9,])) + rbind(
                            bsI.green[4:6,], bsI.red[10:12,]))
        BSII.betas <- bsII.red/(bsII.red + bsII.green)
        apply(rbind(BSI.betas, BSII.betas), 2, median)*100 
    }
)

setMethod(
    f = "pwod",
    signature(object = "gds.class"),
    definition = function(object, mul = 4){
        if(object$readonly) stop("gds in Read-Only mode, please reload!")
        history.submitted <- as.character(Sys.time())
        bet <- betas(object)
        pwod.gdsn(node = bet, mul)
        history.finished <- as.character(Sys.time())
        history.command <- "Filtered with pwod (bigmelon)" 
        h <- data.frame(submitted = history.submitted,  
                        finished = history.finished,
                        command = history.command,
                        stringsAsFactors=FALSE
                        )
        for(j in colnames(h)){
            h_coln <- j
            h_index_str <- paste("history/",j,sep="")
            h_child_n <- index.gdsn(object,h_index_str,silent=TRUE)
            append.gdsn(h_child_n, val=h[,h_coln])
        }
    }
)

#dmrse <- function(betas, idmr=iDMR)
setMethod(
    f = "dmrse",
    signature(betas = "gds.class"),
    definition = function(betas, idmr = iDMR()){
        object <- betas
        betas <- betas(object)[idmr, , name = TRUE]
        dmrse(betas, idmr)
    }
) #OK

setMethod(
    f = "dmrse",
    signature(betas = "gdsn.class"),
    definition = function(betas, idmr = iDMR()){
        object <- betas
        betas <- object[idmr, ,name = TRUE]
        dmrse(betas, idmr)
    }
)

#dmrse_row <- function(betas, idmr=iDMR)
setMethod(
    f = "dmrse_row",
    signature(betas = "gds.class"),
    definition = function(betas, idmr=iDMR()){
        object <- betas
        betas <- betas(object)[idmr, ,name = TRUE]
        dmrse_row(betas, idmr)  
    }
) # AS - OK

setMethod(
    f = "dmrse_row",
    signature(betas = "gdsn.class"),
    definition = function(betas, idmr = iDMR()){
        object <- betas
        betas <- object[idmr, ,name = TRUE]
        dmrse_row( betas, idmr )
    }
)

setMethod(
    f = "dmrse_col",
    signature(betas = "gds.class"),
    definition = function(betas, idmr = iDMR()){
        object <- betas
        betas <- betas(object)[idmr, ,name = TRUE]
        dmrse_col(betas, idmr)  
    }
) # AS - OK

setMethod(
    f = "dmrse_col",
    signature(betas = "gdsn.class"),
    definition=function(betas, idmr = iDMR()){
        object <- betas
        betas <- object[idmr, ,name = TRUE]
        dmrse_col(betas, idmr)
    }
)

#seabi <- function (bn, stop, sex, X){
setMethod( # Not mem eff # Method not working either!
    f = "seabi",
    signature(bn = "gds.class"),
    definition = function(bn, stop = 1, sex=pData(bn)$sex,
                            X = fData(bn)$CHR == "X" ){
        object<- bn
        betasobj  <- betas(object)[,]
        seabi( betasobj, stop, sex, X ) 
    } 
) # AS - OK

#genki <- function(bn, g=getsnp(rownames(bn)), se=TRUE ){
setMethod(
    f = "genki",
    signature(bn = "gds.class"),
    definition = function(bn, se = TRUE){
        object <- bn
        g <- getsnp(rownames(object))
        bn <- betas(object)[g, ,name = TRUE, drop = FALSE]
        g <- rownames(bn)
        genki(bn, g, se) 
    }
) # AS - OK

setMethod(
    f= "genki",
    signature(bn="gdsn.class"),
    definition=function(bn, se = TRUE){
        object <- bn
        g <- getsnp(rownames(object))
        bn <- object[g, ,name = TRUE, drop = FALSE]
        g <- rownames(bn) 
        genki(bn, g, se)
    }
)
