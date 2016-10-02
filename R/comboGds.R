# combo.gds
# Difficult, uncompromising and insulting. Like me in the morning. #TGS
# Forcefully coerces secondary gds node to mirror primarys dimensions forcing 
# NAs(if necessary).
# TODO: Make coherent pData combination.
combo.gds <- function(file, primary, secondary){
    # file = newgds file
    # Primary = gds.class object
    # Secondary = gds.class object
    # Assumes both objects are bigmelon made gds.classes!
    # Preliminary Safety Checks
    # Are they gds objects
    if(!class(primary) == "gds.class") stop("Primary is not a gds object!")
    if(!class(secondary) == "gds.class") stop("Secondary is not a gds object!")
    history.submitted <- as.character(Sys.time())  
    if(!(file_test('-f', file))){
        message(paste(file, "doesn't exist, creating new file"))
        bmln <- newgds(file)
    } else if((file_test('-f', file))){
        stop('File already exists')
    }
    primnod <- ls.gdsn(primary)
    seconod <- ls.gdsn(secondary)
    # Identify Common Nodes
    sharenod <- primnod[primnod %in% seconod]
    names(sharenod) <- sharenod
    # forces row length of primary gds object
    forcerow <- rownames(primary)
    primdim <- objdesp.gdsn(index.gdsn(primary, "betas"))$dim
    secodim <- objdesp.gdsn(index.gdsn(secondary, "betas"))$dim
    for(x in c("betas", "pvals", "methylated", "unmethylated", "NBeads")){
        if(x %in% sharenod){
            # Determines trait class of sharenodes and uses same storage
            if(getTrait(index.gdsn(primary, x)) == 
                getTrait(index.gdsn(secondary, x))
                ){
                trait <- tolower(getTrait(index.gdsn(primary, x))) 
            } else { 
                trait <- "float64"
            }
            n.t <- add.gdsn(bmln, x, valdim = c(primdim[1], 0),
                            val = NULL, storage = trait, replace = TRUE)
            message(paste("Appending", x, "from primary node to", file))
            for(m in 1:primdim[2]){ # This will coerce NA's for missing probes
                val <- index.gdsn(primary, x)[forcerow, m, name = FALSE]
                append.gdsn(n.t, val)
            }
            message(paste("Appending", x, "from secondary node to:", file))
            for(k in 1:secodim[2]){
                val <- index.gdsn(secondary, x)[forcerow, k, name = FALSE]
                append.gdsn(n.t, val)
            }  
        }
    }
    if(!"fData" %in% sharenod){
        # If fData is missing from one node, do not append fData(?). Assumption
        # is if fData is missing, this operation cannot be completed anyway.
    } else {
        fDPrim <- fData(primary)[forcerow, ]
        message("Creating fData node...")
        add.gdsn(bmln, "fData", val = fDPrim)
        message("Appending fData...")
    }
    pD <- data.frame(t(rbind(c(colnames(primary), colnames(secondary)))))
    colnames(pD) <- "barcode"
    #Current pData only combines column names.
    message("Creating pData node...")
    add.gdsn(bmln, "pData", val = pD)
    message("Appending pData...")
    for(x in c("QCmethylated", "QCunmethylated")){
        if(x %in% sharenod){
            message(paste("Creating", x, "node..."))
            n.t <- add.gdsn(bmln, x, valdim = c(length(QCrownames(primary)),0),
                            val = NULL, replace = TRUE, storage="int32")
            message(paste("Appending Primary Node to:", file))
            for(m in 1:primdim[2]){
                val <- readex.gdsn(index.gdsn(primary, x), sel = list(NULL, m))
                append.gdsn(n.t, val, check = FALSE)
            }
            message(paste("Appending Secondary Node to:", file))
            for(k in 1:secodim[2]){
                val <- readex.gdsn(index.gdsn(secondary, x), sel=list(NULL, k))
                append.gdsn(n.t, val, check = FALSE)
            }  
        }
    }
    add.gdsn(bmln, "QCrownames", val = QCrownames(primary))
    history.finished <- as.character(Sys.time())
    history.command <- paste0("Combo of: ",
                                primary[[1]], " and ", secondary[[1]], ".")
    h <- data.frame(submitted = history.submitted, 
                    finished = history.finished, 
                    command = history.command,   
                    stringsAsFactors = FALSE)
    add.gdsn(bmln, "history", val = h)
    add.gdsn(bmln, "paths", val = c("fData/Probe_ID", "pData/barcode"))
    return(bmln)
}

getTrait <- function(x){
    objdesp.gdsn(x)$trait
}
