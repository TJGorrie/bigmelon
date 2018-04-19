getpheno <- function(geo){
    a <- GEOquery::getGEO(geo, destdir = paste0("./",geo), GSEMatrix=T,getGPL=F )
    pd <- phenoData(a[[1]])@data
    write.csv(pd, file=paste0('./',geo,'/',geo,'pheno.csv'))
    message('Removing unpacked files!')
    leftovers <- dir(paste0("./",geo), recursive = TRUE)[grepl(".gz", dir(paste0("./",geo), recursive = TRUE))]
    file.remove(paste0("./",geo,'/',leftovers))
    return(pd)
}

geotogds <- function(geo, gds=NULL, method = "wget", keepidat = F, keeptar = F, ...){
    if(is.null(gds)) gds <- paste0(geo, '.gds')
    if(!(grepl(x = geo, pattern='.tar.gz'))){
        url <- paste0("https://www.ncbi.nlm.nih.gov/geo/download/?acc=", geo, "&format=file")
        tarfile <- paste0("./", geo,".tar.gz")
        message('Downloading \"', geo, '\" to \"', tarfile,'\".')
        download.file(url, destfile = tarfile, method = method, quiet = FALSE)
        dest <- paste0("./", geo)
    } else {
        tarfile <- geo
        geo <- strsplit(x=geo, split=".tar.gz")[[1]][1]
        dest <- paste0("./", geo)
    }
    message('Unpacking file to \'', geo,'/\'')
    untar(tarfile, exdir = dest, compressed = 'gzip')
    message('Enumerating .idats!')
    barsgz <- dir(dest, recursive = TRUE)[grepl(".idat.gz", dir(dest, recursive = TRUE))]
    message('Unpacking ', length(barsgz), " .idats!")
    for(z in barsgz){
        GEOquery::gunzip(filename = paste0(dest, "/", z), destname = paste0(dest, "/", gsub(x = z, pattern = ".gz", replacement='')), remove = TRUE)
    }
    # if(tidy){
    message('Removing unpacked files!')
    leftovers <- dir(dest, recursive = TRUE)[grepl(".gz", dir(dest, recursive = TRUE))]
    file.remove(paste0(dest,'/',leftovers))
    # }

    pda <- getpheno(geo)
    message('Reading Data!')
    gds <- iadd2(path = dest, gds = gds, ...)
    newpd <- try(cbind(pData(gds),pda), silent=T)
    if(!inherits(newpd, 'try-error')){
    add.gdsn(gds, name="pData", valdim=dim(newpd),
                    val=data.frame(lapply(newpd, as.character),
                            stringsAsFactors = FALSE), replace=TRUE)
                        }
    gc() # Release Mem
    if(!keepidat){
    message('Removing idats folder')
    unlink(dest, recursive = T, force = T)
    }
    if(!keeptar){
    message('Removing tar')
    file.remove(tarfile)
    }
    return(gds)
}
