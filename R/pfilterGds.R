# Arguments: See ?wateRmelon::pfilter
# NOTE: One change in how bab is computed, rather than apply beadc to matrix
# supplied we assume that a beadcount matrix is an integer matrix and remove
# values IF the supplied matrix is of integers. This is IF beadcount are
# not supplied, the betas are used instead.
pfilter.gds <- function(mn = NULL, un = NULL, bn = NULL, da = NULL, pn, bc,
                        perCount = NULL, pnthresh = NULL, perc = NULL,
                        pthresh = NULL){
    if(!is.null(list(pn, bc))){
        if(is.null(perCount)){
            perCount = 5
        }
        if(is.null(pnthresh)){
            pnthresh = 0.05
        }
        if(!is.null(pnthresh)){
            pnthresh = pnthresh
        }
        if(is.null(perc)){
            perc = 1
        }
        if(!is.null(perc)){
            perc = perc
        }
        if(is.null(pthresh)){
            pthresh = 1
        }
        if(!is.null(pthresh)){
            pthresh = pthresh
        }
        dim <- objdesp.gdsn(pn)$dim
        goodsamps <- apply.gdsn(node = pn,
                                margin = 2,
                                FUN = function(x, y, z){
                                    (sum(x)>y) < ((length(x)*z)/100)
                                },
                                as.is = "logical",
                                y = pnthresh,
                                z = perc
                                )
        bab <- apply.gdsn(node = bc,
                            margin = 1,
                            as.is = "integer",
                            FUN = function(x, y){
                            if((class(x)=="integer")) x[x<3] <- NA
                            # Incase betas are supplied instead of NBeads
                            length(which(is.na(x[y])=="TRUE"))
                            },
                            y = goodsamps
                            )
        badbead_log <- bab > ((dim[2] * perCount)/100)
        badbead <- which(badbead_log)
        bap <- apply.gdsn(node = pn,
                            margin = 1,
                            as.is = "logical",
                            FUN = function(x, y, z, p){
                                sum(x[y] > z) > ((length(x) * p/100))
                            },
                            y = goodsamps,
                            z = pnthresh,
                            p = pthresh
                            )
        badp <- which(bap)
        message(sum(!goodsamps), " samples having", perc,
            "% of sites with a detection p-value greater than ",
            pnthresh, " were removed.")
        message("Samples removed: ", colnames(bc)[!goodsamps])
        message(length(badbead), " sites were remove as beadcount <3 in ",
            perCount, "% of samples.")
        message(length(badp), " sites having ", pthresh,
            "% of samples with a detection p-value greater than ",
            pnthresh, " were removed.")
        # Only does logical return.
        return(list(probes = (!bap & !badbead_log), samples = goodsamps))
    }
}
