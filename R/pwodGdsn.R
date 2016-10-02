pwod.gdsn <- function(node, mul = 4){ # {{{
    # create temporary file in working directory
    f <- createfn.gds("temp.gds", allow.duplicate = TRUE) 
    dim <- objdesp.gdsn(node)$dim
    # Create new node in "temp.gds"
    n.t <- add.gdsn(node = f, name = "pwod", storage = "float64", 
                    valdim = c(dim[1], 0), val = NULL, replace = TRUE)
    # pwod.R outputs to created node.
    apply.gdsn(node,
                margin = 1,
                FUN = function(x,y){
                    # pwod.R compute row by row
                    quan <- fivenum(x)
                    iqr <- quan[4] - quan[2]
                    bounds <- c(quan[4] + (iqr * y),
                    quan[2] - (iqr * y)) 
                    # Upper Bound is [1], Lower Bound is [2]
                    d <- x > bounds[1] | x < bounds[2]
                    x[d] <- NA 
                    x
                    },
                y = mul,
                as.is = "gdsnode",
                target.node = n.t
                )
    copyto.gdsn(node = getfolder.gdsn(node), source = index.gdsn(f, "pwod"),
                name = paste0("pwod_", objdesp.gdsn(node)$name))
    # Close + Delete temp file
    closefn.gds(f)
    unlink("temp.gds")
} # }}}
