db.gdsn <- function(gds, mns, uns){ # {{{  
    dim <- objdesp.gdsn(mns)$dim
    stopifnot(dim == objdesp.gdsn(uns)$dim)
    # Creating new node for meth + unmeth combination
    n.t <- add.gdsn(node = gds, name = "methunmeth", storage = "float64",
                    valdim = c(dim[1], 0), val = NULL, replace = TRUE)
    # Store meth and unmeth into the same node
    for(x in 1:dim[2]){
        mnval <- readex.gdsn(mns, sel = list(NULL, x))
        append.gdsn(n.t, mnval)
        unval <- readex.gdsn(uns, sel = list(NULL, x))
        append.gdsn(n.t, unval)
    }
    # Normalization
    qn.gdsn(gds = gds, target = index.gdsn(gds, "methunmeth"), newnode = "db")
    # Initializing new nodes for new normalized values
    db.m <- add.gdsn(gds, "db.meth" , storage = "float64", 
                    valdim = c(dim[1], 0), val = NULL, replace = TRUE)
    db.u <- add.gdsn(gds, "db.unmeth" , storage = "float64", 
                    valdim = c(dim[1], 0), val = NULL, replace = TRUE)
    # Split normalized values into new nodes.
    for(x in seq(1,(2*dim[2])-1, 2)){ 
        val <- readex.gdsn(index.gdsn(gds, path = "db"),
                            sel = list(NULL, c(x, x + 1)))
        append.gdsn(db.m, val[ ,1])
        append.gdsn(db.u, val[ ,2])
    }
} # }}}
