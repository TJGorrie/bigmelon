dasen.gds <-   function(gds,
                        node,
                        mns,
                        uns,
                        onetwo,
                        roco,
                        fudge,
                        ret2
                        ){ # {{{  
    # Assuming that mns and uns are 1 element character strings!not gdsn.nodes
    if(length(mns) == 1) mns <- index.gdsn(gds, mns)
    if(length(uns) == 1) uns <- index.gdsn(gds, uns)
    if(length(onetwo) == 1)  onetwo <- read.gdsn(index.gdsn(gds, onetwo))
    if(class(onetwo) == 'gdsn.class') onetwo <- read.gdsn(onetwo)
    f <- createfn.gds("temp.gds", allow.duplicate = TRUE)
    dim <- objdesp.gdsn(mns)$dim
    # NORMALIZING
    dfsfit.gdsn(f, targetnode = mns, roco = roco, newnode = "mnsc",
                onetwo = onetwo)
    dfsfit.gdsn(f, targetnode = uns, roco = NULL, newnode = "unsc",
                onetwo = onetwo)
    ## Splitting arrays by probe Type
    # Initiliazing new nodes
    mI  <- add.gdsn(f, "metI" , storage = "float64", 
                valdim = c(sum(onetwo=='I'), 0),val = NULL, replace = TRUE)
    mII <- add.gdsn(f, "metII" , storage = "float64", 
                valdim = c(sum(onetwo=='II'), 0), val = NULL, replace = TRUE)
    uI  <- add.gdsn(f, "umeI" , storage = "float64", 
                valdim = c(sum(onetwo=='I'), 0), val = NULL, replace = TRUE)
    uII <- add.gdsn(f, "umeII" , storage = "float64", 
                valdim = c(sum(onetwo=='II'), 0), val = NULL, replace = TRUE)
    # Separating probes into type I and II and appending to object col by col.
    for(x in 1:dim[2]){
        append.gdsn(mI , readex.gdsn(index.gdsn(f, "mnsc"),
                    sel = list(onetwo=='I' , x)))
        append.gdsn(mII, readex.gdsn(index.gdsn(f, "mnsc"),
                    sel = list(onetwo=='II', x)))
        append.gdsn(uI , readex.gdsn(index.gdsn(f, "unsc"),
                    sel = list(onetwo=='I' , x)))
        append.gdsn(uII, readex.gdsn(index.gdsn(f, "unsc"),
                    sel = list(onetwo=='II', x)))
    }
    # Normalize seperately using qn.gdsn
    qn.gdsn(f, target = mI , newnode = "metIqn" )
    qn.gdsn(f, target = mII, newnode = "metIIqn")
    qn.gdsn(f, target = uI , newnode = "umeIqn" )
    qn.gdsn(f, target = uII, newnode = "umeIIqn")
    # Relcalculating betas 
    # Creating new node where normalized betas will be stored. / replacement
    n.t <- add.gdsn(gds, name = node, storage = "float64", 
                    valdim = c(dim[1], 0), val = NULL, replace = TRUE)
    if(ret2 == TRUE){
        n.m <- add.gdsn(gds, "methylated", storage = "float64",
                        valdim = c(dim[1], 0), val = NULL, replace = TRUE)
        n.u <- add.gdsn(gds, "unmethylated", storage = "float64",
                        valdim = c(dim[1], 0), val = NULL, replace = TRUE)
    }
    for(x in 1:dim[2]){
        meth   <- rep(x = NA, times = dim[1])
        unmeth <- rep(x = NA, times = dim[1])
        meth[onetwo == 'I']    <- readex.gdsn(index.gdsn(f,   "metIqn"),
                                            sel = list(NULL , x))
        meth[onetwo == 'II']   <- readex.gdsn(index.gdsn(f,  "metIIqn"),
                                            sel = list(NULL , x))
        unmeth[onetwo == 'I']  <- readex.gdsn(index.gdsn(f,   "umeIqn"),
                                            sel = list(NULL , x))
        unmeth[onetwo == 'II'] <- readex.gdsn(index.gdsn(f,  "umeIIqn"),
                                            sel = list(NULL , x))
        beta <- meth/(meth + unmeth + fudge)
        append.gdsn(n.t, beta)
        if(ret2 == TRUE){
            append.gdsn(n.m, meth)
            append.gdsn(n.u, unmeth)
        }
    }
    closefn.gds(f)
    unlink("temp.gds", force = TRUE)
} # }}}

naten.gds <-  function( gds,
                        node,
                        mns,
                        uns,
                        fudge,
                        ret2
                        ){ # {{{
    f <- createfn.gds("temp.gds", allow.duplicate = TRUE)
    dim <- objdesp.gdsn(mns)$dim
    ## NORMALIZING
    qn.gdsn(f, target = mns, newnode = "natenmeth")
    qn.gdsn(f, target = uns, newnode = "natenunmeth")
    ## Recalculating Betas  
    # Creating new node for betas.
    n.t <- add.gdsn(gds, name = node , storage = "float64", 
                    valdim = c(dim[1], 0), val = NULL, replace = TRUE)
    if(ret2 == TRUE){
        n.m <- add.gdsn(gds, "methylated", storage = "float64",
                        valdim = c(dim[1], 0), val = NULL, replace = TRUE)
        n.u <- add.gdsn(gds, "unmethylated", storage = "float64",
                        valdim = c(dim[1], 0), val = NULL, replace = TRUE)
    }
    for(x in 1:dim[2]){
        mn <- readex.gdsn(index.gdsn(f , "natenmeth"), sel = list(NULL, x))
        un <- readex.gdsn(index.gdsn(f,"natenunmeth"), sel = list(NULL, x))
        beta <- mn / ( mn + un + fudge )
        append.gdsn(n.t, beta)
        if(ret2 == TRUE){
            append.gdsn(n.m, mn)
            append.gdsn(n.u, un)
        }
    }  
    closefn.gds(f)
    unlink("temp.gds", force = TRUE) 
} # }}}

danen.gds <-  function( gds,
                        node,
                        mns,
                        uns,
                        onetwo,
                        roco,
                        fudge,
                        ret2
                        ){ # {{{  
    if(length(onetwo) == 1)  onetwo <- read.gdsn(index.gdsn(gds, onetwo))
    if(class(onetwo) == 'gdsn.class') onetwo <- read.gdsn(onetwo)
    f <- createfn.gds("temp.gds", allow.duplicate = TRUE)
    dim <- objdesp.gdsn(mns)$dim
    ## NORMALIZING
    dfsfit.gdsn(f, targetnode = mns, roco = roco,
                newnode = "mnsc", onetwo = onetwo)
    dfsfit.gdsn(f, targetnode = uns, roco = NULL,
                newnode = "unsc", onetwo = onetwo)
    ## Recalculating Betas
    # Creating new node for betas
    n.t <- add.gdsn(gds, name = node , storage = "float64", 
                    valdim = c(dim[1], 0), val = NULL, replace = TRUE)
    if(ret2 == TRUE){
        n.m <- add.gdsn(gds, "methylated", storage = "float64",
                        valdim = c(dim[1], 0), val = NULL, replace = TRUE)
        n.u <- add.gdsn(gds, "unmethylated", storage = "float64",
                        valdim = c(dim[1], 0), val = NULL, replace = TRUE)
    }
    for(x in 1:dim[2]){
        mn <- readex.gdsn(index.gdsn(f, "mnsc"), sel = list(NULL, x))
        un <- readex.gdsn(index.gdsn(f, "unsc"), sel = list(NULL, x))
        beta <- mn / ( mn + un + fudge )
        append.gdsn(n.t, beta)
        if(ret2 == TRUE){
            append.gdsn(n.m, mn)
            append.gdsn(n.u, un)
        }
    }
    closefn.gds(f)
    unlink("temp.gds", force = TRUE) 
} # }}}

daten1.gds <-  function(gds,
                        node,
                        mns,
                        uns,
                        onetwo,
                        roco,
                        fudge,
                        ret2
                        ){ # {{{
    if(length(onetwo) == 1)  onetwo <- read.gdsn(index.gdsn(gds, onetwo))
    if(class(onetwo) == 'gdsn.class') onetwo <- read.gdsn(onetwo)
    f <- createfn.gds("temp.gds", allow.duplicate = TRUE)
    dim <- objdesp.gdsn(mns)$dim
    ## NORMALIZING
    dfsfit.gdsn(f,  targetnode = mns, roco = roco, 
                newnode = "mnsc", onetwo = onetwo)
    dfsfit.gdsn(f,  targetnode = uns, roco = NULL,
                newnode = "unsc", onetwo = onetwo)
    qn.gdsn(f, target = index.gdsn(f,"mnsc"), newnode = "daten1meth")
    qn.gdsn(f, target = index.gdsn(f,"unsc"), newnode = "daten1unmeth") 
    ## Recalculating Betas
    # Creating new node for betas
    n.t <- add.gdsn(gds, name = node , storage = "float64", 
                    valdim = c(dim[1], 0), val = NULL, replace = TRUE)
    if(ret2 == TRUE){
        n.m <- add.gdsn(gds, "methylated", storage = "float64",
                        valdim = c(dim[1], 0), val = NULL, replace = TRUE)
        n.u <- add.gdsn(gds, "unmethylated", storage = "float64",
                        valdim = c(dim[1], 0), val = NULL, replace = TRUE)
    }
    for(x in 1:dim[2]){
        mn <- readex.gdsn(index.gdsn(f, "daten1meth"), sel = list(NULL, x))
        un <- readex.gdsn(index.gdsn(f, "daten1unmeth"), sel = list(NULL, x))
        beta <- mn / ( mn + un + fudge )
        append.gdsn(n.t, beta)
        if(ret2 == TRUE){
            append.gdsn(n.m, mn)
            append.gdsn(n.u, un)
        }
    }
    closefn.gds(f)
    unlink("temp.gds", force = TRUE)
} # }}}

daten2.gds <-  function(gds,
                        node,
                        mns,
                        uns,
                        onetwo,
                        roco,
                        fudge,
                        ret2
                        ){ # {{{
    if(length(onetwo) == 1)  onetwo <- read.gdsn(index.gdsn(gds, onetwo))
    if(class(onetwo) == 'gdsn.class') onetwo <- read.gdsn(onetwo) 
    f <- createfn.gds("temp.gds", allow.duplicate = TRUE)
    dim <- objdesp.gdsn(mns)$dim
    ## NORMALIZATION
    dfsfit.gdsn(f,  targetnode = mns, roco = roco,
                newnode = "mnsc", onetwo = onetwo)
    dfsfit.gdsn(f,  targetnode = uns, roco = roco,
                newnode = "unsc", onetwo = onetwo)
    qn.gdsn(f, target = index.gdsn(f,"mnsc"), newnode = "daten2meth")
    qn.gdsn(f, target = index.gdsn(f,"unsc"), newnode = "daten2unmeth")
    ## Recalculating Betas
    # Adding new node for betas
    n.t <- add.gdsn(gds, name = node, storage = "float64", 
                    valdim = c(dim[1], 0), val = NULL, replace = TRUE)
    if(ret2 == TRUE){
        n.m <- add.gdsn(gds, "methylated", storage = "float64",
                        valdim = c(dim[1], 0), val = NULL, replace = TRUE)
        n.u <- add.gdsn(gds, "unmethylated", storage = "float64",
                        valdim = c(dim[1], 0), val = NULL, replace = TRUE)
    }
    for(x in 1:dim[2]){ 
        mn <- readex.gdsn(index.gdsn(f, "daten2meth"), sel = list(NULL, x))
        un <- readex.gdsn(index.gdsn(f, "daten2unmeth"), sel = list(NULL, x))
        beta <- mn / ( mn + un + fudge )
        append.gdsn(n.t, beta)
        if(ret2 == TRUE){
            append.gdsn(n.m, mn)
            append.gdsn(n.u, un)
        }
    }
    closefn.gds(f)
    unlink("temp.gds", force = TRUE)
} # }}}

nasen.gds <-  function(gds ,
                        node,
                        mns,
                        uns,
                        onetwo,
                        fudge,
                        ret2
                        ){ # {{{
    if(length(onetwo) == 1)  onetwo <- read.gdsn(index.gdsn(gds, onetwo))
    if(class(onetwo) == 'gdsn.class') onetwo <- read.gdsn(onetwo)
    f <- createfn.gds("temp.gds", allow.duplicate = TRUE)
    dim <- objdesp.gdsn(mns)$dim
    ## Normalization
    ## Splitting arrays by probe Type
    # Initiliazing new nodes
    mI  <- add.gdsn(f, "metI" , storage = "float64", 
        valdim = c(sum(onetwo=='I'), 0), val = NULL, replace = TRUE)
    mII <- add.gdsn(f, "metII" , storage = "float64", 
        valdim = c(sum(onetwo=='II'), 0), val = NULL, replace = TRUE)
    uI  <- add.gdsn(f, "umeI" , storage = "float64", 
        valdim = c(sum(onetwo=='I'), 0), val = NULL, replace = TRUE)
    uII <- add.gdsn(f, "umeII" , storage = "float64", 
        valdim = c(sum(onetwo=='II'), 0), val = NULL, replace = TRUE)
    # Separating probes
    for(x in 1:dim[2]){
        append.gdsn(mI , readex.gdsn(mns, sel = list(onetwo == 'I' , x)))
        append.gdsn(mII, readex.gdsn(mns, sel = list(onetwo == 'II', x)))
        append.gdsn(uI , readex.gdsn(uns, sel = list(onetwo == 'I' , x)))
        append.gdsn(uII, readex.gdsn(uns, sel = list(onetwo == 'II', x)))
    }
    # Norm
    qn.gdsn(f, target = mI , newnode = "metIqn" )
    qn.gdsn(f, target = mII, newnode = "metIIqn")
    qn.gdsn(f, target = uI , newnode = "umeIqn" )
    qn.gdsn(f, target = uII, newnode = "umeIIqn")
    ## Recalculating Betas
    # Adding new node
    n.t <- add.gdsn(gds, name = node, storage = "float64", 
                    valdim = c(dim[1], 0), val = NULL, replace = TRUE)
    if(ret2 == TRUE){
        n.m <- add.gdsn(gds, "methylated", storage = "float64",
                        valdim = c(dim[1], 0), val = NULL, replace = TRUE)
        n.u <- add.gdsn(gds, "unmethylated", storage = "float64",
                        valdim = c(dim[1], 0), val = NULL, replace = TRUE)
    }
    for(x in 1:dim[2]){
        meth   <- rep(x = NA, times = dim[1])
        unmeth <- rep(x = NA, times = dim[1])
        meth[onetwo == 'I']    <- readex.gdsn(index.gdsn(f,   "metIqn"),
                                            sel = list(NULL , x))
        meth[onetwo == 'II']   <- readex.gdsn(index.gdsn(f,  "metIIqn"),
                                            sel = list(NULL , x))
        unmeth[onetwo == 'I']  <- readex.gdsn(index.gdsn(f,   "umeIqn"),
                                            sel = list(NULL , x))
        unmeth[onetwo == 'II'] <- readex.gdsn(index.gdsn(f,  "umeIIqn"),
                                            sel = list(NULL , x))
        beta <- meth/(meth + unmeth + fudge)
        append.gdsn(n.t, beta)
        if(ret2 == TRUE){
            append.gdsn(n.m, meth)
            append.gdsn(n.u, unmeth)
        }
    }
    closefn.gds(f)
    unlink("temp.gds", force = TRUE)
} # }}}

nanet.gds <-  function( gds,
                        node,
                        mns,
                        uns,
                        fudge,
                        ret2
                        ){ # {{{
    f <- createfn.gds("temp.gds", allow.duplicate = TRUE)
    dim <- objdesp.gdsn(mns)$dim
    ## Normalization
    db.gdsn(f, mns, uns)
    ## Recalculating Betas
    # Adding new node for betas
    n.t <- add.gdsn(gds, name = node , storage = "float64", 
                    valdim = c(dim[1], 0), val = NULL, replace = TRUE) 
    if(ret2 == TRUE){
        n.m <- add.gdsn(gds, "methylated", storage = "float64",
                        valdim = c(dim[1], 0), val = NULL, replace = TRUE)
        n.u <- add.gdsn(gds, "unmethylated", storage = "float64",
                        valdim = c(dim[1], 0), val = NULL, replace = TRUE)
    }
    for(x in 1:dim[2]){
        mn <- readex.gdsn(index.gdsn(f, "db.meth"), sel = list(NULL, x))
        un <- readex.gdsn(index.gdsn(f, "db.unmeth"), sel = list(NULL, x))
        beta <- mn / ( mn + un + fudge )
        append.gdsn(n.t, beta)
        if(ret2 == TRUE){
            append.gdsn(n.m, mn)
            append.gdsn(n.u, un)
        }
    }
    closefn.gds(f)
    unlink("temp.gds", force = TRUE)
} # }}}

danet.gds <- function(  gds,
                        node,
                        mns,
                        uns,
                        onetwo,
                        roco,
                        fudge,
                        ret2
                        ){ # {{{
    if(length(onetwo) == 1) onetwo <- read.gdsn(index.gdsn(gds, onetwo))
    if(class(onetwo) == 'gdsn.class') onetwo <- read.gdsn(onetwo)
    f <- createfn.gds("temp.gds", allow.duplicate = TRUE)
    dim <- objdesp.gdsn(mns)$dim
    ## Normalization
    dfsfit.gdsn(f, targetnode = mns, roco = roco,
                newnode = "mnsc", onetwo = onetwo)
    dfsfit.gdsn(f, targetnode = uns, roco = roco,
                newnode = "unsc", onetwo = onetwo)
    db.gdsn(f, mns = index.gdsn(f, "mnsc"), uns = index.gdsn(f, "unsc"))
    ## Recalculating Betas
    # Adding new node for betas
    n.t <- add.gdsn(gds, name = node, storage = "float64", 
                    valdim = c(dim[1], 0), val = NULL, replace = TRUE)
    if(ret2 == TRUE){
        n.m <- add.gdsn(gds, "methylated", storage = "float64",
                        valdim = c(dim[1], 0), val = NULL, replace = TRUE)
        n.u <- add.gdsn(gds, "unmethylated", storage = "float64",
                        valdim = c(dim[1], 0), val = NULL, replace = TRUE)
    }
    for(x in 1:dim[2]){
        mn <- readex.gdsn(index.gdsn(gds, "db.meth"), sel = list(NULL, x))
        un <- readex.gdsn(index.gdsn(gds, "db.unmeth"), sel = list(NULL, x))
        beta <- mn / ( mn + un + fudge )
        append.gdsn(n.t, beta)
        if(ret2 == TRUE){
            append.gdsn(n.m, mn)
            append.gdsn(n.u, un)
        }
    }
    closefn.gds(f)
    unlink("temp.gds", force = TRUE)
} # }}}

nanes.gds <- function(  gds,
                        node,
                        mns,
                        uns,
                        onetwo,
                        fudge,
                        ret2
                        ){ # {{{
    if(length(onetwo) == 1) onetwo <- read.gdsn(index.gdsn(gds, onetwo))
    if(class(onetwo) == 'gdsn.class') onetwo <- read.gdsn(onetwo)
    f <- createfn.gds("temp.gds", allow.duplicate = TRUE)
    met <- mns
    ume <- uns
    dim <- objdesp.gdsn(met)$dim
    ## Splitting arrays by probe Type
    # Initiliazing new nodes
    mI  <- add.gdsn(f, "metI" , storage = "float64", 
                valdim = c(sum(onetwo == 'I'), 0), val = NULL, replace = TRUE)
    mII <- add.gdsn(f, "metII" , storage = "float64", 
                valdim = c(sum(onetwo == 'II'), 0), val = NULL,replace = TRUE)
    uI  <- add.gdsn(f, "umeI" , storage = "float64", 
                valdim = c(sum(onetwo == 'I'), 0), val = NULL, replace = TRUE)
    uII <- add.gdsn(f, "umeII" , storage = "float64", 
                valdim = c(sum(onetwo == 'II'), 0), val = NULL,replace = TRUE)
    # Separating probes
    for(x in 1:dim[2]){
        append.gdsn(mI , readex.gdsn(met, sel = list(onetwo == 'I' , x)))
        append.gdsn(mII, readex.gdsn(met, sel = list(onetwo == 'II', x)))
        append.gdsn(uI , readex.gdsn(ume, sel = list(onetwo == 'I' , x)))
        append.gdsn(uII, readex.gdsn(ume, sel = list(onetwo == 'II', x)))
    }
    ## Normalization per nanes method
    qn.gdsn(f, target = index.gdsn(f,"metI"), newnode = "metIqn")
    qn.gdsn(f, target = index.gdsn(f,"umeI"), newnode = "umeIqn")
    db.gdsn(f, mns = index.gdsn(f,"metII"), uns = index.gdsn(f, "umeII"))
    ## Recalculating Betas
    # Initialising nodes for betas, meth and unmeth
    n.t <- add.gdsn(gds, name = node , storage = "float64", 
                    valdim = c(dim[1], 0), val = NULL, replace = TRUE)
    if(ret2 == TRUE){
        n.m <- add.gdsn(gds, "methylated" , storage = "float64", 
                        valdim = c(dim[1], 0), val = NULL, replace = TRUE)
        n.u <- add.gdsn(gds, "unmethylated" , storage = "float64", 
                        valdim = c(dim[1], 0), val = NULL, replace = TRUE)
    }
    for(x in 1:dim[2]){
        meth <- rep(x = NA, times = dim[1])
        unmeth <- rep(x = NA, times = dim[1])
        meth[onetwo == 'I']    <- readex.gdsn(index.gdsn(f,   "metIqn"),
                                                sel = list(NULL , x))
        meth[onetwo == 'II']   <- readex.gdsn(index.gdsn(f,  "db.meth"),
                                                sel = list(NULL , x))
        unmeth[onetwo == 'I']  <- readex.gdsn(index.gdsn(f,   "umeIqn"),
                                                sel = list(NULL , x))
        unmeth[onetwo == 'II'] <- readex.gdsn(index.gdsn(f,"db.unmeth"),
                                                sel = list(NULL , x))
        beta <- meth/(meth + unmeth + fudge)
        append.gdsn(n.t, beta)
        if(ret2 == TRUE){
            append.gdsn(n.m, meth)
            append.gdsn(n.u, unmeth)
        }
    } 
    closefn.gds(f)
    unlink("temp.gds", force = TRUE)
} # }}}

danes.gds <-  function( gds,
                        node,
                        mns,
                        uns,
                        onetwo,
                        fudge,
                        ret2,
                        ...
                        ){ # {{{
    if(length(onetwo) == 1)  onetwo <- read.gdsn(index.gdsn(gds, onetwo))
    if(class(onetwo) == 'gdsn.class') onetwo <- read.gdsn(onetwo)
    f.f <- createfn.gds("tempdanes.gds", allow.duplicate = TRUE)
    ## Normalization
    dfsfit.gdsn(gds = f.f, targetnode = mns, newnode = "mnsc",
                onetwo = onetwo, ...) # roco = NULL?
    dfsfit.gdsn(gds = f.f, targetnode = uns, newnode = "unsc",
                onetwo = onetwo, ...) # roco = NULL?
    nanes.gds(  gds = gds,
                node = node,
                mns = index.gdsn(f.f,"mnsc"),
                uns = index.gdsn(f.f,"unsc"), 
                onetwo = onetwo,
                fudge = fudge,
                ret2 = ret2
                )
    closefn.gds(f.f)
    unlink("tempdanes.gds", force = TRUE)
} # }}}
