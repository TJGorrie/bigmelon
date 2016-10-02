.onLoad <- function(libname = find.package("bigmelon"), pkgname="bigmelon"){
    options(runLast = TRUE)
    .GlobalEnv$.Last <- function(x = ls(.GlobalEnv)){
        for(i in x){
            if(inherits(get(i), "gds.class")){
            a <- try(closefn.gds(get(i)), silent = TRUE)
            }
        }
    }
}


