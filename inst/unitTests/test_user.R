#TODO - unit tests for non-numeric accessors, & all normalisation methods

# unit test stubs, the  idea is that whenever we work on a function
# we create a unit test. any function for strictly internal use should be moved
# to test_internal.R

# TODO: dasenGds.R dbGdsn.R dfsfitGdsn.R
#       ecc.gdsn.R es2gds.R gds2mlumi.R gdsnclass_methods.R
#       GEOtoGDS.R inout.R pfilterGds.R prcompGdsn.R pwodGdsn.R
#       qnGdsn.R ranknorm.R zzz.R
#dbgdsn.R
#dfsfitGdsn.R
#ecc.gdsn.R

# es2gds.R
test_es2gds      <- function(m, file, qc = TRUE){
	data(cantaloupe)
	#t1 - basic tests to ensure function creates object and linked file
	d <- es2gds(cantaloupe,'t0.gds')
	checkTrue(file_test('-f', d$filename))
	checkTrue(inherits(d,'gds.class'))

	#e1 - error condition 1 - use existing filename
	checkException(e <- es2gds(cantaloupe,'t0.gds'))

	#e2 - error condition 2 - use non methylumi set object
	checkException(e <- es2gds(pData(cantaloupe),'t1.gds'))

	closefn.gds(d)
    unlink('t0.gds')
}

# inout.R
test_app2gds 	 <- function(m, bmln){
	data(cantaloupe)

	f <- es2gds(cantaloupe,'t2.gds')
	numsamp <- length(colnames(f))

	#t1 - append to an existing gds.class object
	f <- app2gds(cantaloupe,f)
	checkTrue(length(colnames(f)) == numsamp * 2)

	#close all opened GDS files (for test 3)
	#rv <- showfile.gds()
	#nm <- NULL; rd <- NULL
	#for (i in 1:length(rv)){
	#  names(rv[[i]]) <- c("filename", "id", "root", "readonly")
	#  class(rv[[i]]$root) <- "gdsn.class"
	#  class(rv[[i]]) <- "gds.class"
	#  nm <- c(nm, rv[[i]]$filename)
	#  rd <- c(rd, rv[[i]]$readonly)
	#}
	#for (i in 1:length(rv)){
	#  closefn.gds(rv[[i]])
	#}
	#rm(rv,nm,f)

	#t2 - append to existing file (already linked to by another gds object)
	#test doesn't work: I'm not able to recreate this scenario in test (works for user
	#though)
	#e <- es2gds(cantaloupe,'t1.gds')
	#g <- app2gds(cantaloupe,'t1.gds')
	#checkTrue(file_test('-f', g$filename))
	#checkTrue(inherits(g,'gds.class'))
	#checkTrue(length(colnames(g)) == numsamp * 2)

	#t3 - append to existing file not already linked
	#test doesn't work: I'm not able to recreate this scenario, due to the
	#close.fn function not cleaning up properly, which makes the file look as if it's
	#still open (this wont be the case for the user)
	#h <- app2gds(cantaloupe,'t2.gds')
	#checkTrue(file_test('-f', h$filename))
	#checkTrue(inherits(h,'gds.class'))
	#checkTrue(length(colnames(h)) == numsamp * 3)

	#t4 - append to new file
	i <- app2gds(cantaloupe,'t3.gds')
	checkTrue(file_test('-f', i$filename))
	checkTrue(inherits(i,'gds.class'))
	checkTrue(length(colnames(i)) == numsamp)

	#e1 - error condition 1 - try to use non gds.class object
	checkException(j <- app2gds(cantaloupe,cantaloupe))

	unlink("t1.gds", force=TRUE)
	unlink("t2.gds", force=TRUE)
	unlink("t3.gds", force=TRUE)

}

#test_iadd 		 <- function( bar, ipath=NULL, gds ){
#
#	#commented out as need to have idat files to test
#	#m <- methylumIDAT("barcode1")
#	#n <- es2gds(m,"midat.gds")
#	#o <- iadd("barcode2",n)
#	#unlink("midat.gds", force=TRUE)
#}

test_beta        <- function(object            ){
	data(cantaloupe)
	tempB <- head((betas(cantaloupe)))
	rownames(tempB) <- NULL
	colnames(tempB) <- NULL
	j <- es2gds(cantaloupe,'t4.gds')
	checkEqualsNumeric(head(betas(j)),tempB)
}
test_methylated  <- function(object            ){
	data(cantaloupe)
	tempB <- head((methylated(cantaloupe)))
	rownames(tempB) <- NULL
	colnames(tempB) <- NULL
	k <- es2gds(cantaloupe,'t5.gds')
	checkEqualsNumeric(head(methylated(k)),tempB)
}
test_unmethylated<- function(object            ){
	data(cantaloupe)
	tempB <- head((unmethylated(cantaloupe)))
	rownames(tempB) <- NULL
	colnames(tempB) <- NULL
	l <- es2gds(cantaloupe,'t6.gds')
	checkEqualsNumeric(head(unmethylated(l)),tempB)
}
test_pvals       <- function(object            ){
	data(cantaloupe)
	tempB <- head((pvals(cantaloupe)))
	rownames(tempB) <- NULL
	colnames(tempB) <- NULL
	m <- es2gds(cantaloupe,'t7.gds')
	checkEqualsNumeric(head(pvals(m)),tempB)
}
#test_fData       <- function(object            ){}
#test_pData       <- function(object            ){}
test_QCmethylated<- function(object            ){
	data(cantaloupe)
	tempB <- head(methylated(QCdata(cantaloupe)))
	rownames(tempB) <- NULL
	colnames(tempB) <- NULL
	n <- es2gds(cantaloupe,'t8.gds')
	checkEqualsNumeric(head(QCmethylated(n)),tempB)
}
test_QCunmethylated<- function(object            ){
	data(cantaloupe)
	tempB <- head(unmethylated(QCdata(cantaloupe)))
	rownames(tempB) <- NULL
	colnames(tempB) <- NULL
	o <- es2gds(cantaloupe,'t9.gds')
	checkEqualsNumeric(head(QCunmethylated(o)),tempB)
}
#test_QCrownames  <- function(object            ){}
#test_getHistory    <- function(object            ){}
#test_colnames1   <- function(x, do.NULL=TRUE, prefix=NULL){}
#test_rownames1   <- function(x, do.NULL=TRUE, prefix=NULL){}
#test_exprs		 <- function(object            ){}
#
#test_betaqn      <- function(bn                ){}
#test_naten       <- function(mn, fudge=100     ){}
#test_nanet       <- function(mn, fudge=100     ){}
#test_zot         <- function(x                 ){}
#test_nanes       <- function(mns, fudge=100    ){}
#test_danes       <- function(mn, fudge=100, ...){}
#test_danet       <- function(mn, fudge=100, ...){}
#test_daten1      <- function(mn, fudge=100, ...){}
#test_daten2      <- function(mn, fudge=100, ...){}
#test_nasen       <- function(mns, fudge=100    ){}
#test_dasen       <- function(mns, fudge=100, roco=NULL){}
#test_danen       <- function(mns, fudge=100, ...){}
#test_tost1       <- function(mn                 ){}
#test_fuks1       <- function(data               ){}
#test_swan1       <- function(mn, da=NULL        ){}
#
#test_genki1      <- function(bn, se=TRUE        ){}
#test_dmrse1      <- function(betas, idmr=iDMR() ){}
#test_dmrse_row1  <- function(betas, idmr=iDMR() ){}
#test_dmrse_col1  <- function(betas, idmr=iDMR() ){}
#test_seabi1      <- function( bn, stop=1, sex, X){}
#test_pfilter1    <- function( mn                ){}
#test_BMIQ1       <- function(x                  ){}
#test_normalizeQuantiles.gdsn.class <- function (node, selection){}
#test_dfsfit.gdsn.class.1 <- function (mn, roco  ){}
#test_dasen2      <- function(mns, fudge=100, roco=NULL){}
