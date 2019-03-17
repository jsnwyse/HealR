healpix.nested.write.patch <- function( x = NULL, file = NULL, rows = NULL, cols = NULL, patch = 0, colnum = 1, verbose = TRUE )
{
  if( is.null(x) ) stop("please pass argument x")
  if( class(x) != "HealR" ) stop("argument is not of class 'HealR'") 
  #file = NULL, nside = NULL, rows = NULL, cols = NULL, patch = 0, colnum=1, plot = TRUE, mask = NULL, main="", output.Matrix=FALSE

	nfile <- length( file )
	status <- rep( 0, nfile )
  nside <- x$nside
  
  if( !is.null(rows) | is.null(cols) )
  {
    # get rows and sort
    if( !is.null(rows) ) rows <- sort( rows ) else rows <- 1:nside
    # get cols and sort
    if( !is.null(cols) ) cols <- sort( cols ) else cols <- 1:nside
    nrows <- length( rows )
    ncols <- length( cols )
  }else{
    rows <- 1:nside
    cols <- 1:nside
    nrows <- nside
    ncols <- nside
  }
  
  # need to define the row.idx/col.idx for all in the patch-- is this less efficient (appears best way to do it)
  col.idx <- rep.int( cols, times=rep(nrows,ncols) ) - 1
  row.idx <- rep( rows, ncols ) - 1
  
  nwrite <- nrows * ncols
  
  if( nfile > 0 & length(x$X[[1]]) < nwrite ) 
      stop("length of vector to write to file does not match the\n dimensions passed: check the rows/cols arguments")
	
	if( verbose ) cat("\n ... Searching the current working directory for existing files ... ")
	Files <- list.files()
	
	for( k in 1:nfile )
	{
		if( !any( Files == file[k] ) )
		{
				healpix.nested.create.file( file[k], nside, nfields = 2 )
		    if( verbose ) cat( paste0("\n\t ... file '", file[k] ,"' didn't already exist- it has been created ..."))
		}
	}
	if( verbose ) cat("\n ... Files are now ready for writing ... \n")
	
	
	for( k in 1:nfile )
	{
	  status <- rep(0,2)
		w <- 	.C( "astro_write_nested_healpix_patch", 
					as.character(file[k]), 
					as.double(x$X[[k]]), 
					as.integer(nside), 
					as.integer(nwrite), 
					as.integer(row.idx), #row idx
					as.integer(col.idx), #col idx
					as.integer(colnum), 
					status = as.integer(status), 
					as.integer(patch),
					PACKAGE = "HealR" )
		if( w$status[1] > 0 ) stop( fitsio.match.error.code( w$status[1], file[k] ) )
		status[k] <- w$status[1]
	}
	

	out <- list( file = file, status = status )
	
	return( out )
}
