healpix.nested.write.patch <- function( X = NULL, file = NULL, patch = 0, colnum = 1 )
{

	nfile <- length( file )
	status <- numeric( nfile )
	
	if( nfile > 1 )
	{
		nside <- nrow( X[[1]] )
		Y <- matrix( nrow = nfile, ncol = nside*nside )
		for( k in 1:nfile ) Y[ k , ] <- as.vector( X[[k]] )
	}else{
		nside <- nrow( X )
		Y <- matrix( nrow = 2, ncol = nside*nside )
		Y[1,] <- as.vector( X )
	}
	
	Files <- list.files()
	
	for( k in 1:nfile )
	{
		if( !any( Files == file[k] ) ) 
				healpix.nested.create.file( file[k], nside, nfields = 2 )
	}
	
	
	for( k in 1:nfile )
	{
	
		t <- 	.C( "astro_write_nested_healpix_patch", 
					as.character(file[k]), 
					as.double(Y[ k , ]), 
					as.integer(nside), 
					as.integer(colnum), 
					s = as.integer(0), 
					as.integer(patch),
					PACKAGE = "HealR" )
		status[k] <- t$s
	}
	

	out <- list( file = file, status = status )
	
	return( out )
}
