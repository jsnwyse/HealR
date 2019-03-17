healpix.nested.write.file <- function( x = NULL, file = NULL, colnum = 1 )
{
	
	status <- 0
	
	nside <- sqrt( length(x) / 12 )
	
	.C( 	"astro_write_vector_to_healpix_file_0", 
			as.character(file), 
			as.double(x), 
			as.integer(nside), 
			as.integer(colnum), 
			as.integer(status),
			PACKAGE = "HealR" )

	out <- list( file = file, status = status )
	
	return( out )
}
