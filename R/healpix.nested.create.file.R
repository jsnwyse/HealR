healpix.nested.create.file <- function( file = NULL, nside = NULL, nfields = 2 )
{
	nfile <- length( file )
	status <- numeric( nfile )
	
	for( k in 1:nfile )
	{
		.C( 	"astro_create_nested_healpix_file", 
				as.character(file[k]),
				as.integer(nside),
				as.integer(nfields),
				as.integer(status[k]),
				PACKAGE = "HealR" )	
	}
	
	out <- list( file = file, status = status )
	
	return( out )
}
