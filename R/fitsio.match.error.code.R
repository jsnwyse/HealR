fitsio.match.error.code <- function( errcode, file )
{
  err <- c( 104, 106, 301, 302 )
  msg <- c( paste0( "file '", file ,"' either does not exist\n or could not be opened." ),
            paste0( "file '", file, "' was opened sucessfully\n but error writing to column."),
            paste0( "file '", file, "' was opened sucessfully\n but unable to access the HDU." ),
            paste0( "file '", file, "' did not contain enough\n fields to access the colnum requested."))
  k <- match( errcode, err )
  return( msg[k] )
}