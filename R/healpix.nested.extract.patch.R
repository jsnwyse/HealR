
healpix.nested.extract.patch <- function( file = NULL, nside = NULL, rows = NULL, cols = NULL, patch = 0, colnum=1, plot = TRUE, mask = NULL, main="", output.Matrix=FALSE )
{
  if( is.null(file) ) stop("please provide a file name with the fits extension")
  
  if( is.null(nside) ) stop("please provide the Healpix NESTED ordering resolution e.g. res 9 => nside = 512")
  
  all.rows <- 1:nside - 1
  all.cols <- 1:nside - 1
  
  #nrows <- length(rows)
  #ncols <- length(cols)
  
  x <- numeric( nside * nside )
  
  hduno <- 2
  
  nfile <- length(file)
  
  status.files <- numeric(nfile)
  
  # if subsetting, determine index sets corresponding
  #    to the specified rows and columns
  if( !is.null(rows) | !is.null(cols) )
  {
    # get rows and sort
    if( !is.null(rows) ) rows <- sort( rows ) else rows <- 1:nside
    # get cols and sort
    if( !is.null(cols) ) cols <- sort( cols ) else cols <- 1:nside
    
    nrows <- length( rows )
    ncols <- length( cols )
    
    a <- rep.int( (cols-1) * nside, times = rep( nrows, ncols ) )
    b <- rep( rows, ncols ) 
    idx <- a + b  
    
    sub.arr <- TRUE
  }else{
    nrows <- nside
    ncols <- nside
    sub.arr <- FALSE
  }
  
  if( !is.null(mask) )
  {
    status <- rep(0,2)
    w0 <- .C( "astro_extract_nested_healpix_patch", 		
              as.character(mask),
              as.integer(nside),							
              as.integer(patch),
              as.integer(nside), 		
              as.integer(all.rows),
              as.integer(nside), 							
              as.integer(all.cols),
              x = as.double(x),							
              status = as.integer(status),
              as.integer(hduno), 
              as.integer(1),
              PACKAGE="HealR" )	
    if( w0$status[1] == 104 ) stop(paste0("file '", mask ,"' either does not exist or could not be opened"))
    msk <- w0$x
  }else{
    msk <- rep( 1, nside * nside )
  }
  
  #the list to return
  out <- list( filenames = file, nside=nside, X = vector( length=nfile, mode = "list" ) )
  
  for( k in 1:nfile )
  {
    if(!is.na(file[k]))
    {
      status <- rep(0,2)
      w <- .C( "astro_extract_nested_healpix_patch", 		
               as.character(file[k]),
               as.integer(nside),							
               as.integer(patch),
               as.integer(nside), 		
               as.integer(all.rows),
               as.integer(nside), 							
               as.integer(all.cols),
               x = as.double(x),							
               status = as.integer(status),
               as.integer(hduno),
               as.integer(colnum),
               PACKAGE="HealR" )	
      status.files[k] <- w$status[1]
      if( w$status[1] == 104 ) stop(paste0("file '", file[k] ,"' either does not exist or could not be opened"))
      if( w$status[1] == 302 ) stop(paste0("file '", file[k] ,": requested column number exceeds number of fields in file"))
      
      if( sub.arr ) 
        out$X[[k]] <- w$x[idx] * msk[idx]
      else
        out$X[[k]] <- w$x * msk #vector in column major form
    }
  }
  
  if( sub.arr ) out$msk <- msk[idx] else out$msk <- msk
  
  if( output.Matrix )
  {
    for( k in 1:nfile ) out$X[[k]] <- Matrix( out$X[[k]], nrow=nrows, ncol=ncols )
    out$msk <- Matrix( out$msk, nrow=nrows, ncol=ncols )
  }
  
  if( sub.arr )
  {
    out$subim <- TRUE
    out$subim.idx <- list( rows=rows, cols=cols )
    out$subim.dim <- c( nrows, ncols )
  }else{
    out$subim <- FALSE
  }
  
  out$output.Matrix <- output.Matrix
  
  class(out) <- "HealR"
  
  return( out )
  
}
