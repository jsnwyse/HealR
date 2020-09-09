plot.HealR <- function( x, show=1, main="", colorgradient=c("darkblue","white","darkred"), mask.color = "grey25", ... )
{
 	# check if cast as type 'Matrix' 
  if( x$output.Matrix )
  {
    stop("cant plot this a.t.m.")
  }
  
  if( x$subim )
  {
    nrows <- x$subim.dim[1] 
    ncols <- x$subim.dim[2] 
  }else{
    nrows <- ncols <- x$nside
  }
  
  u <- x$X[[show]] #rescale( x$X[[show]], to=c(0,1) )
  z <- u
  z[ u <= 0 ] <- rescale( u[ u <= 0 ], to=c(0,.5) )
  z[ u >= 0 ] <- rescale( u[ u >= 0 ], to=c(.5,1) )
  # set mased pixels to NA to utilise col.na
  z[ x$msk == 0 ] <- NA
  im.mat <- matrix( z, nrow=nrows )
  im <- as.cimg( im.mat )
  
  mai.orig <- par()$mai
  par( mai=rep(0,4) )
  cscale <- scales::gradient_n_pal(colorgradient,c(0,.5,1))
  plot( im, colorscale = cscale, rescale=FALSE, col.na=mask.color, axes=FALSE )
  par( mai=mai.orig )
}
