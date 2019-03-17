
astroFITS.scale.template.source <- function( source=NULL, nu0, nu, beta.s, beta.d )
{

	if( is.null( source ) ) warning("\n no source argument given, only scaling matrix will be returned")

	log10Planck <- -33.1787441114855639057
	log10Boltzmann <- -22.859916308396282858
	
	#compute two constants
	c0 <- 10^( log10Planck + log10(nu) + 9 - log10Boltzmann - log10(2.725) )
	c1 <- 10^(  log10Planck + log10(nu) + 9 - log10Boltzmann - log10(18.1) )

	scaling <- matrix( nrow=4, ncol=length(nu) )

	scaling[1,] <- exp(c0) * ( c0/expm1(c0) )^2
	scaling[1,] <- scaling[1,]/scaling[1,1]
	scaling[2,] <- (nu/nu[1])^beta.s
	scaling[3,] <- ( expm1(c1[5])/expm1(c1) ) * ( nu/nu[5] )^(1+beta.d)
	scaling[4,] <- (nu/nu[1])^(-2.14)
	
	out <- list()
	
	out$X <- list()
	
	if( !is.null(source) )
	{
	
		nside <- nrow(source$X[[1]])
	
		for( k in 1:length(nu) )
		{
			X <- Matrix( 0, nrow = nside, ncol=nside )  
			for( l in 1:4 )
			{
				X <- X + scaling[l,k] * source$X[[l]]
			}
			out$X[[k]] <- X
		}
		
		out$scalemat <- scaling
	
	}else{
	
		out$scalemat <- scaling
	
	}
	
	return(out)
	
}


HealR.scale.template.source <- function( source=NULL, nu0, nu, beta.s, beta.d )
{

	if( is.null( source ) ) warning("\n no source argument given, only scaling matrix will be returned")
	
	#compute two constants

	scaling <- matrix( nrow=4, ncol=length(nu) )
	
	scaling[1,] <- c( 1, 1, 1, 1, 1 )
	scaling[2,] <- c( 2, 1, 3, 4, 5 )
	scaling[3,] <- c( 5, 4, 3, 2, 1 )
	scaling[4,] <- c( 2, 2, 1, 2, 2 )

	scaling[2,] <- (scaling[2,])^beta.s
	scaling[3,] <- ( scaling[3,])^(1+beta.d)
	scaling[4,] <- (scaling[4,])^(-2)
	
	out <- list()
	
	out$X <- list()
	
	if( !is.null(source) )
	{
	
		nside <- nrow(source$X[[1]])
	
		for( k in 1:length(nu) )
		{
			X <- Matrix( 0, nrow = nside, ncol=nside )  
			for( l in 1:4 )
			{
				X <- X + scaling[l,k] * source$X[[l]]
			}
			out$X[[k]] <- X
		}
		
		out$scalemat <- scaling
	
	}else{
	
		out$scalemat <- scaling
	
	}
	
	return(out)
	
}

