\name{healpix.nested.extract.patch}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Extract a base patch from a Healpix nested FITS file}

\description{
Extract a base patch from a Healpix FITS file in the nested ordering.
}
\usage{
healpix.nested.extract.patch( file = NULL, nside = NULL, rows = NULL,
 cols = NULL, patch = 0, plot = TRUE, mask = NULL, main="" )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{
The name of the file(s) from which the patch(es) is (are) to be extracted.
}
\item{nside}{The nside parameter giving the resolution of the file.}
  \item{rows}{
Extract only the rows that are specified (useful for sub-images of the base patch} 
\item{cols}{
Extract only the columns that are specified (useful for sub-images of the base patch}
\item{patch}{The base patch to be extracted. Detault 0. Patches are indexed 0-11.}
\item{plot}{Whether the (first) base patch should be plotted or not.}
\item{mask}{The name of the mask file if some of the pixels are to be masked out. This can improve visualisation if there is a mask routinely applied.}
\item{main}{The title to be given to the plot if \code{plot} is \code{TRUE}.}
}
 


\value{
This function returns a list whose number of entries is the number of files read. The list is of class \code{HealR}, which has an S3 method for plotting.  
}



\author{
Jason Wyse
}


