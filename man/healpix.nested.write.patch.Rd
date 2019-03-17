
\name{healpix.nested.write.patch}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Write an image to a Healpix nested FITS file.
}
\description{
Write an image (images) to a Healpix nested FITS file. The image can be passed as a matrix (list of matrices), and this is then re-ordered appropriately according to the Healpix nested ordering. 
}
\usage{
healpix.nested.write.patch( X = NULL, file = NULL, patch = 0, colnum = 1 )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
A matrix or list of matrices representing the image which are to be written to file. The nside parameter is determined by taking the number of rows in X (first entry of X).
}
\item{file}{A file name or vector of file names in the case of multiple images.}
  \item{patch}{
Which base patch to write to in the file.} 
\item{colnum}{
Which field to write to in the FITS file.}
}

 


\author{
Jason Wyse
}


