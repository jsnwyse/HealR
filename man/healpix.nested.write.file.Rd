\name{healpix.nested.write.file}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Write a vector to a FITS file in nested ordering.}

\description{
Write a vector to a FITS file, when the vector is already in the correct ordering. Note that this function does no re-ordering of pixels and is meant only for directly accessing the FITS files.
}
\usage{
healpix.nested.write.file( x = NULL, file = NULL, colnum = 1 )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{x}{The vector to write to file.}
  \item{file}{
The name of the file to be created.
}
\item{colnum}{The column (or field) in the data table to write to.}
}
 





\author{
Jason Wyse
}


