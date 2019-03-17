\name{healpix.nested.create.file}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create a FITS file in Healpix nested ordering}

\description{
Creates a FITS file in Healpix nested ordering with specified nside and number of fields.
}
\usage{
healpix.nested.create.file( file = NULL, nside = NULL, nfields = 2 )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{
The name of the file to be created.
}
\item{nside}{The nside parameter giving the resolution of the file.}
  \item{nfields}{
The number of fields to create in the binary table. 
}
 
}

\value{
This function creates the file in the current working directory or in the path specified in the \code{file} argument.
}



\author{
Jason Wyse
}


