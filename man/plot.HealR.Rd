\name{plot.HealR}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot a \code{HealR} object.
}
\description{
Plot a \code{HealR} object.
}
\usage{
plot.HealR( x, show=1, main="", colorkey=TRUE, ... )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A \code{HealR} object.
}
\item{show}{Which entry of the \code{HealR} object to visualise.}
\item{main}{The title to be given to the plot if \code{plot} is \code{TRUE}.}
\item{colorkey}{If \code{TRUE} then a colorkey is given at the side of the plot.}
}
 


\value{
This function uses a direct call to \code{image.Matrix} to create a visualisation of the Healpix file.
}



\author{
Jason Wyse
}


