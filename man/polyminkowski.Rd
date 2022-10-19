\name{polyminkowski}
\alias{polyminkowski}
\title{
  Minkowski Sum of Polygon with Path
}
\description{
  Compute the Minkowski Sum of a polygon with
  a path or paths. 
}
\usage{
  polyminkowski(A, B, \dots, eps, x0, y0,
               closed=TRUE)
}
\arguments{
  \item{A}{
    Data specifying a polygon or polygons. See Details.
  }
  \item{B}{
    Data specifying a path or paths. See Details.
  }
  \item{\dots}{
    Ignored.
  }
  \item{eps}{Spatial resolution for coordinates.}
  \item{x0,y0}{Spatial origin for coordinates.}
  \item{closed}{
    Logical value indicating whether the resulting path
    should be interpreted as closed (last vertex joined to first vertex
    and interior filled in).
  }
}
\details{
  This is an interface to the function \code{MinkowskiSum} in
  Angus Johnson's \code{C++} library \pkg{Clipper}.

  It computes the Minkowski Sum of the polygon \code{A}
  (including its interior) with the path or paths \code{B} 
  (\emph{excluding} their interior). 

  The argument \code{A} should be 
  a list containing two components \code{x} and \code{y}
  giving the coordinates of the vertices of a single polygon.
  The last vertex should
  not repeat the first vertex.

  The argument \code{B} should be either
  \itemize{
    \item
    a list containing two components \code{x} and \code{y}
    giving the coordinates of the vertices of a single polygon or path.
    The last vertex should
    not repeat the first vertex.
    \item
    a \code{list} of  \code{list(x,y)} structures giving
    the coordinates of the vertices of several polygons or paths.
  }

  \bold{Calculations are performed in integer arithmetic}
  after subtracting \code{x0,y0} from the coordinates,
  dividing by \code{eps}, and rounding to the nearest integer.
  Thus, \code{eps} is the effective spatial resolution.
  The default values ensure reasonable accuracy.
}
\value{
  Data specifying polygons, in the same format as \code{A}.
}
\author{Angus Johnson.
  Ported to \R by Adrian Baddeley
  \email{Adrian.Baddeley@curtin.edu.au}.
}
\references{
  Clipper Website: \url{http://www.angusj.com}
}
\seealso{
  \code{\link{polyclip}},
  \code{\link{polyoffset}},
  \code{\link{polylineoffset}},
  \code{\link{polysimplify}}
}
\examples{
  A <- list(list(x=c(-3,3,3,-3),y=c(-3,-3,3,3)))
  B <- list(list(x=c(-1,1,1,-1),y=c(-1,-1,1,1)))
  C <- polyminkowski(A, B)
  opa <- par(mfrow=c(1,3))
  rr <- c(-4, 4)
  plot(rr, rr, type="n", axes=FALSE, xlab="", ylab="", main="A")
  polygon(A[[1]], col="blue")
  plot(rr,rr, type="n", axes=FALSE, xlab="", ylab="", main="B")
  polygon(B[[1]], border="red", lwd=4)
  plot(rr,rr, type="n", axes=FALSE, xlab="", ylab="", main="A+B")
  polygon(C[[1]], col="green")
  polygon(C[[2]], col="white")
  par(opa)
}
\keyword{spatial}
\keyword{manip}

