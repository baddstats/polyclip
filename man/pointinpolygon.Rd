\name{pointinpolygon}
\alias{pointinpolygon}
\title{Test Whether a Point Lies Inside a Polygon}
\description{
  Test whether each point lies inside a specified polygon.
}
\usage{
 pointinpolygon(P, A, eps, x0, y0)
}
\arguments{
  \item{P}{
    Spatial coordinates of the points to be tested.
    A list of two vectors named \code{x} and \code{y}.
  }
  \item{A}{
    A single polygon, specified as a list of two vectors
    named \code{x} and \code{y}.
  }
  \item{eps}{Spatial resolution for coordinates.}
  \item{x0,y0}{Spatial origin for coordinates.}
}
\value{
  An integer vector with one entry for each point in \code{P}.
  The result is 0 if the point lies outside \code{A},
  1 if the point lies inside \code{A}, and -1 if it lies on the
  boundary.
}
\details{
  This is part of an interface to the polygon-clipping library
  \code{Clipper} written by Angus Johnson.

  The argument \code{A} represents a closed polygon.
  It should be 
  a list containing two components \code{x} and \code{y}
  giving the coordinates of the vertices.
  The last vertex should
  not repeat the first vertex.

  \bold{Calculations are performed in integer arithmetic}
  after subtracting \code{x0,y0} from the coordinates,
  dividing by \code{eps}, and rounding to the nearest integer.
  Thus, \code{eps} is the effective spatial resolution.
  The default values ensure reasonable accuracy.
}
\seealso{
  \code{\link{polyclip}}.
}
\author{Angus Johnson.
  Ported to \R by Adrian Baddeley
  \email{Adrian.Baddeley@curtin.edu.au}.
}
\examples{
  A <- list(x=1:10, y=c(1:5,5:1))
  P <- list(x=4, y=2)
  pointinpolygon(P, A)
}
\references{
  Clipper Website: \url{http://www.angusj.com}

  Vatti, B. (1992) A generic solution to polygon clipping.
  \emph{Communications of the ACM} \bold{35} (7) 56--63.
  \url{https://dl.acm.org/doi/10.1145/129902.129906}

  Agoston, M.K. (2005)
  \emph{Computer graphics and geometric modeling:
    implementation and algorithms.}
  Springer-Verlag.
  \url{http://books.google.com/books?q=vatti+clipping+agoston}
}
\keyword{spatial}
\keyword{math}
