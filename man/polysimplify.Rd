\name{polysimplify}
\alias{polysimplify}
\title{
  Remove Self-Intersections from a Polygon
}
\description{
  This function attempts to remove self-intersections and duplicated
  vertices from the given polygon. 
}
\usage{
  polysimplify(A, \dots, eps, x0, y0,
               filltype = c("evenodd", "nonzero", "positive", "negative"))
}
\arguments{
  \item{A}{
    Data specifying a polygon or polygons. See Details.
  }
  \item{\dots}{
    Ignored.
  }
  \item{eps}{Spatial resolution for coordinates.}
  \item{x0,y0}{Spatial origin for coordinates.}
  \item{filltype}{Polygon-filling rule. See Details.}
}
\details{
  This is an interface to the function \code{SimplifyPolygons} in
  Angus Johnson's \code{C++} library \pkg{Clipper}.
  It tries to remove self-intersections from the supplied polygon,
  by performing a boolean union operation using the nominated
  \code{filltype}. The result may be one or several polygons.

  The argument \code{A} should be either
  \itemize{
    \item
    a list containing two components \code{x} and \code{y}
    giving the coordinates of the vertices of a single polygon.
    The last vertex should
    not repeat the first vertex.
    \item
    a \code{list} of  \code{list(x,y)} structures giving
    the coordinates of the vertices of several polygons.
  }

  The argument \code{filltype} determines the polygon fill type.
  \describe{
    \item{Even-Odd:}{
      The default rule is \emph{even-odd} filling,
      in which every polygon edge demarcates a boundary between
      the inside and outside of the region. 
      It does not matter whether a polygon is traversed in
      clockwise or anticlockwise order. Holes are determined
      simply by their locations relative to other polygons such that
      outers contain holes and holes contain outers. 
    }
    \item{Non-Zero:}{
      Under the \emph{nonzero} filling rule, an outer boundary must be
      traversed in clockwise order, while a hole must be traversed
      in anticlockwise order.
    }
    \item{Positive:}{
      Under the \code{positive} filling rule, the filled region
      consists of all points with positive winding number.
    }
    \item{Negative:}{
      Under the \code{negative} filling rule, the filled region
      consists of all points with negative winding number.
    }
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
  \code{\link{polyminkowski}}
}
\examples{
  theta <- 2 * pi * (0:5) * 2/5

  A <- list(list(x=sin(theta), y=cos(theta)))
  B <- polysimplify(A, filltype="nonzero")

  opa <- par(mfrow=c(1,2))
  plot(c(-1,1),c(-1,1), type="n", axes=FALSE, xlab="", ylab="")
  with(A[[1]], segments(x[-6], y[-6], x[-1], y[-1], col="red"))
  points(A[[1]], col="red")

  with(A[[1]], text(x[1:5], y[1:5], labels=1:5, cex=2))
  plot(c(-1,1),c(-1,1), type="n", axes=FALSE, xlab="", ylab="")
  polygon(B[[1]], lwd=3, col="green")
  with(B[[1]], text(x,y,labels=seq_along(x), cex=2))
  par(opa)
}
\keyword{spatial}
\keyword{manip}

