\name{polyclip}
\alias{polyclip}
\title{Polygon Clipping}
\description{
  Find intersection, union or set difference of two polygonal regions.
}
\usage{
 polyclip(A, B, op=c("intersection", "union", "minus", "xor"),
         \dots,
         eps, x0, y0,
         fillA=c("evenodd", "nonzero", "positive", "negative"),
         fillB=c("evenodd", "nonzero", "positive", "negative"),
         closed = TRUE)
}
\arguments{
  \item{A,B}{
    Data specifying polygons. See Details.
  }
  \item{op}{Set operation to be performed to combine \code{A} and \code{B}.}
  \item{\dots}{Ignored.}
  \item{eps}{Spatial resolution for coordinates.}
  \item{x0,y0}{Spatial origin for coordinates.}
  \item{fillA,fillB}{Polygon-filling rule for \code{A} and \code{B}.}
  \item{closed}{Whether \code{A} is a closed shape or an open shape (a line).}
}
\value{
  Data specifying polygons, in the same format as \code{A} and \code{B}.
}
\details{
  This is an interface to the polygon-clipping library
  \code{Clipper} written by Angus Johnson.

  Given two polygonal regions \code{A} and \code{B}
  the function \code{polyclip} performs one of the following
  geometrical operations:
  \itemize{
    \item \code{op="intersection"}: set intersection of \code{A} and \code{B}.
    \item \code{op="union"}: set union of \code{A} and \code{B}.
    \item \code{op="minus"}: set subtraction (sometimes called set difference):
    the region covered by \code{A} that is not covered by \code{B}.
    \item \code{op="xor"}: exclusive set difference (sometimes called
    exclusive-or): the region covered by exactly one of the sets
    \code{A} and \code{B}.
  }

  Each of the arguments \code{A} and \code{B} represents a region in the 
  Euclidean plane bounded by closed polygons. The format of these
  arguments is either
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
  Note that calculations are performed in integer arithmetic: see below.
  
  The interpretation of the polygons
  depends on the \emph{polygon-filling rule} for \code{A} and \code{B}
  that is specified by the arguments \code{fillA} and \code{fillB}
  respectively.
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
\seealso{
  \code{\link{polysimplify}},
  \code{\link{polyoffset}},
  \code{\link{polylineoffset}},
  \code{\link{polyminkowski}}
}
\author{Angus Johnson.
  Ported to \R by Adrian Baddeley
  \email{Adrian.Baddeley@curtin.edu.au}.
}
\examples{
  A <- list(list(x=1:10, y=c(1:5,5:1)))
  B <- list(list(x=c(2,8,8,2),y=c(0,0,10,10)))
  plot(c(0,10),c(0,10), type="n", axes=FALSE, xlab="", ylab="")
  polygon(A[[1]])
  polygon(B[[1]])
  C <- polyclip(A, B)
  polygon(C[[1]], lwd=3, col=3)
}
\references{
  Clipper Website: \url{http://www.angusj.com}

  Vatti, B. (1992) A generic solution to polygon clipping.
  \emph{Communications of the ACM} \bold{35} (7) 56--63.
  \url{http://portal.acm.org/citation.cfm?id=129906}

  Agoston, M.K. (2005)
  \emph{Computer graphics and geometric modeling:
    implementation and algorithms.}
  Springer-Verlag.
  \url{http://books.google.com/books?q=vatti+clipping+agoston}

  Chen, X. and McMains, S. (2005)
  Polygon Offsetting by Computing Winding Numbers.
  Paper no. DETC2005-85513 in \emph{Proceedings of IDETC/CIE 2005}
  (ASME 2005 International Design Engineering Technical Conferences
  and Computers and Information in Engineering Conference),
  pp. 565--575
  \url{http://www.me.berkeley.edu/~mcmains/pubs/DAC05OffsetPolygon.pdf}
}
\keyword{spatial}
\keyword{math}
