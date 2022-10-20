\name{polyoffset}
\alias{polyoffset}
\title{Polygon Offset}
\description{
  Given a polygonal region, compute the offset region (aka: guard region,
  buffer region, morphological dilation) formed by shifting the
  boundary outwards by a specified distance. 
}
\usage{
 polyoffset(A, delta,
         \dots,
         eps, x0, y0,
         miterlim=2, arctol=abs(delta)/100,
         jointype=c("square", "round", "miter"))
}
\arguments{
  \item{A}{Data specifying polygons. See Details.}
  \item{delta}{Distance over which the boundary should be shifted.}
  \item{\dots}{Ignored.}
  \item{eps}{Spatial resolution for coordinates.}
  \item{x0,y0}{Spatial origin for coordinates.}
  \item{miterlim,arctol}{Tolerance parameters: see Details.}
  \item{jointype}{
    Type of join operation to be performed at each vertex.
    See Details.
  }
}
\value{
  Data specifying polygons, in the same format as \code{A}.
}
\details{
  This is part of an interface to the polygon-clipping library
  \code{Clipper} written by Angus Johnson.

  Given a polygonal region \code{A},
  the function \code{polyoffset} computes the offset region
  (also known as the morphological dilation, guard region,
  buffer region, etc) obtained by shifting the boundary of \code{A}
  outward by the distance \code{delta}.

  The argument \code{A} represents a region in the 
  Euclidean plane bounded by closed polygons. The format is either
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
  
  The argument \code{jointype} determines what happens at the convex vertices
  of \code{A}. See the Examples for illustrations.
  \itemize{
    \item \code{jointype="round"}: a circular arc is generated.
    \item \code{jointype="square"}: the circular arc is
    replaced by a single straight line.
    \item \code{jointype="miter"}: the circular arc is
    omitted entirely, or replaced by a single straight line.
  }

  The arguments \code{miterlim} and \code{arctol} are tolerances.
  \itemize{
    \item if \code{jointype="round"}, then \code{arctol} is the maximum
    permissible distance between the true circular arc and its
    discretised approximation. 
    \item if \code{jointype="miter"}, then \code{miterlimit * delta}
    is the maximum permissible displacement between the original vertex and the
    corresponding offset vertex if the circular arc were to be
    omitted entirely. The default is \code{miterlimit=2}
    which is also the minimum value.
  }

  \bold{Calculations are performed in integer arithmetic}
  after subtracting \code{x0,y0} from the coordinates,
  dividing by \code{eps}, and rounding to the nearest 64-bit integer.
  Thus, \code{eps} is the effective spatial resolution.
  The default values ensure reasonable accuracy.
}
\author{Angus Johnson.
  Ported to \R by Adrian Baddeley
  \email{Adrian.Baddeley@curtin.edu.au}.
}
\seealso{
  \code{\link{polylineoffset}},
  \code{\link{polyclip}},
  \code{\link{polysimplify}},
  \code{\link{polyminkowski}}
}
\examples{
  A <- list(list(x=c(4,8,8,2,6), y=c(3,3,8,8,6)))
  plot(c(0,10),c(0,10), type="n", main="jointype=square", axes=FALSE, xlab="", ylab="")
  polygon(A[[1]], col="grey")
  C <- polyoffset(A, 1, jointype="square")
  polygon(C[[1]], lwd=3, border="blue")
  plot(c(0,10),c(0,10), type="n", main="jointype=round", axes=FALSE, xlab="", ylab="")
  polygon(A[[1]], col="grey")
  C <- polyoffset(A, 1, jointype="round")
  polygon(C[[1]], lwd=3, border="blue")
  plot(c(0,10),c(0,10), type="n", main="jointype=miter", axes=FALSE, xlab="", ylab="")
  polygon(A[[1]], col="grey")
  C <- polyoffset(A, 1, jointype="miter")
  polygon(C[[1]], lwd=3, border="blue")
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

  Chen, X. and McMains, S. (2005)
  Polygon Offsetting by Computing Winding Numbers.
  Paper no. DETC2005-85513 in \emph{Proceedings of IDETC/CIE 2005}
  (ASME 2005 International Design Engineering Technical Conferences
  and Computers and Information in Engineering Conference),
  pp. 565--575
  \url{https://mcmains.me.berkeley.edu/pubs/DAC05OffsetPolygon.pdf}
}
\keyword{spatial}
\keyword{math}
