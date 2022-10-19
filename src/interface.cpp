#include "clipper.h"
#include <R.h>
#include <Rdefines.h>

using namespace std;
using namespace ClipperLib;

void CopyToPath(int *x, int *y, int n, ClipperLib::Path &p)
{
  p.clear();
  p.reserve(n);
  for (int i = 0; i < n; i++)
    p.push_back(IntPoint(x[i], y[i]));
}

void CopyFromPath(ClipperLib::Path &p, int *x, int *y, int nmax, int *n)
{
  int N;
  *n = N = p.size();
  if(N <= nmax) {
    for (int i = 0; i < N; i++) 
      {
	x[i] = p[i].X;
	y[i] = p[i].Y;
      }
  }
}

void ScaleToPath(double *x, double *y, int n, ClipperLib::Path &p,
                 double x0, double y0, double eps)
{
  int i;
  cInt cxi, cyi;
  p.clear();
  p.reserve(n);
  for (i = 0; i < n; i++) {
    cxi = (cInt) ((x[i] - x0)/eps);
    cyi = (cInt) ((y[i] - y0)/eps);
    p.push_back(IntPoint(cxi, cyi));
  }
}

void ScaleFromPath(ClipperLib::Path &p, double *x, double *y, int nmax, int *n,
		  double x0, double y0, double eps)
{
  int N;
  *n = N = p.size();
  if(N <= nmax) {
    for (int i = 0; i < N; i++) 
      {
	x[i] = x0 + eps * ((double) p[i].X);
	y[i] = y0 + eps * ((double) p[i].Y);
      }
  }
}

void ScaleToPoint(double x, double y, ClipperLib::IntPoint &p,
		  double x0, double y0, double eps)
{
  p.X = (cInt) ((x - x0)/eps);
  p.Y = (cInt) ((y - y0)/eps);
}

extern "C" {
  SEXP Csimplify(SEXP A,
		 SEXP pft,
		 SEXP X0,
		 SEXP Y0,
		 SEXP Eps) {
    int nA, i, nAi, m, mi, mitrue;
    double *x, *y, *xx, *yy;
    SEXP Ai = R_NilValue;
    SEXP out, outi, xouti, youti;
    int pftcode;
    PolyFillType filltype;
    double x0, y0, eps;

    // protect arguments from garbage collector    
    PROTECT(A   = AS_LIST(A));
    PROTECT(pft  = AS_INTEGER(pft));
    PROTECT(X0  = AS_NUMERIC(X0));
    PROTECT(Y0  = AS_NUMERIC(Y0));
    PROTECT(Eps = AS_NUMERIC(Eps));
    // that's 5 arguments

    // number of polygons
    nA = LENGTH(A);

    // Initialise object containing n polygons
    Paths polyA(nA);

    // Get scale parameters
    x0 = *(NUMERIC_POINTER(X0));
    y0 = *(NUMERIC_POINTER(Y0));
    eps = *(NUMERIC_POINTER(Eps));

    // copy data
    for(i = 0; i < nA; i++) {
      Ai = VECTOR_ELT(A, i);
      nAi = LENGTH(VECTOR_ELT(Ai, 0));
      x = NUMERIC_POINTER(VECTOR_ELT(Ai, 0));
      y = NUMERIC_POINTER(VECTOR_ELT(Ai, 1));
      ScaleToPath(x, y, nAi, polyA[i], x0, y0, eps);
    }

    // interpret clipping parameters
    pftcode = *(INTEGER_POINTER(pft));
    switch(pftcode) {
    case 1: 
      filltype = pftEvenOdd; 
      break;
    case 2:
      filltype = pftNonZero;
      break;
    case 3:
      filltype = pftPositive;
      break;
    case 4:
      filltype = pftNegative;
      break;
    default: 
      error("polyclip: unrecognised code for fill type A");
    }

    // simplify polygon;
    Paths result;
    SimplifyPolygons(polyA, result, filltype);

    // number of polygons
    m = result.size();
    
    // initialise output list
    PROTECT(out  = NEW_LIST(m));
    
    // copy data
    if(m > 0) {
      for(i = 0; i < m; i++) {
	mi = result[i].size();
	// Allocate space for output
	PROTECT(outi = NEW_LIST(2));
	PROTECT(xouti = NEW_NUMERIC(mi));
	PROTECT(youti = NEW_NUMERIC(mi));
	xx = NUMERIC_POINTER(xouti);
	yy = NUMERIC_POINTER(youti);
	// copy to output space
	ScaleFromPath(result[i], xx, yy, mi, &mitrue, x0, y0, eps);
	// Put vectors into list
	SET_VECTOR_ELT(outi, 0, xouti);
	SET_VECTOR_ELT(outi, 1, youti);
	SET_VECTOR_ELT(out, i, outi);
      }
    }

    UNPROTECT(6 + 3*m); // 5 arguments + out + m * (outi, xouti, youti)
    return(out);
  }
}

// -----------------------------------------------------------------

extern "C" {
  SEXP Cclipbool(SEXP A,
		 SEXP B,
		 SEXP pftA,
		 SEXP pftB,
		 SEXP ct,
		 SEXP X0,
                 SEXP Y0,
                 SEXP Eps,
                 SEXP clo // whether paths in A are closed
		 ){ 
    int nA, nB, i, n, m, mi, mitrue;
    bool closed;
    double *x, *y, *xx, *yy;
    SEXP Ai = R_NilValue, Bi = R_NilValue;
    SEXP out, outi, xouti, youti;
    ClipType cliptype;
    PolyFillType filltypeA, filltypeB;
    int ctcode, pftAcode, pftBcode;
    double x0, y0, eps;
    
    // protect arguments from garbage collector    
    PROTECT(A   = AS_LIST(A));
    PROTECT(B   = AS_LIST(B));
    PROTECT(clo = AS_LOGICAL(clo));
    PROTECT(ct  = AS_INTEGER(ct));
    PROTECT(pftA  = AS_INTEGER(pftA));
    PROTECT(pftB  = AS_INTEGER(pftB));
    PROTECT(X0  = AS_NUMERIC(X0));
    PROTECT(Y0  = AS_NUMERIC(Y0));
    PROTECT(Eps = AS_NUMERIC(Eps));

    // lengths of lists
    nA = LENGTH(A);
    nB = LENGTH(B);

    // Initialise object containing n polygons
    Paths polyA(nA), polyB(nB);
    closed = *(LOGICAL_POINTER(clo));

    // Get scale parameters
    x0 = *(NUMERIC_POINTER(X0));
    y0 = *(NUMERIC_POINTER(Y0));
    eps = *(NUMERIC_POINTER(Eps));

    // copy data
    for(i = 0; i < nA; i++) {
      Ai = VECTOR_ELT(A, i);
      n = LENGTH(VECTOR_ELT(Ai, 0));
      x = NUMERIC_POINTER(VECTOR_ELT(Ai, 0));
      y = NUMERIC_POINTER(VECTOR_ELT(Ai, 1));
      ScaleToPath(x, y, n, polyA[i], x0, y0, eps);
    }
    for(i = 0; i < nB; i++) {
      Bi = VECTOR_ELT(B, i);
      n = LENGTH(VECTOR_ELT(Bi, 0));
      x = NUMERIC_POINTER(VECTOR_ELT(Bi, 0));
      y = NUMERIC_POINTER(VECTOR_ELT(Bi, 1));
      ScaleToPath(x, y, n, polyB[i], x0, y0, eps);
    }

    // interpret clipping parameters
    ctcode = *(INTEGER_POINTER(ct));
    pftAcode = *(INTEGER_POINTER(pftA));
    pftBcode = *(INTEGER_POINTER(pftB));
    switch(ctcode) {
    case 1: 
      cliptype = ctIntersection; 
      break;
    case 2:
      cliptype = ctUnion;
      break;
    case 3:
      cliptype = ctDifference;
      break;
    case 4:
      cliptype = ctXor;
      break;
    default: 
      error("polyclip: unrecognised code for cliptype");
    }
    switch(pftAcode) {
    case 1: 
      filltypeA = pftEvenOdd; 
      break;
    case 2:
      filltypeA = pftNonZero;
      break;
    case 3:
      filltypeA = pftPositive;
      break;
    case 4:
      filltypeA = pftNegative;
      break;
    default: 
      error("polyclip: unrecognised code for fill type A");
    }
    switch(pftBcode) {
    case 1: 
      filltypeB = pftEvenOdd; 
      break;
    case 2:
      filltypeB = pftNonZero;
      break;
    case 3:
      filltypeB = pftPositive;
      break;
    case 4:
      filltypeB = pftNegative;
      break;
    default: 
      error("polyclip: unrecognised code for fill type B");
    }

    // perform clipping operation
    Clipper c;
    Paths result;
    c.AddPaths(polyA, ptSubject, closed);
    c.AddPaths(polyB, ptClip, true);

    if (closed) {
        c.Execute(cliptype, result, filltypeA, filltypeB);
    } else {
        PolyTree polyTreeResult;
        c.Execute(cliptype, polyTreeResult, filltypeA, filltypeB);
        OpenPathsFromPolyTree(polyTreeResult, result);
    }

    // number of polygons
    m = result.size();
    
    // initialise output list
    PROTECT(out  = NEW_LIST(m));
    
    // copy data
    if(m > 0) {
      for(i = 0; i < m; i++) {
	mi = result[i].size();
	// Allocate space for output
	PROTECT(outi = NEW_LIST(2));
	PROTECT(xouti = NEW_NUMERIC(mi));
	PROTECT(youti = NEW_NUMERIC(mi));
	xx = NUMERIC_POINTER(xouti);
	yy = NUMERIC_POINTER(youti);
	// copy to output space
	ScaleFromPath(result[i], xx, yy, mi, &mitrue, x0, y0, eps);
	// Put vectors into list
	SET_VECTOR_ELT(outi, 0, xouti);
	SET_VECTOR_ELT(outi, 1, youti);
	SET_VECTOR_ELT(out, i, outi);
      }
    }

    UNPROTECT(10 + 3*m); // 9 arguments + out + m * (outi, xouti, youti)
    return(out);
  }
}

// offset (dilation) operation for closed polygons

extern "C" {
  SEXP Cpolyoffset(SEXP A,
		   SEXP del,
		   SEXP jt,
		   SEXP mlim,
		   SEXP atol,
		   SEXP X0,
		   SEXP Y0,
		   SEXP Eps
		 ){ 
    int nA, i, n, m, mi, mitrue;
    double *x, *y, *xx, *yy;
    SEXP Ai = R_NilValue;
    SEXP out, outi, xouti, youti;
    JoinType jointype;
    int jtcode;
    double delta, miterlimit, arctolerance;
    double x0, y0, eps;
    
    // protect arguments from garbage collector    
    PROTECT(A   = AS_LIST(A));
    PROTECT(del = AS_NUMERIC(del));
    PROTECT(jt  = AS_INTEGER(jt));
    PROTECT(mlim = AS_NUMERIC(mlim));
    PROTECT(atol = AS_NUMERIC(atol));
    PROTECT(X0  = AS_NUMERIC(X0));
    PROTECT(Y0  = AS_NUMERIC(Y0));
    PROTECT(Eps = AS_NUMERIC(Eps));

    // length of list
    nA = LENGTH(A);

    // Initialise object containing nA polygons
    Paths polyA(nA);

    // Get scale parameters
    x0 = *(NUMERIC_POINTER(X0));
    y0 = *(NUMERIC_POINTER(Y0));
    eps = *(NUMERIC_POINTER(Eps));

    // copy data
    for(i = 0; i < nA; i++) {
      Ai = VECTOR_ELT(A, i);
      n = LENGTH(VECTOR_ELT(Ai, 0));
      x = NUMERIC_POINTER(VECTOR_ELT(Ai, 0));
      y = NUMERIC_POINTER(VECTOR_ELT(Ai, 1));
      ScaleToPath(x, y, n, polyA[i], x0, y0, eps);
    }

    // interpret offset parameters
    jtcode = *(INTEGER_POINTER(jt));
    switch(jtcode) {
    case 1: 
      jointype = jtSquare; 
      break;
    case 2:
      jointype = jtRound;
      break;
    case 3:
      jointype = jtMiter;
      break;
    default: 
      error("polyclip: unrecognised code for jointype");
    }

    // get parameters
    delta = *(NUMERIC_POINTER(del));   // absolute distance
    miterlimit = *(NUMERIC_POINTER(mlim));   // multiple of 'delta'
    arctolerance = *(NUMERIC_POINTER(atol));   // absolute distance
    // rescale
    delta = delta/eps;
    arctolerance = arctolerance/eps;

    // perform offset operation
    ClipperOffset co;
    Paths result;
    co.AddPaths(polyA, jointype, etClosedPolygon);
    co.MiterLimit = miterlimit;
    co.ArcTolerance = arctolerance;
    co.Execute(result, delta);

    // number of polygons
    m = result.size();
    
    // initialise output list
    PROTECT(out  = NEW_LIST(m));
    
    // copy data
    if(m > 0) {
      for(i = 0; i < m; i++) {
	mi = result[i].size();
	// Allocate space for output
	PROTECT(outi = NEW_LIST(2));
	PROTECT(xouti = NEW_NUMERIC(mi));
	PROTECT(youti = NEW_NUMERIC(mi));
	xx = NUMERIC_POINTER(xouti);
	yy = NUMERIC_POINTER(youti);
	// copy to output space
	ScaleFromPath(result[i], xx, yy, mi, &mitrue, x0, y0, eps);
	// Put vectors into list
	SET_VECTOR_ELT(outi, 0, xouti);
	SET_VECTOR_ELT(outi, 1, youti);
	SET_VECTOR_ELT(out, i, outi);
      }
    }

    UNPROTECT(9 + 3*m); // 8 arguments + out + m * (outi, xouti, youti)
    return(out);
  }
}


// offset (dilation) operation for polygonal lines

extern "C" {
  SEXP Clineoffset(SEXP A,
		   SEXP del,
		   SEXP jt,
		   SEXP et,
		   SEXP mlim,
		   SEXP atol,
		   SEXP X0,
		   SEXP Y0,
		   SEXP Eps
	 ){ 
    int nA, i, n, m, mi, mitrue;
    double *x, *y, *xx, *yy;
    SEXP Ai = R_NilValue;
    SEXP out, outi, xouti, youti;
    JoinType jointype;
    EndType endtype;
    int jtcode, etcode;
    double delta, miterlimit, arctolerance;
    double x0, y0, eps;
    
    // protect arguments from garbage collector    
    PROTECT(A   = AS_LIST(A));
    PROTECT(del = AS_NUMERIC(del));
    PROTECT(jt  = AS_INTEGER(jt));
    PROTECT(et  = AS_INTEGER(et));
    PROTECT(mlim = AS_NUMERIC(mlim));
    PROTECT(atol = AS_NUMERIC(atol));
    PROTECT(X0  = AS_NUMERIC(X0));
    PROTECT(Y0  = AS_NUMERIC(Y0));
    PROTECT(Eps = AS_NUMERIC(Eps));

    // length of list
    nA = LENGTH(A);

    // Initialise object containing nA polygonal lines
    Paths polyA(nA);

    // Get scale parameters
    x0 = *(NUMERIC_POINTER(X0));
    y0 = *(NUMERIC_POINTER(Y0));
    eps = *(NUMERIC_POINTER(Eps));

    // copy data
    for(i = 0; i < nA; i++) {
      Ai = VECTOR_ELT(A, i);
      n = LENGTH(VECTOR_ELT(Ai, 0));
      x = NUMERIC_POINTER(VECTOR_ELT(Ai, 0));
      y = NUMERIC_POINTER(VECTOR_ELT(Ai, 1));
      ScaleToPath(x, y, n, polyA[i], x0, y0, eps);
    }

    // interpret offset parameters
    jtcode = *(INTEGER_POINTER(jt));
    switch(jtcode) {
    case 1: 
      jointype = jtSquare; 
      break;
    case 2:
      jointype = jtRound;
      break;
    case 3:
      jointype = jtMiter;
      break;
    default: 
      error("polyclip: unrecognised code for jointype");
    }
    etcode = *(INTEGER_POINTER(et));
    switch(etcode) {
    case 1: 
      endtype = etClosedPolygon; 
      break;
    case 2:
      endtype = etClosedLine;
      break;
    case 3:
      endtype = etOpenButt;
      break;
    case 4:
      endtype = etOpenSquare;
      break;
    case 5:
      endtype = etOpenRound;
      break;
    default: 
      error("polyclip: unrecognised code for endtype");
    }

    // get parameters
    delta = *(NUMERIC_POINTER(del));   // absolute distance
    miterlimit = *(NUMERIC_POINTER(mlim));   // multiple of 'delta'
    arctolerance = *(NUMERIC_POINTER(atol));   // absolute distance
    // rescale
    delta = delta/eps;
    arctolerance = arctolerance/eps;

    // perform offset operation
    ClipperOffset co;
    Paths result;
    co.AddPaths(polyA, jointype, endtype);
    co.MiterLimit = miterlimit;
    co.ArcTolerance = arctolerance;
    co.Execute(result, delta);

    // number of polygons
    m = result.size();
    
    // initialise output list
    PROTECT(out  = NEW_LIST(m));
    
    // copy data
    if(m > 0) {
      for(i = 0; i < m; i++) {
	mi = result[i].size();
	// Allocate space for output
	PROTECT(outi = NEW_LIST(2));
	PROTECT(xouti = NEW_NUMERIC(mi));
	PROTECT(youti = NEW_NUMERIC(mi));
	xx = NUMERIC_POINTER(xouti);
	yy = NUMERIC_POINTER(youti);
	// copy to output space
	ScaleFromPath(result[i], xx, yy, mi, &mitrue, x0, y0, eps);
	// Put vectors into list
	SET_VECTOR_ELT(outi, 0, xouti);
	SET_VECTOR_ELT(outi, 1, youti);
	SET_VECTOR_ELT(out, i, outi);
      }
    }

    UNPROTECT(10 + 3*m); // 9 arguments + out + m * (outi, xouti, youti)
    return(out);
  }
}

// Minkowski sum of polygon with **path(s)** 

extern "C" {
  SEXP Cminksum(SEXP A,            // list(list(x,y)) : polygon
		SEXP B,            // list(list(x,y), list(x,y), ....)
		SEXP clo,          // whether paths in B are closed
		SEXP X0,
		SEXP Y0,
		SEXP Eps) {
    int nB, i, nBi, nA0, m, mi, mitrue;
    double *x, *y, *xx, *yy;
    SEXP A0 = R_NilValue;
    SEXP Bi = R_NilValue;
    SEXP out, outi, xouti, youti;
    bool closed;
    double x0, y0, eps;
    Path pathA;

    // protect arguments from garbage collector    
    PROTECT(A   = AS_LIST(A));
    PROTECT(B   = AS_LIST(B));
    PROTECT(clo = AS_LOGICAL(clo));
    PROTECT(X0  = AS_NUMERIC(X0));
    PROTECT(Y0  = AS_NUMERIC(Y0));
    PROTECT(Eps = AS_NUMERIC(Eps));
    // that's 6 arguments

    // Get scale parameters
    x0 = *(NUMERIC_POINTER(X0));
    y0 = *(NUMERIC_POINTER(Y0));
    eps = *(NUMERIC_POINTER(Eps));

    // logical value specifying whether paths in B should be closed
    closed = *(LOGICAL_POINTER(clo));

    // copy data from A
    A0 = VECTOR_ELT(A, 0);
    nA0 = LENGTH(VECTOR_ELT(A0, 0));
    x = NUMERIC_POINTER(VECTOR_ELT(A0, 0));
    y = NUMERIC_POINTER(VECTOR_ELT(A0, 1));
    ScaleToPath(x, y, nA0, pathA, x0, y0, eps);

    // number of polygons in B
    nB = LENGTH(B);
    // Initialise object representing nB polygons
    Paths pathsB(nB);

    // copy data from B
    for(i = 0; i < nB; i++) {
      Bi = VECTOR_ELT(B, i);
      nBi = LENGTH(VECTOR_ELT(Bi, 0));
      x = NUMERIC_POINTER(VECTOR_ELT(Bi, 0));
      y = NUMERIC_POINTER(VECTOR_ELT(Bi, 1));
      ScaleToPath(x, y, nBi, pathsB[i], x0, y0, eps);
    }

    // hit it
    Paths result;
    MinkowskiSum(pathA, pathsB, result, closed);

    // number of polygons
    m = result.size();
    
    // initialise output list
    PROTECT(out  = NEW_LIST(m));

    // adjust origin: (x0,y0) were subtracted from both A and B
    x0 = 2.0 * x0;
    y0 = 2.0 * y0;

    // copy data
    if(m > 0) {
      for(i = 0; i < m; i++) {
	mi = result[i].size();
	// Allocate space for output
	PROTECT(outi = NEW_LIST(2));
	PROTECT(xouti = NEW_NUMERIC(mi));
	PROTECT(youti = NEW_NUMERIC(mi));
	xx = NUMERIC_POINTER(xouti);
	yy = NUMERIC_POINTER(youti);
	// copy to output space
	ScaleFromPath(result[i], xx, yy, mi, &mitrue, x0, y0, eps);
	// Put vectors into list
	SET_VECTOR_ELT(outi, 0, xouti);
	SET_VECTOR_ELT(outi, 1, youti);
	SET_VECTOR_ELT(out, i, outi);
      }
    }

    UNPROTECT(7 + 3*m); // 6 arguments + out + m * (outi, xouti, youti)
    return(out);
  }
}

// point in polygon test ///////////////////

extern "C" {
  SEXP Cpiptest(SEXP P,  // test points, list(x, y)
		SEXP A,  // single polygon, list(x,y) 
		SEXP X0,
		SEXP Y0,
		SEXP Eps
		){ 
    int na, np, i;
    double *xa, *ya, *xp, *yp;
    double x0, y0, eps;
    int *result;
    SEXP out;

    // protect arguments from garbage collector    
    PROTECT(P   = AS_LIST(P));
    PROTECT(A   = AS_LIST(A));
    PROTECT(X0  = AS_NUMERIC(X0));
    PROTECT(Y0  = AS_NUMERIC(Y0));
    PROTECT(Eps = AS_NUMERIC(Eps));

    // Get test point coordinates
    np = LENGTH(VECTOR_ELT(P, 0));
    xp = NUMERIC_POINTER(VECTOR_ELT(P, 0));
    yp = NUMERIC_POINTER(VECTOR_ELT(P, 1));

    // Get polygon coordinates
    na = LENGTH(VECTOR_ELT(A, 0));
    xa = NUMERIC_POINTER(VECTOR_ELT(A, 0));
    ya = NUMERIC_POINTER(VECTOR_ELT(A, 1));
    
    // Get scale parameters
    x0 = *(NUMERIC_POINTER(X0));
    y0 = *(NUMERIC_POINTER(Y0));
    eps = *(NUMERIC_POINTER(Eps));

    // copy and scale polygon
    ClipperLib::Path pathA;
    ScaleToPath(xa, ya, na, pathA, x0, y0, eps);

    // Allocate space for output
    PROTECT(out = NEW_INTEGER(np));
    result = INTEGER_POINTER(out);

    // handle one point at a time
    ClipperLib::IntPoint pti;
    for(i = 0; i < np; i++) {
      ScaleToPoint(xp[i], yp[i], pti, x0, y0, eps);
      result[i] = PointInPolygon(pti, pathA);
    }
    
    UNPROTECT(6);
    return(out);
  }
}
