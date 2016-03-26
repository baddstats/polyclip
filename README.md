polyclip
========

[![Travis-CI Build Status](https://travis-ci.org/baddstats/polyclip.png?branch=master)](https://travis-ci.org/baddstats/polyclip)
[![codecov.io](https://codecov.io/github/baddstats/polyclip/coverage.svg?branch=master)](https://codecov.io/github/baddstats/polyclip?branch=master)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/polyclip)](http://cran.r-project.org/web/packages/polyclip)
[![Research software impact](http://depsy.org/api/package/cran/polyclip/badge.svg)](http://depsy.org/package/r/polyclip)

This repository holds the contributed R-package `polyclip`, which is
an R port of Angus Johnson's library 
[Clipper](http://angusj.com/delphi/clipper.php) for polygon clipping.

## Version

The code in this repository is a development version, which may be
newer than the official release of `polyclip` on 
[CRAN](http://cran.r-project.org).

This version of `polyclip` is derived from 
Clipper C++ library version 6.4.0 which was obtained from the
[Sourceforge repository](https://sourceforge.net/projects/polyclipping).
Minor changes have been made to the C++ code to satisfy the
requirements for R packages (namely, data type declarations must be portable,
and error messages must go through R's error handler). 

## Installation

The easiest way to install the development version of `polyclip` 
from github is through the `devtools` package:

```R
require(devtools)
install_github('baddstats/polyclip')
```

If you don't have `devtools` installed you should first run

```R
install.packages('devtools')
```

## Bug reports 

Users of `polyclip` are encouraged to report bugs here 
(go to *issues* in the menu above, 
and press *new issue* to start a new bug report
or feature request).

## Making your own changes

Feel free to fork `polyclip`, make changes to the code,
and ask us to include them in the package by making a github *pull request*. 

