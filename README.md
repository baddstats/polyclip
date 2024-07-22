polyclip
========

[![Travis-CI Build Status](https://travis-ci.org/baddstats/polyclip.png?branch=master)](https://travis-ci.org/baddstats/polyclip)
[![codecov.io](https://codecov.io/github/baddstats/polyclip/coverage.svg?branch=master)](https://codecov.io/github/baddstats/polyclip?branch=master)
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/polyclip)](http://cran.r-project.org/web/packages/polyclip)
[![Research software impact](http://depsy.org/api/package/cran/polyclip/badge.svg)](http://depsy.org/package/r/polyclip)

This repository holds the contributed R-package `polyclip`, which is
an R port of Angus Johnson's library 
[Clipper](http://angusj.com/delphi/clipper.php) for polygon clipping.

## Version of Clipper Library

This version of `polyclip` is derived from 
Clipper C++ library version `6.4.0 [r496]` which was obtained from the
[Sourceforge repository](https://sourceforge.net/projects/polyclipping)
(click `Code` then `Download snapshot`).
Minor changes have been made to the C++ code to satisfy the
requirements for R packages (namely, data type declarations must be portable,
and error messages must go through R's error handler). 

**Note:** If your system already includes the `polyclipping` library
(another derivative of `clipper`)
then **that version of the library will be used**.
That is, the R package `polyclip` will be compiled against
the executable library `polyclipping` on your system,
rather than using the bundled source code of `clipper 6.4.0`
that comes with the `polyclip` sources.

## Installation

The current official release of `polyclip` is available
on [CRAN](http://cran.r-project.org/web/packages/polyclip)
and can be downloaded and installed automatically
using the R command `install.packages`. 

The code in this repository is the development version,
which may be newer than the official release.
The easiest way to install the development version of `polyclip` 
from github is through the `remotes` package:

```R
require(remotes)
install_github('baddstats/polyclip')
```

If you don't have `remotes` installed you should first run

```R
install.packages('remotes')
```

## Bug reports 

Users of `polyclip` are encouraged to report bugs here 
(go to *issues* in the menu above, 
and press *new issue* to start a new bug report
or feature request).

## Making your own changes

Feel free to fork `polyclip`, make changes to the code,
and ask us to include them in the package by making a github *pull request*. 

