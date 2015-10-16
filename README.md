<a name="logo"/>
<div align="center">
<a href="https://github.com/QTCAT" target="_blank">
<img src="https://github.com/QTCAT/qtcat_logo/blob/master/logo/QTCAT_l.png?raw=true" width="600"></img>
</a>
</div>

--------------------------------------------------------------------------------

# Quantitative Trait Cluster Association Test
The current built and test status for Linux (Mac)
[![Travis-CI Build Status](https://travis-ci.org/QTCAT/qtcat.png?branch=master)](https://travis-ci.org/QTCAT/qtcat) and for Windows 
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/hx1pvqer9flugwew/branch/master?svg=true)](https://ci.appveyor.com/project/jrklasen/qtcat)
.

## Description:
An association mapping method which jointly analyses all SNPs at once and at 
the same time accounts for the correlation between them. This makes correction 
for population structure unnecessary and therefore increases power compared to 
classical methods like the mixed model.

## Install:
The package can be installed from an R console via [`devtools`](https://github.com/hadley/devtools#updating-to-the-latest-version-of-devtools)
. If you haven't yet installed `devtools` please do so first. One of the R 
package dependency of `qtcat` is `hit` which is currently only available from 
github and has therefore likewise to be installed before `qtcat`s installation.

```R
# install.packages("devtools")
# devtools::install_github("QTCAT/hit")
devtools::install_github("QTCAT/qtcat")
```

## Example:
The `qtcatQtc`-function example gives an overview of the functionality of 
the package and can be accessed once the package is loaded (please be aware, 
that running the example can take a few seconds).

```R
library(qtcat)
example("qtcatQtc")
```

--------------------------------------------------------------------------------

[![License](http://img.shields.io/badge/license-GPL%20%28%3E=%202%29-brightgreen.svg?style=flat)](http://www.gnu.org/licenses/gpl-2.0.html)
&copy; 2015 JR Klasen
