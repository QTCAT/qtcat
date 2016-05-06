<a name="logo"/>
<div align="center">
<a href="https://github.com/QTCAT" target="_blank">
<img src="https://github.com/QTCAT/qtcat_logo/blob/master/logo/QTCAT_l.png?raw=true" width="600"></img>
</a>
</div>

--------------------------------------------------------------------------------

# Quantitative Trait Cluster Association Test
The current built and test status for Linux (Mac)
[![Build Status](https://travis-ci.org/QTCAT/qtcat.svg)](https://travis-ci.org/QTCAT/qtcat)
and for Windows 
[![Build status](https://ci.appveyor.com/api/projects/status/hx1pvqer9flugwew/branch/master?svg=true)](https://ci.appveyor.com/project/jrklasen/qtcat/branch/master)


## Description:
An association mapping method which jointly analyses all SNPs at once and at 
the same time accounts for the correlation between them. This makes correction 
for population structure unnecessary and therefore increases power compared to 
classical methods like the mixed model.

## Installation:
The package can be installed from an R console via [`devtools`](https://github.com/hadley/devtools#updating-to-the-latest-version-of-devtools)
(If you haven't yet installed `devtools` please do so first). 

```R
# install.packages("devtools")
devtools::install_github("QTCAT/qtcat") 

```

## Example:
The `qtcatQtc`-function example gives an overview of the functionality of 
the package and can be accessed once the package is loaded.

```R
library(qtcat)
example(qtcatQtc, run.dontrun = TRUE)

```

There is also a [Arabidopsis example](https://github.com/QTCAT/qtcat.data) available.

--------------------------------------------------------------------------------

[![License](https://img.shields.io/badge/license-GPL%20%28%3E=%202%29-brightgreen.svg)](https://www.gnu.org/licenses/gpl-2.0.html)
&copy; 2015 JR Klasen
