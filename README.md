<a name="logo"/>
<div align="center">
<a href="https://github.com/QTCAT" target="_blank">
<img src="https://github.com/QTCAT/qtcat_logo/blob/master/logo/QTCAT_l.png?raw=true" width="600"></img>
</a>
</div>

--------------------------------------------------------------------------------

# Quantitative Trait Cluster Association Test
The current built and test status for Linux (Mac)
[![Build Status](https://travis-ci.org/QTCAT/qtcat.svg?branch=master)](https://travis-ci.org/QTCAT/qtcat)
and for Windows 
[![Build status](https://ci.appveyor.com/api/projects/status/hx1pvqer9flugwew/branch/master?svg=true)](https://ci.appveyor.com/project/jrklasen/qtcat/branch/master)


## Description:
All SNPs are jointly associated to the phenotype and at the same time correlation among 
them is considered. Thus, correction for population structure becomes unnecessary, which in 
many cases results in a power advantages compared to other methods.

**Klasen, J. R. et al. (2016)**. *A multi-marker association method for genome-wide 
association studies without the need for population structure correction*. Nature 
Communications. [Paper](http://www.nature.com/articles/ncomms13299)

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
