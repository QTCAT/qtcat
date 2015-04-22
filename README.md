# Quantitative Trait Cluster Association Test

The package is built and tested under continuous integration. Current status for the main branch is:
Linux 
[![Travis-CI Build Status](https://travis-ci.org/QTCAT/qtcat.png?branch=master)]
(https://travis-ci.org/QTCAT/qtcat), and Windows
[![AppVeyor Build Status]
  (https://ci.appveyor.com/api/projects/status/github/jrklasen/qtcat?branch=master&svg=true)]
(https://ci.appveyor.com/project/jrklasen/qtcat).
The licence is: 
[![License]
  (http://img.shields.io/badge/license-GPL%20%28%3E=%202%29-brightgreen.svg?style=flat)]
(http://www.gnu.org/licenses/gpl-2.0.html)

## Description:

A association mapping method which jointly analysis all SNPs at ones and at the 
same time accounts for the correlation among them. This makes population 
structure correction unnecessary and therefore increases power compared to 
classical methods like the mixed model.

## Install:

The package can be installed from an R console via `devtools`. If you haven't 
yet `devtools` installed, you have to do so first.

```R
# install.packages("devtools")
devtools::install_github("jrklasen/hit")
devtools::install_github("QTCAT/qtcat")
```
