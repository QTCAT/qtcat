# Quantitative Trait Cluster Association Test

The current built and test status under Linux 
[![Travis-CI Build Status](https://travis-ci.org/QTCAT/qtcat.png?branch=master)]
(https://travis-ci.org/QTCAT/qtcat) and Windows 
[![AppVeyor Build status]
  (https://ci.appveyor.com/api/projects/status/hx1pvqer9flugwew?branch=master&svg=true)]
(https://ci.appveyor.com/project/jrklasen/qtcat).

## Description:

A association mapping method which jointly analyses all SNPs at ones and at the 
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

[![License]
  (http://img.shields.io/badge/license-GPL%20%28%3E=%202%29-brightgreen.svg?style=flat)]
(http://www.gnu.org/licenses/gpl-2.0.html)
&copy; 2015 JR Klasen
