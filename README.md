
<a name="logo"/>
<div align="center">
<a href="https://qtcat.github.io/" target="_blank">
<img src="https://github.com/QTCAT/qtcat_logo/blob/master/logo/QTCAT_l.png?raw=true" width="500"></img>
</a>
</div>

# Quantitative Trait Cluster Association Test

The current built and test status for Linux 
[![Travis-CI Build Status](https://travis-ci.org/QTCAT/qtcat.png?branch=master)]
(https://travis-ci.org/QTCAT/qtcat) and for Windows 
[![AppVeyor Build status]
  (https://ci.appveyor.com/api/projects/status/hx1pvqer9flugwew?branch=master&svg=true)]
(https://ci.appveyor.com/project/jrklasen/qtcat).


## Description:
An association mapping method which jointly analyses all SNPs at once and at the 
same time accounts for the correlation among them. This makes correction for 
population structure unnecessary and therefore increases power compared to 
classical methods like the mixed model.

## Install:

The package can be installed from an R console via `devtools`. If you haven't 
yet installed `devtools`, you have to do so first. In addition `hit` is a 
dependency, which has to be installed first.

```R
# install.packages("devtools")
# devtools::install_github("jrklasen/hit")
devtools::install_github("QTCAT/qtcat")
```

[![License]
  (http://img.shields.io/badge/license-GPL%20%28%3E=%202%29-brightgreen.svg?style=flat)]
(http://www.gnu.org/licenses/gpl-2.0.html)
&copy; 2015 JR Klasen
