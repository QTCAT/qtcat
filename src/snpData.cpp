#include <Rcpp.h>
#include <algorithm>
#include <string>

using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix design(RawMatrix x, IntegerVector inx1, IntegerVector inx2) {
     size_t n = x.nrow();
     size_t p = inx1.size();
     NumericMatrix out(n, p);
     for (size_t i = 0; i < p; ++i) {
        if (inx1[i] == inx2[i]) {
            for (size_t j = 0; j < n; ++j) {
                if (x(j, inx1[i]) == 0x01) {
                    out(j, i) = 0;
                } else if (x(j, inx1[i]) == 0x05) {
                    out(j, i) = 1;
                }else if (x(j, inx1[i]) == 0x02 || 
                           x(j, inx1[i]) == 0x03 ||
                           x(j, inx1[i]) == 0x04) {
                    out(j, i) = .5;
                }
            }
        } else {
            for (size_t j = 0; j < n; ++j) {
                if (x(j, inx1[i]) == 0x01 || 
                    x(j, inx2[i]) == 0x01) {
                    out(j, i) = 0;
                } else if (x(j, inx1[i]) == 0x05 && 
                           x(j, inx2[i]) == 0x05) {
                    out(j, i) = 1;
                } else if ((x(j, inx1[i]) == 0x02 || 
                            x(j, inx1[i]) == 0x03 ||
                            x(j, inx1[i]) == 0x04) &&
                           x(j, inx2[i]) == 0x05) {
                    out(j, i) = .5;
                } else if (x(j, inx1[i]) == 0x05 &&
                           (x(j, inx2[i]) == 0x02 || 
                            x(j, inx2[i]) == 0x03 ||
                            x(j, inx2[i]) == 0x04)) {
                    out(j, i) = .5;
                } else if ((x(j, inx1[i]) == 0x02 ||
                            x(j, inx1[i]) == 0x03 ||
                            x(j, inx1[i]) == 0x04) &&
                           (x(j, inx2[i]) == 0x02 ||
                            x(j, inx2[i]) == 0x03 ||
                            x(j, inx2[i]) == 0x04)) {
                    out(j, i) = .25;
                }
            }
        }
    }
    return out;
}

// [[Rcpp::export]]
List freqs2(RawMatrix x) {
    size_t p = x.ncol();
    size_t n = x.nrow();
    double aCount = 0;
    double af = 0;
    double hetCount = 0;
    NumericVector maf(p), hetf(p);
    for (size_t i = 0; i < p; i ++) {
        for (size_t j = 0; j < n; j ++) {
            if (x(j, i) == 0x01) {
                aCount ++;
            } else if (x(j, i) == 0x02 || 
                       x(j, i) == 0x03 || 
                       x(j, i) == 0x04) {
                aCount += .5;
                hetCount ++;
            }
        }
        af = aCount/n;
        if (af <= .5) {
            maf[i] = af;
        } else {
            maf[i] = 1-af;
        }
        hetf[i] = hetCount/n;
        aCount = 0;
        hetCount = 0;
    }
    return List::create(maf, hetf);
}

// [[Rcpp::export]]
NumericVector freq1(RawMatrix x) {
    size_t p = x.ncol();
    size_t n = x.nrow();
    double hetCount = 0;
    NumericVector hetf(n);
    for (size_t i = 0; i < n; i ++) {
        for (size_t j = 0; j < p; j ++) {
            if (x(i, j) == 0x02 || 
                x(i, j) == 0x03 || 
                x(i, j) == 0x04) {
                hetCount ++;
            }
        }
        hetf[i] = hetCount/p;
        hetCount = 0;
    }
    return hetf;
}
