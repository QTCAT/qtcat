#include <Rcpp.h>
#include <algorithm>
#include <string>

using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix design(RawMatrix x) {
    size_t n = x.nrow();
    size_t p = x.ncol();
    NumericMatrix out(n, p);
    for (size_t i = 0; i < p; i ++) {
        for (size_t j = 0; j < n; j ++) {
            if (x(j, i) == 0x01) {
                out(j, i) = 0;
            } else if (x(j, i) == 0x03) {
                out(j, i) = 1;
            } else if (x(j, i) == 0x02) {
              out(j, i) = .5;
            }
        }
    }
    return out;
}


// [[Rcpp::export]]
NumericVector afreq(RawMatrix x, bool maf) {
    size_t p = x.ncol();
    size_t n = x.nrow();
    size_t n_na = 0;
    double a = 0;
    double af = 0;
    NumericVector allfreq(p);
    for (size_t i = 0; i < p; i ++) {
        for (size_t j = 0; j < n; j ++) {
            if (x(j, i) == 0x00) {
                n_na ++;
            } else if (x(j, i) == 0x01) {
                a ++;
            } else if (x(j, i) == 0x02) {
                a += .5;
            }
        }
        af = a / (n - n_na);
        if (maf) {
            if (af <= .5) {
                allfreq[i] = af;
            } else {
                allfreq[i] = 1 - af;
            }
        } else {
            allfreq[i] = af;
        }
        n_na = 0;
        a = 0;
    }
    return allfreq;
}


// [[Rcpp::export]]
NumericVector hetfreq(RawMatrix x, int dim) {
    size_t k, l;
    if (dim == 1) {
        k = x.nrow();
        l = x.ncol();
    } else {
        k = x.ncol();
        l = x.nrow();
    }
    size_t l_na = 0;
    double het = 0;
    NumericVector hetf(k);
    for (size_t i = 0; i < k; i ++) {
        for (size_t j = 0; j < l; j ++) {
            if (dim == 1) {
                if (x(i, j) == 0x00) {
                    l_na ++;
                } else if (x(i, j) == 0x02) {
                    het ++;
                }
            } else {
                if (x(j, i) == 0x00) {
                    l_na ++;
                } else if (x(j, i) == 0x02) {
                    het ++;
                }
            }
        }
        hetf[i] = het / (l - l_na);
        het = 0;
        l_na = 0;
    }
    return hetf;
}


// [[Rcpp::export]]
NumericVector nafreq(RawMatrix x, int dim) {
  size_t k, l;
    if (dim == 1) {
        k = x.nrow();
        l = x.ncol();
    } else {
        k = x.ncol();
        l = x.nrow();
    }
    double na = 0;
    NumericVector naf(k);
    for (size_t i = 0; i < k; i ++) {
        for (size_t j = 0; j < l; j ++) {
            if (dim == 1) {
                if (x(i, j) == 0x00)
                    na ++;
            } else {
                if (x(j, i) == 0x00)
                    na ++;
            }
        }
        naf[i] = na / l;
        na = 0;
    }
    return naf;
}
