#include <Rcpp.h>

using namespace Rcpp;

// correlation
double cor(RawVector x, RawVector y) {
    int n = x.size();
    IntegerVector X(n), Y(n);
    double ex = 0, ey = 0, xt = 0, yt = 0,
           sxx = 0, syy = 0, sxy = 0;
    // setup variables and find the mean
    for (int i = 0; i < n; i++) {
        if (x[i] == 0x05)  {
            X[i] = 4;
        } else if (x[i] == 0x02 |
                   x[i] == 0x03 |
                   x[i] == 0x04) {
            X[i] = 2;
        }
        if (y[i] == 0x05)  {
            Y[i] = 4;
        } else if (y[i] == 0x02 |
                   y[i] == 0x03 |
                   y[i] == 0x04) {
            Y[i] = 2;
        }
        ex += X[i];
        ey += Y[i];
    }
    ex /= n;
    ey /= n;
    // correlation coefficent
    for (int i = 0; i < n; i++) {
        xt = X[i] - ex;
        yt = Y[i] - ey;
        sxx += xt * xt;
        syy += yt * yt;
        sxy += xt * yt;
    }
    double cor = sxy / (sqrt(sxx * syy) + 1e-16);
    return cor;
} // cor

//  1 x 1 correlation distance
// [[Rcpp::export]]
double corDist(RawVector x, RawVector y) {
   double dcor = 1 - std::fabs(cor(x, y));
   return dcor;
}

// distance from all marker
// [[Rcpp::export]]
NumericVector corDists(RawMatrix x) {
    int n = x.ncol();
    NumericVector dist(n * (n -1) / 2);
    int count = 0;
    for (int i = 0; i < n - 1; i++) {
        for (int j = i + 1; j < n; j++) {
            dist[count] = corDist(x(_, i), x(_, j));
            count ++;
        }
    }
   return dist;
} // corDists

//pre cluster for identical search
// [[Rcpp::export]]
List corPreIdenticals(RawMatrix x, const int step) {
    int n = x.ncol();
    IntegerVector clusters(n);
    std::map< int, std::vector<int> > preclust;
    std::vector<int> premedoInx;
    int k = 0;
    int i0, j0;
    if (n >= step * 2) {
        // medoids equally distributed over the search space
        bool smallDist = false;
        int startStep = (int)step * .5;
        double medodist = 0;
        premedoInx.push_back(startStep);
        for (int i = startStep+step; i < n; i += step) {
            for(int j = startStep; j < i; j += step) {
              medodist = corDist(x(_, i), x(_, j));
              if (medodist < .05) {
                  smallDist = true;
                  break;
              }
            }
            if (smallDist) {
              smallDist = false;
              continue;
            }
            premedoInx.push_back(i);
        }
        k = premedoInx.size();
    }
    if (k > 1) {
        // cluster all variables to the nearest medoid
        NumericVector predist(k);
        for (int i = 0; i < n; i ++) {
            for (int j = 0; j < k; j ++) {
                j0 = premedoInx[j];
                predist[j] = corDist(x(_, i), x(_, j0));
            }
            preclust[which_min(predist)].push_back(i);
        }
    } else {
        for (int i = 0; i < n; i ++) {
            preclust[0].push_back(i);
        }
        k=1;
    }
    return List::create(preclust);
} // preClustIdenticals

// find identicals
// [[Rcpp::export]]
List corIdenticals(RawMatrix x, IntegerVector clustIdx) {
    int nCluster = clustIdx.size();
    double dist = 0;
    IntegerVector clusters(nCluster);
    int clustCount = 1;
    std::vector<int> grp, medoInx;
    int gMedo, i0, j0;
    for (int i = 0; i < nCluster; i ++) {
        i0 = clustIdx[i];
        if (clusters[i] == 0) {
            clusters[i] = clustCount;
            grp.push_back(i);
            for (int j = i+1; j < nCluster; j ++) {
                j0 = clustIdx[j];
                if (clusters[j] == 0) {
                    dist = corDist(x(_, i0), x(_, j0));
                    if (dist <= 1e-7) {
                        clusters[j] = clustCount;
                        grp.push_back(j);
                    }
                }
            }
            gMedo = grp.size() / 2;
            medoInx.push_back(grp[gMedo]);
            clustCount ++;
            grp.clear();
        }
    }
    return List::create(clusters, medoInx);
} // identicals

// join data to one object
// [[Rcpp::export]]
List joinCorIdenticals(int n, List preclust, List ClustMedo) {
     IntegerVector clusters(n);
     std::vector<int> medoInx;
     List SubClustMedo;
     IntegerVector clustIdx, clust, medo;
     int count = 0, j0 = 0;
     for (int i = 0; i < preclust.size(); i ++) {
         clustIdx = preclust[i];
         SubClustMedo = ClustMedo[i];
         clust = SubClustMedo[0];
         medo = SubClustMedo[1];
         for (int j = 0; j < clustIdx.size(); j ++) {
             j0 = clustIdx[j];
             clusters[j0] = clust[j] + count;
         }
         for (int j = 0; j < medo.size(); j ++) {
             j0 = medo[j];
             medoInx.push_back(clustIdx[j0]);
         }
         count = medoInx.size();
     }
     return List::create(clusters, medoInx);
} // joinIdenticals

// clustering of IDB by correlation dictamce with CLARANS
// [[Rcpp::export]]
List corClarans(RawMatrix x, const int k, const int maxNeigbours) {
    RNGScope scope;
    Environment base("package:base");
    Function sample_int = base["sample.int"];
    int n = x.ncol();
    NumericMatrix medoDist(n, k);
    IntegerVector medoInx(k);
    double dist = 0;
    IntegerVector clusters(n);
    // starting medoid and clusters
    NumericVector pdist(n, 1.0);
    int m = 0;
    for (int c = 0; c < k; c ++) {
        // starting medoid
        m = as<int>(sample_int(n, 1, false, pdist))-1;
        medoInx[c] = m;
        for (int i= 0; i < n; i ++) {
            dist = corDist(x(_, m), x(_, i));
            medoDist(i, c) = dist;
        }
    }
    int clust = 0;
    double costs = 0, objective = 0;
    for (int i= 0; i < n; i ++) {
        // starting cluster membership
        clust = which_min(medoDist(i, _));
        clusters[i] = clust;
        costs += medoDist(i, clust);
    }
    objective = costs / n;
    IntegerVector ranVar(maxNeigbours);
    int i = 0, o = 0;
    double minDist = 0;
    NumericVector objectDist(n);
    NumericVector objectMinDist(n);
    double object_cost = 0;
    // iteration of inner clarans loop
    int iter = 0;
    while (iter < maxNeigbours) {
        if ( iter == 0) {
            ranVar = floor(runif(maxNeigbours, 0, n-1e-15));
        }
        // replacing medoids with objects if costs are smaller
        i = ranVar[iter];
        o = clusters[i];
        for (int j = 0; j < n; j ++) {
            dist = corDist(x(_, i), x(_, j));
            objectDist[j] = dist;
            minDist = dist;
            for (int c = 0; c < k; c ++) {
                if (c != o) {
                    minDist = std::min(minDist, medoDist(j, c));
                }
            }
            objectMinDist[j] = minDist;
        }
        object_cost = sum(objectMinDist);
        if (object_cost < costs) {
            // new medoids
            medoDist(_, o) = objectDist;
            medoInx[o] = i;
            //  new costs
            costs = object_cost;
            objective = costs / n;
            // new membership to cluster
            for (int j = 0; j < n; j ++) {
                clusters[j] = which_min(medoDist(j, _));
            }
            // restarting loop
            iter = 0;
        } else {
            iter ++;
        }
    }
    // Result as list
    return List::create(clusters+1, medoInx, objective);
} // clarans

// medoids of SNP clusters
// [[Rcpp::export]]
IntegerVector corMedoids(RawMatrix x, IntegerVector clusters) {
    int n = clusters.size();
    // map of vectors with cluster indexes
    std::map< int, std::vector<int> > indexMap;
    for (int i = 0; i < n; ++i) {
        indexMap[clusters[i]].push_back(i);
    }
    // for each cluster a dist matrix were the col with min dist is the medoid
    std::vector<int> clustInx;
    int nClust = 0;
    double dist = 0;
    int dist_j, dist_i, medoInx;
    IntegerVector medoids(n);
    for (std::map< int, std::vector<int> >::iterator it = indexMap.begin();
         it != indexMap.end(); ++ it) {
        clustInx = it->second;
        nClust = clustInx.size();
        NumericMatrix dist_mat(nClust, nClust);
        NumericVector dist_sums(nClust);
        for (int j = 0; j < nClust-1; ++ j) {
          dist_j = clustInx[j];
            for (int i = j; i < nClust; ++ i) {
                if ( i == j) {
                    dist_mat(j, i) = 0;
                } else {
                   dist_i = clustInx[i];
                   dist = corDist(x(_, dist_j), x(_, dist_i));
                   dist_mat(j, i) = dist;
                   dist_mat(i, j) = dist;
                }
            }
        }
        for (int i = 0; i < nClust; ++ i) {
            dist_sums[i] = sum(dist_mat(_, i));
        }
        medoInx = which_min(dist_sums);
        medoids[clustInx[medoInx]] = 1;
    }
    return medoids;
}
