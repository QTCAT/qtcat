#include <Rcpp.h>
#include <iostream>
#include <fstream>
#include <string>

using namespace std;



// split string in to vector of strings
void split(const string &s, char delim, vector<string> &elems) {
    stringstream ss(s);
    string item;
    while (getline(ss, item, delim)) {
        elems.push_back(item);
    }
}

// [[Rcpp::export]]
Rcpp::List read_snpData(Rcpp::CharacterVector file, char sep, char quote, 
                        bool rowNames, Rcpp::CharacterVector na_str, 
                        int nrows) {
    string oneLine;
    vector<string> lineElements;
    const string fname = Rcpp::as<string>(file);
    const string naStr = Rcpp::as<string>(na_str);
    // Input file stream
    ifstream fileIn(fname.c_str());
    
    // Column names (Individual names)
    getline(fileIn, oneLine);
    // Delete qoutes if exist
    if (quote != ' ') {
        oneLine.erase(remove(oneLine.begin(), oneLine.end(), quote), 
                      oneLine.end());
    }
    // split string in to vector of strings
    split(oneLine, sep, lineElements);
    // Individual names
    Rcpp::CharacterVector indivNames = Rcpp::wrap(lineElements);
    
    // 
    vector<string> lociNames;
    vector<int> pos;
    set<char> allelesSet;
    vector<char> allelepos;
    string AA;
    string AB;
    string BA;
    string BB;
    vector<char> snpData;
    int col = 0;
    int row = 0;
    int posStart = 0;
    int dataStart = 2;
    if (rowNames) {
        ++ dataStart;
        ++ posStart;
    }
    bool unphased = true;
    while(getline(fileIn, oneLine)) {
        lineElements.clear();
        allelesSet.clear();
        // check if sep is part of first line
        if (col == 0 && oneLine.find(sep) == string::npos) {
            Rcpp::stop("The separator character 'sep' doesn't exist in line 1");
        }
        // Delete qoutes if exist
        if (quote != ' ') {
            oneLine.erase(remove(oneLine.begin(), oneLine.end(), quote), 
                          oneLine.end());
        }
        // split string in to vector of strings
        split(oneLine, sep, lineElements);
        // length of line
        if (col == 0) {
            row = lineElements.size();
        } else if (row != lineElements.size()){
            Rcpp::Rcerr << "Error: Length of line " << col +1 << " is " << 
                lineElements.size() << " instead of " << row <<  "\n" << endl;
            Rcpp::stop("Reading stopt");
        }
        // detect Nucleotides in this line (SNP)
        for (int i = dataStart; i < lineElements.size(); i++) {
            if (lineElements[i] != naStr) {
                allelesSet.insert(lineElements[i].begin(), 
                                  lineElements[i].end());
            }
        }
        vector<char> allele(allelesSet.begin(), allelesSet.end());
        // if more not two alleles skip
        if (allele.size() != 2) {
            Rcpp::Rcerr << "warning:\nline " << col +1 << 
                " has more or less than two alleles and therefore the line is skipt\n" 
                << endl;
            continue;
        }
        allelepos.push_back(allele[0]);
        allelepos.push_back(allele[1]);
        // if two alleles make new coding
        AA = allele[0]; AA += allele[0];
        AB = allele[0]; AB += allele[1];
        BA = allele[1]; BA += allele[0];
        BB = allele[1]; BB += allele[1];
        // raw coding of line (SNP) 
        for (int i = dataStart; i < lineElements.size(); ++i) {
            if (lineElements[i] == AA) {
                snpData.push_back(0x01);
            } else if (lineElements[i] == AB) {
                snpData.push_back(0x02);
            } else if (lineElements[i] == BA) {
                unphased = false;
                snpData.push_back(0x04);
            } else if (lineElements[i] == BB) {
                snpData.push_back(0x05);
            } else {
                snpData.push_back(0x00);
            }
        }
                // SNP names
        if (rowNames) {
            lociNames.push_back(lineElements[0]);
        }
        // genetic positions
         pos.push_back(atoi(lineElements[posStart].c_str()));
         pos.push_back(atoi(lineElements[posStart +1].c_str()));
        ++col;
        if (nrows > 0 && nrows == col) {
            break;
        }
    }
    // Rcpp conversioan and vector as matrix
    // position
    Rcpp::IntegerVector position = Rcpp::wrap(pos);
    position.attr("dim") = Rcpp::Dimension(2, col);
    
    // alleles
    Rcpp::CharacterVector alleles = Rcpp::wrap(allelepos);
    alleles.attr("dim") = Rcpp::Dimension(2, col);
    
    // snpData
    Rcpp::RawVector snpOutData(snpData.size());
    copy(snpData.begin(), snpData.end(), snpOutData.begin());
    if (unphased) {
        for (int i = 0; i < snpOutData.size(); i++) {
            if (snpOutData[i] == 0x02) {
                snpOutData[i] = 0x03;
            }
        }
    }
    row  = snpOutData.size() / col;
    snpOutData.attr("dim") = Rcpp::Dimension(row, col);
    
    // individual names
    indivNames.erase(indivNames.begin(), indivNames.begin() +2);
    
    // Results as List
    Rcpp::List  out = Rcpp::List::create(Rcpp::Named("snpData", snpOutData),
                                         Rcpp::Named("alleles", alleles),
                                         Rcpp::Named("position", position),
                                         Rcpp::Named("lociNames", lociNames),
                                         Rcpp::Named("indivNames", indivNames));
    return out;
}
