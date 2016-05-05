#include <Rcpp.h>
#include <iostream>
#include <fstream>
#include <string>

using namespace std;

void split(const string &s, char delim, vector<string> &elems);

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
        oneLine.erase(remove(oneLine.begin(), oneLine.end(), quote), oneLine.end());
    }
    // check if sep is part of first line
    if (oneLine.find(sep) == string::npos) {
        Rcpp::stop("In line 1 the separator character 'sep' doesn't exist");
    }
    // split string in to vector of strings
    split(oneLine, sep, lineElements);
    unsigned int row = lineElements.size();
    // Individual names
    Rcpp::CharacterVector indivNames = Rcpp::wrap(lineElements);
    //
    vector<string> lociNames;
    vector<string> chr;
    vector<int> pos;
    set<char> allelesSet;
    vector<char> allelepos;
    string AA;
    string AB;
    string BA;
    string BB;
    vector<char> snpData;
    int col = 0;
    int posStart = 0;
    int dataStart = 2;
    if (rowNames) {
        ++ dataStart;
        ++ posStart;
    }
    while(getline(fileIn, oneLine)) {
        lineElements.clear();
        allelesSet.clear();
        // Delete qoutes if exist
        if (quote != ' ') {
            oneLine.erase(remove(oneLine.begin(), oneLine.end(), quote), oneLine.end());
        }
        // split string in to vector of strings
        split(oneLine, sep, lineElements);
        // length of line
        if (row != lineElements.size()){
            Rcpp::stop("Error: Length of line ", col + 2, " is ",
                lineElements.size(), " instead of ", row);
        }
        // detect Nucleotides in this line (SNP)
        for (unsigned int i = dataStart; i < lineElements.size(); ++ i) {
            if (lineElements[i] != naStr) {
                allelesSet.insert(lineElements[i].begin(), lineElements[i].end());
            }
        }
        vector<char> allele(allelesSet.begin(), allelesSet.end());
        // if more not two alleles skip
        if (allele.size() != 2) {
            Rcpp::Rcerr << "warning: Line " << col + 2 <<
                " has more or less than two alleles and therefore the line is skipt"
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
        for (unsigned int i = dataStart; i < lineElements.size(); ++ i) {
            if (lineElements[i] == AA) {
                snpData.push_back(0x01);
            } else if (lineElements[i] == BB) {
                snpData.push_back(0x03);
            } else if ((lineElements[i] == AB) |
                       (lineElements[i] == BA)) {
                snpData.push_back(0x02);
            } else {
                snpData.push_back(0x00);
            }
        }
        // SNP names
        if (rowNames) {
            lociNames.push_back(lineElements[0]);
        }
        // genetic positions
         chr.push_back(lineElements[posStart]);
         pos.push_back(atoi(lineElements[posStart + 1].c_str()));
        ++ col;
        if (nrows > 0 && nrows <= col + 1) {
            break;
        }
    }
    if (snpData.size() == 0)
        Rcpp::stop("No valid SNP found");
    // Rcpp conversioan and vector as matrix
    // alleles
    Rcpp::CharacterVector alleles = Rcpp::wrap(allelepos);
    alleles.attr("dim") = Rcpp::Dimension(2, col);
    // snpData
    Rcpp::RawVector snpOutData(snpData.size());
    copy(snpData.begin(), snpData.end(), snpOutData.begin());
    row  = snpOutData.size() / col;
    snpOutData.attr("dim") = Rcpp::Dimension(row, col);
    // individual names
    indivNames.erase(indivNames.begin(), indivNames.begin() + dataStart);
    // Results as List
    Rcpp::List  out = Rcpp::List::create(Rcpp::Named("snpData", snpOutData),
                                         Rcpp::Named("alleles", alleles),
                                         Rcpp::Named("chr", chr),
                                         Rcpp::Named("pos", pos),
                                         Rcpp::Named("lociNames", lociNames),
                                         Rcpp::Named("indivNames", indivNames));
    return out;
}

// split string in to vector of strings
void split(const string &s, char delim, vector<string> &elems) {
    stringstream ss(s);
    string item;
    while (getline(ss, item, delim)) {
        elems.push_back(item);
    }
}
