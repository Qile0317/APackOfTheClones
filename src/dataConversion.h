#pragma once

#include <Rcpp.h>
#include <vector>
#include <string>
#include <unordered_map>

// assumes list names are unique
// assumes all elements of the list are single integers
std::unordered_map<std::string, int> namedIntListToHashMap(Rcpp::List l) {
    int n = l.size();
    if (n == 0) return std::unordered_map<std::string, int>();

    std::unordered_map<std::string, int> outputHashMap (n);
    Rcpp::CharacterVector listNames = l.names();
   
    for (int i = 0; i < n; i++) {
        outputHashMap[listNames[i]] = Rcpp::as<int>(l[i]);
    }

    return outputHashMap;
}

// for each element in the list, same assumptions as function above
std::vector<std::unordered_map<std::string, int>> listOfNamedIntListToHashMapVector(
    Rcpp::List l
) {
    int n = l.size();
    if (n == 0) return std::unordered_map<std::string, int>();

    std::vector<std::unordered_map<std::string, int>> outputHashMapVector (n);

    for (int i = 0; i < n; i++) {
        outputHashMapVector[i] = namedIntListToHashMap(
            Rcpp::as<Rcpp::List>(l[i])
        );
    }

    return outputHashMapVector;
}
