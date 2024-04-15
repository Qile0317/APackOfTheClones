#include <Rcpp.h>
#include <unordered_map>
#include <string>
#include "math.h"

// assumes x > 1
// [[Rcpp::export]]
std::vector<std::vector<int>> rcppGetUniquePairsUpTo(int x, bool oneIndexed) {

    std::vector<std::vector<int>> uniquePairList (x * (x - 1) / 2);

    int index = 0;
    for (int i = oneIndexed; i < x - 1 + oneIndexed; i++)  {
        for (int j = i + 1; j < x + oneIndexed; j++) {
            uniquePairList[index] = {i, j};
            index++;
        }
    }

    return uniquePairList;
}

class AddOnlyOrderedStringIntMap {
protected:
    std::unordered_map<std::string, unsigned> keyIndex;
    std::vector<int> values;

public:
    AddOnlyOrderedStringIntMap() {
        keyIndex = std::unordered_map<std::string, unsigned>();
        values = std::vector<int>();
    }

    bool containsKey(std::string key) {
        return keyIndex.find(key) != keyIndex.end();
    }

    AddOnlyOrderedStringIntMap& AddIntAtKey(Rcpp::NumericVector namedNumeric) {

        if (namedNumeric.size() == 0) return *this;

        std::vector<std::string> keys = namedNumeric.names();

        for (int i = 0; i < (int) namedNumeric.size(); i++) {
            AddIntAtKey(keys[i], namedNumeric[i]);
        }

        return *this;
    }

    AddOnlyOrderedStringIntMap& AddIntAtKey(std::string key, int value) {

        if (containsKey(key)) {
            values[keyIndex[key]] += value;
            return *this;
        }

        keyIndex[key] = values.size();
        values.push_back(value);
        return *this;
    }

    std::vector<int> getRawValues() {
        return values;
    }

    std::vector<std::string> getKeys() {
        std::vector<std::string> orderedNames (values.size());
        for (auto keyValPair : keyIndex) {
            orderedNames[keyValPair.second] = keyValPair.first;
        }
        return orderedNames;
    }

    Rcpp::NumericVector asNamedNumeric() {
        Rcpp::NumericVector output = Rcpp::wrap(getRawValues());
        output.names() = getKeys();
        return output;
    }

};

// input should be a list of named numeric that represent tables. some of which may be empty
// [[Rcpp::export]]
Rcpp::NumericVector rcppUnionListOfNamedNumericsHelper(Rcpp::List listOfNamedNumerics) {

    AddOnlyOrderedStringIntMap unionedNamedNumericMap;

    for (int i = 0; i < (int) listOfNamedNumerics.size(); i++) {
        Rcpp::NumericVector currNamedNumeric = listOfNamedNumerics[i];
        unionedNamedNumericMap.AddIntAtKey(currNamedNumeric);
    }

    return unionedNamedNumericMap.asNamedNumeric();
}
