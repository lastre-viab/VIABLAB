#ifndef MODELPARAMETERS_H
#define MODELPARAMETERS_H

#include <cstdlib>
#include <string>
#include <unordered_map>
#include <vector>
#include <sstream>
#include "viablab_export.h"

// En lowerCase par cohérence avec les autres paramètres
class VIABLAB_LIBRARY_EXPORT modelParams {
public:
    std::string &operator[](const std::string &str);
    const std::string &operator[](const std::string &str) const;
    const std::string &at(const std::string &str) const;
    
    void addToList(const std::string &key, const std::string &value);
    
    double getDouble(const std::string &str) const;
    int getInt(const std::string &str) const;
    const std::string &getString(const std::string &str) const;
    bool getBool(const std::string &str) const;
    
    void copyDoubles(const std::string &str, double *array) const;
    void copyInts(const std::string &str, int *array) const;
    void copyStrings(const std::string &str, std::string *array) const;
    void copyBools(const std::string &str, bool *array) const;
private:
    std::unordered_map<std::string, std::string> modelParameters;
    std::unordered_map<std::string, std::vector<std::string>> arrayModelParameters;
};

#endif /* MODELPARAMETERS_H */
