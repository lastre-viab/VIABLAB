#include "../include/ModelParameters.h"

double modelParams::getDouble(const std::string &str) const {
    return atof(modelParameters.at(str).c_str());
}

int modelParams::getInt(const std::string &str) const {
    return atoi(modelParameters.at(str).c_str());
}

bool modelParams::getBool(const std::string &str) const {
    bool b;
    std::istringstream(modelParameters.at(str)) >> std::boolalpha >> b;
    return b;
}
const std::string &modelParams::getString(const std::string &str) const {
    return modelParameters.at(str);
}

const std::string &modelParams::operator[](const std::string &str) const {
    return modelParameters.at(str);
}

std::string &modelParams::operator[](const std::string &str) {
    return modelParameters[str];
}

const std::string &modelParams::at(const std::string &str) const {
    return modelParameters.at(str);
}

void modelParams::addToList(const std::string &key, const std::string &value) {
    using vs = std::vector<std::string>;
    using uss = std::unordered_map<std::string, vs>;

    uss::iterator keyIt = arrayModelParameters.find(key);
    if (keyIt == arrayModelParameters.end()) {
        vs newParamsList{};
        newParamsList.push_back(value);
        arrayModelParameters[key] = newParamsList;
    }
    else {
        keyIt->second.push_back(value);
    }
}

void modelParams::copyDoubles(const std::string &key, double *array) const {
    int i = 0;
    for (const std::string &str : arrayModelParameters.at(key)) {
        array[i] = atof(str.c_str());
    }
}

void modelParams::copyInts(const std::string &key, int *array) const {
    int i = 0;
    for (const std::string &str : arrayModelParameters.at(key)) {
        array[i] = atoi(str.c_str());
    }
}

void modelParams::copyStrings(const std::string &key, std::string *array) const {
    int i = 0;
    for (const std::string &str : arrayModelParameters.at(key)) {
        array[i] = str;
    }
}

void modelParams::copyBools(const std::string &key, bool *array) const {
    int i = 0;
    for (const std::string &str : arrayModelParameters.at(key)) {
        bool b;
        std::istringstream(modelParameters.at(str)) >> std::boolalpha >> b;
        array[i] = b;
    }
}
