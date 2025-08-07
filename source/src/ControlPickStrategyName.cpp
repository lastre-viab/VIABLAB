#include "../include/ControlPickStrategyName.h"

ControlPickStrategyName ControlPickStrategyName::create(const std::string &str) {

    ControlPickStrategyName res;
    
    STRATEGY_VALUES(STRING_TO_ENUM) {
        res = ControlPickStrategyName(str);
    }
    return res;
}

PredefinedStrategyName ControlPickStrategyName::getPredefinedStrategyName() const {
    return name;
}

const std::string & ControlPickStrategyName::getUserStrategyName() const {
    return userStrategyName;
}

bool ControlPickStrategyName::operator==(PredefinedStrategyName n) const {
    return n == name;
} 

ControlPickStrategyName::ControlPickStrategyName(PredefinedStrategyName predefinedName) :
    userStrategyName(""),
    name(predefinedName) {
    if (predefinedName == USER_STRATEGY) {
        throw std::invalid_argument("ControlPickStrategyName defined as custom in constructor that expects a predefined strategy");
    }
}

ControlPickStrategyName::ControlPickStrategyName(const std::string &str) :
    userStrategyName(str),
    name(USER_STRATEGY) {}

std::ostream &operator<<(std::ostream &os, const ControlPickStrategyName &name) {
    if (name == USER_STRATEGY) {
        return os << name.getUserStrategyName();
    }
    else {
        return os << name.getPredefinedStrategyName();
    }
}

using Translator = ControlPickStrategyNameTranslator;

boost::optional<Translator::external_type> Translator::get_value(const internal_type &str) {
    boost::optional<external_type> res = ControlPickStrategyName::create(str);
    return res;
}

boost::optional<Translator::internal_type> Translator::put_value(external_type value) {
    boost::optional<internal_type> res;
    STRATEGY_VALUES(ENUM_TO_STRING) {
        res = value.getUserStrategyName();
    }
    return res;
}  
