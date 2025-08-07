#ifndef CONTROLPICKSTRATEGY_H
#define CONTROLPICKSTRATEGY_H

#include <string>

#include "ControlPickCriteria.h"

class ControlPicker;
class ControlPickCriteria;

class ControlPickStrategy {    
public:
    ControlPickStrategy() = default;
    /*
      pickControl doit renvoyer soit un contrôle viable de grille, soit une raison pour laquelle il n'a pas réussi à en trouver un
     */
    virtual OptionalCu pickControl(const OptionalCu &previousStrategyResult,
                                   ControlPickCriteria &criteria) = 0;
    virtual ~ControlPickStrategy() = default;

    virtual const std::string &getName() const = 0;
};

class UserPickStrategy : public ControlPickStrategy {    
public:
    UserPickStrategy() = default;
    virtual ~UserPickStrategy() = default;

    const std::string &getName() const override final;
    void setName(const std::string &name);
private:
    std::string name;
};

#endif /* CONTROLPICKSTRATEGY_H */
