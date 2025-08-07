#ifndef TYCHASTICCONTROLPICKSTRATEGY_H
#define TYCHASTICCONTROLPICKSTRATEGY_H

#include <string>

#include "TychasticControlPickCriteria.h"

class TychasticControlPicker;
class TychasticControlPickCriteria;

class TychasticControlPickStrategy {    
public:
    TychasticControlPickStrategy() = default;
    /*
      pickControl doit renvoyer soit un contrôle viable de grille, soit une raison pour laquelle il n'a pas réussi à en trouver un
     */
    virtual OptionalCu pickControl(const OptionalCu &previousStrategyResult,
                                   TychasticControlPickCriteria &criteria) = 0;
    virtual ~TychasticControlPickStrategy() = default;

    virtual const std::string &getName() const = 0;
};

class TychasticUserPickStrategy : public TychasticControlPickStrategy {    
public:
    TychasticUserPickStrategy() = default;
    virtual ~TychasticUserPickStrategy() = default;

    const std::string &getName() const override final;
    void setName(const std::string &name);
private:
    std::string name;
};

#endif /* TYCHASTICCONTROLPICKSTRATEGY_H */
