#include <numeric>

#include "../include/ControlPicker.h"
#include "../include/ControlPickStrategies.h"
#include "../include/ControlPicker.h"
#include "../include/TrajectoryHelpers.h"
#include "../include/TrajectoryPointsSimulation.h"
#include "../include/TrajectorySimulation.h"

OptionalCu HeavyPickStrategy::pickControl(const OptionalCu &opt, ControlPickCriteria &criteria) {
        
    if (opt.isLeft()) {return opt;}
    
    using ull = unsigned long long int;

    double rho = criteria.getTimeStep();
    const Grid *grid = criteria.getGrid();
    const int dim = criteria.getDim();    
    const double *xCoordsDouble = criteria.getCurrentPos();
    
    Trajectory *traj = criteria.getTrajectory();
    TrajectoryPoints *trajDiscrete = criteria.getDiscreteTrajectory();
    
    double **controlCoords = criteria.getControlCoords();  
    int lastControlIndex = traj->getLastControlIndex(initCu);
    const double *lastControl = controlCoords[lastControlIndex];

    double imageVect[dim];
    pickedControl pickedControl{rho, lastControlIndex};
    
    criteria.findViableTimeStep(pickedControl);

    // Si, pour un pas de temps ayant décrû jusqu'à ce qu'il n'y ait pas de différence
    // de position discrète aucun pas de temps ne permet d'atteindre une position viable,
    // on peut supposer que le contrôle lourd ne marche plus
    if (criteria.isSameGridPosition(pickedControl)) {
        spdlog::warn("Unsatisfied heavy pick strategy with control n°{}", lastControlIndex);
        return OptionalCu(UNSATISFIED_STRATEGY);
    }
    else {
        spdlog::info("Heavy control found");                
        return OptionalCu(pickedControl);
    }    
}

OptionalCu FirstPickStrategyBitSet::pickControl(const OptionalCu &opt, ControlPickCriteria &criteria) {

    if (opt.isLeft()) {return opt;}
    
    const Grid *grid = criteria.getGrid();
    // Conversion du vecteur vers un array
    const double *xCoordsDouble = criteria.getCurrentPos();
    unsigned long long int currentPos = criteria.getCurrentPosGridIndex();
    const double time = criteria.getTime();
    double rho = criteria.getTimeStep();
    const int dim = criteria.getDim();

    double **controlCoords = criteria.getControlCoords();

    double imageVect[dim];
    
    int usedCu;

    const unsigned long long viabSuccessor = viabiHelper->findViableDiscreteSuccessor(currentPos, time, rho, usedCu);

    double discreteRho = rho;
    
    bool succes = false;
    if (viabSuccessor > grid->nbTotalPoints) {
        spdlog::error("Optimal successor not found with first viable strategy");
        return OptionalCu(UNSATISFIED_STRATEGY);
    }
    else if (const int realCu = viabiHelper->findViabControl_bis(xCoordsDouble, viabSuccessor, rho, 25, 0.5, imageVect, succes);
             succes) {
        spdlog::info("First strategy control found");
        pickedControl c {rho, realCu};
        return OptionalCu(c);
    }
    else {
        spdlog::warn("No viable real control found to match with discrete position with first viable strategy. Returning grid control n°{}", usedCu);
        
        pickedControl c {discreteRho, usedCu};
        return OptionalCu(c);
    }
}

OptionalCu FirstOnlyPickStrategy::pickControl(
    const OptionalCu &opt, ControlPickCriteria &criteria) {

    if (opt.isLeft()) {return opt;}
    
    const Grid *grid = criteria.getGrid();
    // Conversion du vecteur vers un array
    const double rho = criteria.getTimeStep();
    
    double **controlCoords = criteria.getControlCoords();
    const int *preferedControlIndexes = criteria.sortPreferedControlIndexes();
    const double *control = controlCoords[preferedControlIndexes[0]];

    const int dim = criteria.getDim();
    double imageVect[dim];

    pickedControl c {rho, preferedControlIndexes[0]};
    
    bool viable = criteria.isViableControl(c, imageVect);
    unsigned long long viabSuccessor;

    bool succes = false;
    
    if (!viable ||
        (viabSuccessor = grid->getNearestPointInSet(imageVect)) >= grid->nbTotalPoints) {
        spdlog::warn("Optimal successor not found with first viable strategy");
        return OptionalCu(UNSATISFIED_STRATEGY);
    }
    else {
        spdlog::info("First only strategy control found");
        return OptionalCu(c);
    }    
}

OptionalCu ShuffleStrategy::pickControl(const OptionalCu &opt, ControlPickCriteria &criteria) {
    criteria.setIndexSorter(shuffleControlIndexes);
    return opt;
}

OptionalCu SortStrategy::pickControl(const OptionalCu &opt, ControlPickCriteria &criteria) {
    criteria.setIndexSorter(sortControlIndexesByWeight);
    return opt;
}

OptionalCu ResetOrderStrategy::pickControl(const OptionalCu &opt, ControlPickCriteria &criteria) {
    criteria.setIndexSorter(noSort);
    // noSort laisse l'ordre tel quel
    // On doit donc, si jamais l'ordre initial a été modifié, le restaurer
    criteria.resetPreferedControlIndexes();
    
    return opt;
}

OptionalCu ClosestPickStrategy::pickControl(const OptionalCu &previousCu, ControlPickCriteria &criteria) {
    using ull = unsigned long long int;

    // Aucun plus proche voisin possible si les strategies précédentes n'ont pas choisies de voisin
    if (previousCu.isRight()) {return previousCu;}

    double rho = criteria.getTimeStep();
    
    pickedControl pickedControl = {rho, 0};
    pickedControl = previousCu.fromLeft(pickedControl);
    const int closestTo = pickedControl.controlIndex;
    
    const ull nbTotalC = criteria.getNbControlCoords();

    double **controlCoords = criteria.getControlCoords();
    // Conversion du vecteur vers un array
    const double *xCoordsDouble = criteria.getCurrentPos();
    
    std::vector<int> indexes(nbTotalC);
    std::iota(indexes.begin(), indexes.end(), 0);
    
    // On trie les indices par distance dans le tableau d'indices avec le point voulu
    std::sort(indexes.begin(), indexes.end(), [closestTo](int i1, int i2){
        return abs(i1 - closestTo) < abs(i2 - closestTo);
    });

    ull cu = 0;
    while (cu < nbTotalC && !criteria.isViableControl(pickedControl)) {
        ++cu;
        pickedControl.controlIndex = cu;
    }
    if (cu < nbTotalC) {
        spdlog::info("Closest strategy control found");
        return OptionalCu(pickedControl);
    }
    else {
        spdlog::error("No control found through closest control strategy");
        return OptionalCu(UNSATISFIED_STRATEGY);
    }
}

OptionalCu PreferedPickStrategy::pickControl(const OptionalCu &previousCu, ControlPickCriteria &criteria) {
    using ull = unsigned long long int;

    // Aucun plus proche voisin possible si les strategies précédentes n'ont pas choisies de voisin
    if (previousCu.isRight()) {return previousCu;}

    double rho = criteria.getTimeStep();
    pickedControl pickedControl = {rho, 0};
    pickedControl = previousCu.fromLeft(pickedControl);
    int pickedCu = pickedControl.controlIndex;
    
    const ull nbTotalC = criteria.getNbControlCoords();
    double **controlCoords = criteria.getControlCoords();
    
    // Conversion du vecteur vers un array
    const double *xCoordsDouble = criteria.getCurrentPos();
    const int *preferedControlIndexes = criteria.sortPreferedControlIndexes();
    const ull indexInPref = std::find(preferedControlIndexes, preferedControlIndexes+nbTotalC, pickedCu) - preferedControlIndexes;

    int prefCu = indexInPref;
    pickedControl.controlIndex = preferedControlIndexes[prefCu];
    
    bool viableControlFound = criteria.isViableControl(pickedControl);       
    while (prefCu >= 0 && !viableControlFound) {
        pickedControl.controlIndex = preferedControlIndexes[prefCu];
        viableControlFound = criteria.isViableControl(pickedControl);
        --prefCu;
    }
    prefCu = (viableControlFound) ? prefCu : indexInPref;
    while (!viableControlFound && (ull) prefCu < nbTotalC) {        
        pickedControl.controlIndex = preferedControlIndexes[prefCu];
        viableControlFound = criteria.isViableControl(pickedControl);
        ++prefCu;
    }
    
    if ((ull) prefCu < nbTotalC) {
        spdlog::info("Closest strategy control found");
        return OptionalCu(pickedControl);
    }
    else {
        spdlog::error("No control found through closest control strategy");
        return OptionalCu(UNSATISFIED_STRATEGY);
    }
}

OptionalCu BubblePickStrategy::findBorderControl(
    ControlPickCriteria &criteria,
    unsigned long long int currentPos,
    unsigned long long int *currentPosIntCoords, double *doubleCoordsOnDiscreteTraj,
    double *imageVect) {

    int nbIter = 0;
    bool onBorder = false;

    const SysDyn *sysDyn = criteria.getSysDyn();
    ControlPicker *picker = criteria.getPicker();
    Trajectory *traj = criteria.getTrajectory();
    TrajectoryPoints *trajDiscrete = criteria.getDiscreteTrajectory();
    double startingTime = criteria.getTime();
    double rho = criteria.getTimeStep();
    
    const int strategyIndex = bubble.getStrategyIndex();
    const Grid *grid = sysDyn->getGrid();
    const int dim = sysDyn->getDim();
    const int dimC = sysDyn->getDimC();
    
    TrajectorySimulation trajCopy = traj->asSimulation();
    TrajectoryPointsSimulation trajDiscreteCopy = trajDiscrete->asSimulation();

    double **controlCoords = sysDyn->getControlCoords();
    
    const double *xCoordsDouble = &(trajCopy.getLastPoint()[0]);
    const double *startingDoubleCoordsOnDiscreteTraj = &(trajDiscreteCopy.getLastPoint()[0]);

    ControlPickCriteria::StrategyIndexBitFlag flag = 0;

    const auto checkViability = [&](pickedControl cu) {
        double pickedRho = cu.timeStep;
        int pickedCu = cu.controlIndex;
                
        if (!criteria.isViableControl(cu, imageVect)) {
            return OptionalCu(UNSATISFIED_STRATEGY);
        }
        else {
            return OptionalCu(cu);
        }
    };

    unsigned long long int imagePos = grid->nbTotalPoints;

    const auto updateTrajectoriesAndCoords = [&](pickedControl cu){
        double pickedRho = cu.timeStep;
        int pickedCu = cu.controlIndex;

        startingTime += pickedRho;
        trajCopy.addPoints(pickedCu, realTimeStepsPerDiscreteStep, startingTime, imageVect, flag);        
        imagePos = grid->getNearestPointInSet(imageVect);

        grid->numToIntAndDoubleCoords(imagePos, currentPosIntCoords, doubleCoordsOnDiscreteTraj);
        trajDiscreteCopy.addPoint(doubleCoordsOnDiscreteTraj, startingTime);
        return OptionalCu(cu);
    };
    
    OptionalCu opt(UNSATISFIED_STRATEGY);

    double time = startingTime;
    
    do {
        flag = 0;
        
        xCoordsDouble = &(trajCopy.getLastPoint()[0]);
        rho = sysDyn->calculRho_local(doubleCoordsOnDiscreteTraj);
        opt = picker
            ->pickControlFromSubPickerList(strategyIndex+1,
                                          trajCopy, trajDiscreteCopy,
                                          rho, flag)
            .applyLeft(checkViability)
            .applyLeft(updateTrajectoriesAndCoords);

            
        while (opt.isRight()) {
            flag = 0;
            xCoordsDouble = &(trajCopy.getLastPoint()[0]);
            rho /= 2;
            opt = picker
                ->pickControlFromSubPickerList(strategyIndex+1,
                                              trajCopy, trajDiscreteCopy,
                                              rho, flag)
                .applyLeft(checkViability)
                .applyLeft(updateTrajectoriesAndCoords);
        }

        onBorder = (bubble.isOnBorder(currentPosIntCoords) || (imagePos == currentPos));

        currentPos = imagePos;
        nbIter++;
        
    } while (bubble.isInBubble(currentPosIntCoords)
             && !onBorder && nbIter < NB_MAX_TRAJ_SIMULATIONS);

    spdlog::info("End of bubble simulation");
    
    if (nbIter >= NB_MAX_TRAJ_SIMULATIONS) {
        spdlog::warn("Simulation never left bubble after {} iterations, delegating to next strategy", nbIter);
        opt = OptionalCu(UNSATISFIED_STRATEGY);
    }
    else if (!onBorder) {
        spdlog::warn("Simulation left bubble, delegating to next strategy");        
        opt = OptionalCu(UNSATISFIED_STRATEGY);
    }
    else if (opt.isRight()) {
        spdlog::warn("No viable control found before reaching border, reusing previous");
        pickedControl c {rho, trajCopy.getLastControlIndex(0)};
        opt = OptionalCu(c);
    }
    else {
        spdlog::info("Border found");
        rho = sysDyn->calculRho_local(doubleCoordsOnDiscreteTraj);
        ControlPickCriteria::StrategyIndexBitFlag flag = 0;
        opt = picker->pickControlFromSubPickerList(strategyIndex+1,
                                                  trajCopy, trajDiscreteCopy,
                                                  rho, flag);
        if (opt.isRight()) {
            spdlog::warn("No viable border control found, reusing previous");
            pickedControl c {rho, trajCopy.getLastControlIndex(0)};
            opt = OptionalCu(c);
        }
        else {
            // Il est intéressant pour la stratégie de bulle de dire quelles
            // stratégies ont contribué au choix de contrôle au bord
            criteria.addContributionFromStrategies(flag);
        }
    }
    
    return opt;
}

OptionalCu BubblePickStrategy::pickControl(const OptionalCu &opt, ControlPickCriteria &criteria) {
    
    using ull = unsigned long long int;
    if (opt.isLeft()) {return opt;}
    const Grid *grid = criteria.getGrid();
    const int dim = criteria.getDim();

    unsigned long long int currentPos = criteria.getCurrentPosGridIndex();
    
    ull *currentPosIntCoords = new ull[dim];
    double *doubleCoordsOnDiscreteTraj = new double[dim];
    double *imageVect = new double[dim];
    grid->numToIntAndDoubleCoords(currentPos, currentPosIntCoords, doubleCoordsOnDiscreteTraj);
    
    bubble.setCenter(currentPosIntCoords);
    OptionalCu borderCu(UNSATISFIED_STRATEGY);
    
    if (!bubble.isTouchingBoundaryInwards()) {
        spdlog::info("Bubble not touching");
    }
    else {
        spdlog::info("Bubble touching - Start of bubble simulation");
        borderCu = findBorderControl(
            criteria, currentPos, currentPosIntCoords, doubleCoordsOnDiscreteTraj,
            imageVect);
    }

    delete[] currentPosIntCoords;
    delete[] doubleCoordsOnDiscreteTraj;
    delete[] imageVect;

    return borderCu;
}

double SmoothPickStrategy::angleBetween(const double *originA, const double *destinationA, const double *originB, const double *destinationB, int dim) const {

    double scalarProductAB = 0;
    double normA = 0;
    double normB = 0;
    
    for (int i = 0; i < dim; ++i) {
        const double vA = destinationA[i] - originA[i];
        const double vB = destinationB[i] - originB[i];
        scalarProductAB += vA*vB;

        normA += vA*vA;
        normB += vB*vB;
    }

    // Deux racines carrées car les normes peuvent être très proches de 0
    const double operand = scalarProductAB / (sqrt(normA) * sqrt(normB));


    if (operand < -1 || 1 < operand) {
        // Si operand n'est pas dans les bornes, et puisque aucune valeur
        // interdite n'existe pour scalarProduct, l'erreur doit provenir de normA ou normB
        //
        // La seule possibilité d'échec est que la norme de A ou de B soit nulle
        // Dans quel cas, on peut supposer que l'angle est respecté, puisque si norme de A
        // ou norme de B valent nul, alors tout angle est autorisé
        // donc on renvoie un angle nul, angle qui passera toujours le test "<= maxAngleRadians"
        spdlog::warn("Angle between vector was invalid, ||A|| = {}, ||B|| = {}, A.B = {}. Returning angle 0", scalarProductAB, normA, normB);
        return 0.0;
    }
    else {
        return acos(operand);
    }
}

double SmoothPickStrategy::squaredDistanceBetween(const double *p1, const double *p2, int dim) const {
    double sum = 0;
    for (int i = 0; i < dim; ++i) {
        double diff = p1[i] - p2[i];
        sum += diff*diff;

        // On a que
        // (x - y)^2 = x^2 + y^2 -2xy,
        // cette formule est plus stable numériquement pour des points p1 et p2 proches
        // Mais elle est (probablement) moins efficace à calculer
        //sum += p1[i]*p1[i] + p2[i]*p2[i] - 2*p1[i]*p2[i];
    }
    return sum;
}

OptionalCu SmoothPickStrategy::smoothControlClosestTo(pickedControl &p,
                                                      const double *previousToLastPosition, const double *lastPosition,
                                                      ControlPickCriteria &criteria, double *imageVect) {

    int chosenCu = p.controlIndex;
    double rho = p.timeStep;

    SysDyn *sysDyn = criteria.getSysDyn();
    
    double **controlCoords = sysDyn->getControlCoords();
    unsigned long long int nbTotalC = sysDyn->getTotalNbPointsC();
    const int dim = sysDyn->getDim();
    
    int bestCu = 0;
    double smallestSquaredDistance = PLUS_INF;
    // On échange les deux pour ne pas avoir à vérifier une condition supplémentaire dans le for
    std::swap(controlCoords[0], controlCoords[chosenCu]);    

    double *chosenNextPosition = new double[dim];
    std::copy(imageVect, imageVect+dim, chosenNextPosition);

    
    for (unsigned long long int cu = 1; cu < nbTotalC; ++cu) {
        p.controlIndex = cu;
        if (criteria.isViableControl(p, imageVect)) {
                
            const double squaredDistance = squaredDistanceBetween(imageVect, chosenNextPosition, dim);
                
            if (squaredDistance < smallestSquaredDistance
                && angleBetween(previousToLastPosition, lastPosition,
                                lastPosition, imageVect, dim) <= maxAngleRadians) {
                        
                bestCu = cu;
                smallestSquaredDistance = squaredDistance;
            }
        }
    }
    std::swap(controlCoords[0], controlCoords[chosenCu]);

    if (bestCu == 0) {
        spdlog::warn("Unsatisfied smooth pick strategy");
        return OptionalCu(UNSATISFIED_STRATEGY);
    }
    // Comme on a échangé ces deux indices, il faut faire attention à renvoyer
    // l'indice correct
    else if (bestCu == chosenCu) {
        bestCu = 0;
    }

    spdlog::info("Smooth control found");

    delete[] chosenNextPosition;

    p.controlIndex = bestCu;
        
    return OptionalCu(p);
}

OptionalCu SmoothPickStrategy::pickControl(const OptionalCu &opt, ControlPickCriteria &criteria) {

    const std::vector<Trajectory::point> &points = criteria.getPoints();
    
    // Une trajectoire contenant uniquement un point de départ sera toujours lisse
    if (opt.isRight() || points.size() <= 1) { return opt; }
    
    const int dim = criteria.getDim();
    double **controlCoords = criteria.getControlCoords();

    double imageVect[dim];
    
    // La syntaxe &(v[0]) permet de convertir un vecteur (stockage contigu)
    // en pointeur C
    const double *lastPosition = criteria.getCurrentPos();    
    const double *previousToLastPosition = &(points.end()[-2][0]);
    double rho = criteria.getTimeStep();
    
    pickedControl p = {rho, 0};
    p = opt.fromLeft(p);

    const int chosenCu = p.controlIndex;
    rho = p.timeStep;

    // Si le contrôle est viable et que l'angle respecte déjà la contrainte,
    // on accepte le contrôle proposé
    if (criteria.isViableControl(p, imageVect)
        && angleBetween(previousToLastPosition, lastPosition,
                        lastPosition, imageVect, dim) <= maxAngleRadians) {
        spdlog::info("Control produces valid smooth trajectory, continuing");
        return opt;
    }
    else {
        return smoothControlClosestTo(p,
                                      previousToLastPosition, lastPosition,
                                      criteria, imageVect);
    }
}

double TemporalControlPickStrategy::squaredDistanceBetween(const double *point1, const double *point2, int dim) const {
    double sum = 0;
    for (int i = 0; i < dim; ++i) {
        double diff = point1[i] - point2[i];
        sum += diff*diff;

        // On a que
        // (x - y)^2 = x^2 + y^2 -2xy,
        // cette formule est plus stable numériquement pour des points p1 et p2 proches
        // Mais elle est (probablement) moins efficace à calculer
        //sum += point1[i]*point1[i] + point2[i]*point2[i] - 2*point1[i]*point2[i];
    }
    return sum;
}

int TemporalControlPickStrategy::getClosestControlTo(double *userControl, double **controlCoords, unsigned long long int nbTotalC, int dimC) const {
    int closest = 0;
    double smallestDistance = squaredDistanceBetween(userControl, controlCoords[0], dimC);

    for (unsigned long long int i = 1; i < nbTotalC; ++i) {
        double tmpDistance = squaredDistanceBetween(userControl, controlCoords[i], dimC);
        if (tmpDistance < smallestDistance) {
            smallestDistance = tmpDistance;
            closest = i;
        }
    }
    return closest;
}

OptionalCu TemporalControlPickStrategy::pickControl(const OptionalCu &opt,
                                                ControlPickCriteria &criteria) {

    if (opt.isLeft()) {return opt;}
    
    double **controlCoords = criteria.getControlCoords();
    const unsigned long long int nbTotalC = criteria.getNbControlCoords();
    const double time = criteria.getTime();
    const int trajIndex = criteria.getTrajectoryIndex();
    const int stratIndex = criteria.getStrategyIndex();
    const int dimC = criteria.getDimC();
    
    double *userControl = new double[dimC];
    
    double rho = criteria.getGridTimeStep();
    
    getTemporalControl(time, trajIndex, stratIndex, userControl);
    /* Il se pourrait que le contrôle choisi en t ne soit pas
       le même que celui en t+rho

       On pourrait chercher une valeur approché du pas de temps pendant
       lequel le contrôle choisi en t reste le même et renvoyer ce pas
       de temps avec le contrôle. Mais on va supposer rho "assez petit"
       pour que ce calcul ne soit pas nécessaire à la satisfaction de
       l'utilisateur */
    int controlIndex = getClosestControlTo(userControl, controlCoords, nbTotalC, dimC);
    pickedControl pickedControl {rho, controlIndex};

    criteria.findViableTimeStep(pickedControl);
    
    if (criteria.isSameGridPosition(pickedControl)) {
        spdlog::warn("Temporal control not viable, delegating to next strategy");
        return OptionalCu(UNSATISFIED_STRATEGY);
    }
    else {
        spdlog::info("Viable temporal control found");
        return OptionalCu(pickedControl);        
    }
}

OptionalCu BubbleBorderPickStrategy::pickControl(const OptionalCu &opt,
                                                 ControlPickCriteria &criteria) {
    using ull = unsigned long long int;
    if (opt.isLeft()) {return opt;}
    const Grid *grid = criteria.getGrid();
    const int dim = criteria.getDim();

    OptionalCu res(UNSATISFIED_STRATEGY);
    
    unsigned long long int currentPos = criteria.getCurrentPosGridIndex();
    
    ull *currentPosIntCoords = new ull[dim];
    double *doubleCoordsOnDiscreteTraj = new double[dim];
    double *imageVect = new double[dim];
    grid->numToIntAndDoubleCoords(currentPos, currentPosIntCoords, doubleCoordsOnDiscreteTraj);
    
    bubble.setCenter(currentPosIntCoords);

    if (!bubble.isTouchingBoundaryInwards()) {
        spdlog::info("Bubble not touching");
    }
    else {
        spdlog::info("Bubble touching");
        bubble.getClosestPointOnBorder(currentPosIntCoords);
        grid->intCoordsToDoubleCoords(currentPosIntCoords, doubleCoordsOnDiscreteTraj);
        
        const SysDyn *sysDyn = criteria.getSysDyn();
        ControlPicker *picker = criteria.getPicker();
        const int strategyIndex = bubble.getStrategyIndex();

        // Création de trajectoires contenant comme seul point le point de bord de bulle
        // Ici, on n'utilise pas une simulation car certaines stratégies décident du
        // contrôle en fonction de l'historique des positions (exemple: SMOOTH) et ce
        // déplacement étant impossible, cela risque de fausser leurs valeurs de retour
        TrajectoryStorage traj = TrajectoryStorage::createNonFlagSavingStorage(doubleCoordsOnDiscreteTraj, PLUS_INF, sysDyn);
        TrajectoryPointsStorage trajDiscrete(doubleCoordsOnDiscreteTraj, PLUS_INF, sysDyn);
        double rho = sysDyn->calculRho_local(doubleCoordsOnDiscreteTraj);
        ControlPickCriteria::StrategyIndexBitFlag flag = 0;
        double time = criteria.getTime();
        
        // On demande au pickers suivants de choisir un contrôle "comme si" il étaient au bord
        res = picker->pickControlFromSubPickerList(strategyIndex+1,
                                                   traj, trajDiscrete,
                                                   rho, flag);

        pickedControl picked;
        picked.timeStep = rho;
        picked.controlIndex = 0;
        
        if (res.isRight()) {
            spdlog::warn("No viable border control found.");
            res = OptionalCu(UNSATISFIED_STRATEGY);
        }
        else if (!criteria.isViableControl((picked = res.fromLeft(picked)))) {
            spdlog::warn("Border control not viable in current position");
            res = OptionalCu(UNSATISFIED_STRATEGY);
        }
        else {
            spdlog::info("Border control found");
            // Il est intéressant pour la stratégie de bulle de dire quelles
            // stratégies ont contribué au choix de contrôle au bord
            criteria.addContributionFromStrategies(flag);
        }
    }

    delete[] currentPosIntCoords;
    delete[] doubleCoordsOnDiscreteTraj;
    delete[] imageVect;

    return res;
}



void UserPickStrategy::setName(const std::string &newName) {
    name = newName;
}

FirstPickStrategyBitSet::FirstPickStrategyBitSet(ViabiBitSetTrajectoryHelper *viabiHelper) :
    viabiHelper(viabiHelper) {}

HeavyPickStrategy::HeavyPickStrategy(int initCu) :
    initCu(initCu) {}

BubblePickStrategy::BubblePickStrategy(int strategyIndex, const TrajectoryParametersManager *params, const Grid *grid) :
    bubble(params, grid, strategyIndex),
    realTimeStepsPerDiscreteStep(params->getTrajectoryParameters()->REAL_TIME_STEPS_PER_DISCRETE_STEP) {}

BubbleBorderPickStrategy::BubbleBorderPickStrategy(int strategyIndex, const TrajectoryParametersManager *params, const Grid *grid) :
    bubble(params, grid, strategyIndex) {}

SmoothPickStrategy::SmoothPickStrategy(double maxAngleRadians) :
    maxAngleRadians(maxAngleRadians) {}

TemporalControlPickStrategy::TemporalControlPickStrategy(const TrajectoryParametersManager *params) :
    getTemporalControl(params->getTrajectoryParameters()->TEMPORAL_CONTROL)
{}

const std::string HeavyPickStrategy::name = toString(HEAVY);
const std::string &HeavyPickStrategy::getName() const {
    return name;
}

const std::string FirstPickStrategyBitSet::name = toString(FIRST);
const std::string &FirstPickStrategyBitSet::getName() const {
    return name;
}

const std::string FirstOnlyPickStrategy::name = toString(FIRST_ONLY);
const std::string &FirstOnlyPickStrategy::getName() const {
    return name;
}

const std::string ShuffleStrategy::name = toString(SHUFFLE);
const std::string &ShuffleStrategy::getName() const {
    return name;
}

const std::string SortStrategy::name = toString(SORT);
const std::string &SortStrategy::getName() const {
    return name;
}

const std::string ResetOrderStrategy::name = toString(RESET_ORDER);
const std::string &ResetOrderStrategy::getName() const {
    return name;
}

const std::string ClosestPickStrategy::name = toString(CLOSEST);
const std::string &ClosestPickStrategy::getName() const {
    return name;
}

const std::string PreferedPickStrategy::name = toString(PREFERED);
const std::string &PreferedPickStrategy::getName() const {
    return name;
}

const std::string BubblePickStrategy::name = toString(BUBBLE);
const std::string &BubblePickStrategy::getName() const {
   return name;
}

const std::string SmoothPickStrategy::name = toString(SMOOTH);
const std::string &SmoothPickStrategy::getName() const {
    return name;
}

const std::string TemporalControlPickStrategy::name = toString(TEMPORAL_CONTROL);
const std::string &TemporalControlPickStrategy::getName() const {
    return name;
}

const std::string BubbleBorderPickStrategy::name = toString(BUBBLE_BORDER);
const std::string &BubbleBorderPickStrategy::getName() const {
    return name;
}

const std::string &UserPickStrategy::getName() const {
    return name;
}
