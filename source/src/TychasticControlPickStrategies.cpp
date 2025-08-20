#include "../include/TychasticControlPickStrategies.h"

#include "../include/TrajectoryPointsSimulation.h"
#include "../include/TychasticTrajectorySimulation.h"

TychasticFirstGuaranteedPickStrategy::TychasticFirstGuaranteedPickStrategy(ViabiTrajectoryHelper *viabiHelper) :
    trajectoryHelper(viabiHelper) {}

OptionalCu TychasticFirstGuaranteedPickStrategy::findViableDiscreteControl(const double *xCoordsDouble, TychasticControlPickCriteria &criteria) {
    
    SysDyn *sysDyn = criteria.getSysDyn();
    double **controlCoords = sysDyn->getControlCoords();
    unsigned long long int nbTotalC = sysDyn->getTotalNbPointsC();

    TychasticControlPicker *picker = criteria.getPicker();
    
    int *preferedControlIndexes = picker->sortPreferedControlIndexes(xCoordsDouble,
                                                                     criteria.getTime(),
                                                                     criteria.getStrategyIndex());
    double rho = criteria.getTimeStep();

    unsigned long long int prefCu = 0;
    bool viableFound = false;
    do {
        viableFound = trajectoryHelper->isViableGuaranteedControl(xCoordsDouble, controlCoords[preferedControlIndexes[prefCu]], rho);
    } while (!viableFound && (++prefCu) < nbTotalC);
    
    if (!viableFound) {
        return OptionalCu(UNSATISFIED_STRATEGY);
    }
    else {
        pickedControl p;
        p.timeStep = rho;
        p.controlIndex = prefCu;
        return OptionalCu(p);
    }
}

OptionalCu TychasticFirstGuaranteedPickStrategy::pickControl(const OptionalCu &opt,
                                                   TychasticControlPickCriteria &criteria) {

    if (opt.isLeft()) { return opt; }
    
    pickedControl pick = {criteria.getTimeStep(), 0};
    const double *xCoordsDouble = criteria.getCurrentDiscretePos();
    
    OptionalCu res = findViableDiscreteControl(xCoordsDouble, criteria);

    if (res.isRight()) {
        spdlog::error("No grid control found with first pick strategy");
    }
    else if (criteria.isViableGuaranteedControl((pick = res.fromLeft(pick)), 0)) {
        spdlog::info("First strategy control found");
    }
    else {
        // Transformation contrôle discret -> continu pas évidente...
        // Pour l'instant on joue seulement sur le pas de temps
        criteria.findViableGuaranteedTimeStep(pick, 0);
        if (criteria.isGuaranteedSameGridPosition(pick, 0)) {
            spdlog::warn("No viable real control found to match with discrete position with first viable strategy. Returning grid control n°{}", pick.controlIndex);
            res = OptionalCu(UNSATISFIED_STRATEGY);
        }
        else {
            res = OptionalCu(pick);
        }
    }
    return res;
}

OptionalCu TychasticHeavyPickStrategy::pickControl(const OptionalCu &opt, TychasticControlPickCriteria &criteria) {
        
    if (opt.isLeft()) {return opt;}
    
    using ull = unsigned long long int;

    double rho = criteria.getTimeStep();
    const Grid *grid = criteria.getGrid();
    const int dim = criteria.getDim();    
    const double *xCoordsDouble = criteria.getCurrentPos();
    
    TychasticTrajectory *traj = criteria.getTrajectory();
    TrajectoryPoints *trajDiscrete = criteria.getDiscreteTrajectory();
    
    double **controlCoords = criteria.getControlCoords();  
    int lastControlIndex = traj->getLastControlIndex(initCu);
    const double *lastControl = controlCoords[lastControlIndex];

    double imageVect[dim];
    pickedControl pickedControl{rho, lastControlIndex};

    int tyIndex = criteria.getCurrentTycheIndex();
    std::invoke(criteria.findViableTimeStep, criteria, pickedControl, tyIndex);

    // Si, pour un pas de temps ayant décrû jusqu'à ce qu'il n'y ait pas de différence
    // de position discrète aucun pas de temps ne permet d'atteindre une position viable,
    // on peut supposer que le contrôle lourd ne marche plus
    if (std::invoke(criteria.isSameGridPosition, criteria, pickedControl, tyIndex)) {
        spdlog::warn("Unsatisfied heavy pick strategy with control n°{}", lastControlIndex);
        return OptionalCu(UNSATISFIED_STRATEGY);
    }
    else {
        spdlog::info("Heavy control found");                
        return OptionalCu(pickedControl);
    }    
}

OptionalCu TychasticFirstNonGuaranteedPickStrategyBitSet::pickControl(const OptionalCu &opt, TychasticControlPickCriteria &criteria) {

    if (opt.isLeft()) {return opt;}
    
    const Grid *grid = criteria.getGrid();
    // Conversion du vecteur vers un array
    const double *xCoordsDouble = criteria.getCurrentPos();
    unsigned long long int currentPos = criteria.getCurrentPosGridIndex();
    int tycheIndex = criteria.getCurrentTycheIndex();
    const double time = criteria.getTime();
    double rho = criteria.getTimeStep();
    const int dim = criteria.getDim();

    double **controlCoords = criteria.getControlCoords();

    double imageVect[dim];
    
    int usedCu;

    const unsigned long long viabSuccessor = viabiHelper->findViableDiscreteSuccessor_tych(currentPos, time, rho, tycheIndex, usedCu);

    double discreteRho = rho;
    
    bool succes = false;
    if (viabSuccessor > grid->nbTotalPoints) {
        spdlog::error("Optimal successor not found with first viable strategy");
        return OptionalCu(UNSATISFIED_STRATEGY);
    }
    else if (const int realCu = viabiHelper->findViabControl_bis_tych(xCoordsDouble, viabSuccessor, rho, tycheIndex, 25, 0.5, imageVect, succes); succes) {
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

OptionalCu TychasticFirstOnlyPickStrategy::pickControl(
    const OptionalCu &opt, TychasticControlPickCriteria &criteria) {

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
    
    int tyIndex = criteria.getCurrentTycheIndex();
    bool viable = std::invoke(criteria.isViableControl, criteria, c, tyIndex);
    
    if (!viable) {
        spdlog::warn("Optimal successor not found with first viable strategy");
        return OptionalCu(UNSATISFIED_STRATEGY);
    }
    else {
        spdlog::info("First only strategy control found");
        return OptionalCu(c);
    }    
}

OptionalCu TychasticShuffleStrategy::pickControl(const OptionalCu &opt, TychasticControlPickCriteria &criteria) {
    criteria.setIndexSorter(shuffleControlIndexes);
    return opt;
}

OptionalCu TychasticSortStrategy::pickControl(const OptionalCu &opt, TychasticControlPickCriteria &criteria) {
    criteria.setIndexSorter(sortControlIndexesByWeight);
    return opt;
}

OptionalCu TychasticResetOrderStrategy::pickControl(const OptionalCu &opt, TychasticControlPickCriteria &criteria) {
    criteria.setIndexSorter(noSort);
    // noSort laisse l'ordre tel quel
    // On doit donc, si jamais l'ordre initial a été modifié, le restaurer
    criteria.resetPreferedControlIndexes();
    
    return opt;
}

OptionalCu TychasticClosestPickStrategy::pickControl(const OptionalCu &previousCu, TychasticControlPickCriteria &criteria) {
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

    int tyIndex = criteria.getCurrentTycheIndex();
    
    std::vector<int> indexes(nbTotalC);
    std::iota(indexes.begin(), indexes.end(), 0);
    
    // On trie les indices par distance dans le tableau d'indices avec le point voulu
    std::sort(indexes.begin(), indexes.end(), [closestTo](int i1, int i2){
        return abs(i1 - closestTo) < abs(i2 - closestTo);
    });

    ull cu = 0;
    while (cu < nbTotalC && !std::invoke(criteria.isViableControl, criteria, pickedControl, tyIndex)) {
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

OptionalCu TychasticPreferedPickStrategy::pickControl(const OptionalCu &previousCu, TychasticControlPickCriteria &criteria) {
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

    int tyIndex = criteria.getCurrentTycheIndex();
    
    bool viableControlFound = std::invoke(criteria.isViableControl, criteria, pickedControl, tyIndex);
    while (prefCu >= 0 && !viableControlFound) {
        pickedControl.controlIndex = preferedControlIndexes[prefCu];
        viableControlFound = std::invoke(criteria.isViableControl, criteria, pickedControl, tyIndex);
        --prefCu;
    }
    prefCu = (viableControlFound) ? prefCu : indexInPref;
    while (!viableControlFound && (ull) prefCu < nbTotalC) {        
        pickedControl.controlIndex = preferedControlIndexes[prefCu];
        viableControlFound = std::invoke(criteria.isViableControl, criteria, pickedControl, tyIndex);
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

double TychasticSmoothPickStrategy::angleBetween(const double *originA, const double *destinationA, const double *originB, const double *destinationB, int dim) const {

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

double TychasticSmoothPickStrategy::squaredDistanceBetween(const double *p1, const double *p2, int dim) const {
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

OptionalCu TychasticSmoothPickStrategy::smoothControlClosestTo(pickedControl &p,
                                                      const double *previousToLastPosition, const double *lastPosition,
                                                      TychasticControlPickCriteria &criteria, double *imageVect) {

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

    int tyIndex = criteria.getCurrentTycheIndex();
    
    for (unsigned long long int cu = 1; cu < nbTotalC; ++cu) {
        p.controlIndex = cu;
        if (criteria.isViableControlForTyche(p, tyIndex)) {

            criteria.getTychasticImage(p, tyIndex, imageVect);            
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

OptionalCu TychasticSmoothPickStrategy::pickControl(const OptionalCu &opt, TychasticControlPickCriteria &criteria) {

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
    int tyIndex = criteria.getCurrentTycheIndex();
    rho = p.timeStep;

    // Si le contrôle est viable et que l'angle respecte déjà la contrainte,
    // on accepte le contrôle proposé
    if (criteria.isViableControlForTyche(p, tyIndex)) {
        criteria.getTychasticImage(p, tyIndex, imageVect);
        if (angleBetween(previousToLastPosition, lastPosition,
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
    else {
        return smoothControlClosestTo(p,
                                      previousToLastPosition, lastPosition,
                                      criteria, imageVect);
    }
}

double TychasticTemporalControlPickStrategy::squaredDistanceBetween(const double *point1, const double *point2, int dim) const {
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

int TychasticTemporalControlPickStrategy::getClosestControlTo(double *userControl, double **controlCoords, unsigned long long int nbTotalC, int dimC) const {
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

OptionalCu TychasticTemporalControlPickStrategy::pickControl(const OptionalCu &opt,
                                                             TychasticControlPickCriteria &criteria) {

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
    
    int tyIndex = criteria.getCurrentTycheIndex();    

    std::invoke(criteria.findViableTimeStep, criteria, pickedControl, tyIndex);
    
    if (std::invoke(criteria.isSameGridPosition, criteria, pickedControl, tyIndex)) {
        spdlog::warn("Temporal control not viable, delegating to next strategy");
        return OptionalCu(UNSATISFIED_STRATEGY);
    }
    else {
        spdlog::info("Viable temporal control found");
        return OptionalCu(pickedControl);        
    }
}

OptionalCu TychasticBubbleBorderPickStrategy::pickControl(const OptionalCu &opt,
                                                          TychasticControlPickCriteria &criteria) {
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
        TychasticControlPicker *picker = criteria.getPicker();
        const int strategyIndex = bubble.getStrategyIndex();
        
        // Création de trajectoires contenant comme seul point le point de bord de bulle
        // Ici, on n'utilise pas une simulation car certaines stratégies décident du
        // contrôle en fonction de l'historique des positions (exemple: SMOOTH) et ce
        // déplacement étant impossible, cela risque de fausser leurs valeurs de retour
        TychasticTrajectoryStorage traj = TychasticTrajectoryStorage::createNonFlagSavingStorage(doubleCoordsOnDiscreteTraj, PLUS_INF, sysDyn);
        TrajectoryPointsStorage trajDiscrete(doubleCoordsOnDiscreteTraj, PLUS_INF, sysDyn);
        double rho = sysDyn->calculRho_local(doubleCoordsOnDiscreteTraj);
        TychasticControlPickCriteria::StrategyIndexBitFlag flag = 0;
        double time = criteria.getTime();
        
        // On demande au pickers suivants de choisir un contrôle "comme si" il étaient au bord
        res = picker->pickControlFromSubPickerList(strategyIndex+1,
                                                   traj, trajDiscrete,
                                                   rho, flag);

        int tyIndex = criteria.getCurrentTycheIndex();
        pickedControl picked;
        picked.timeStep = rho;
        picked.controlIndex = 0;
        
        if (res.isRight()) {
            spdlog::warn("No viable border control found");
            res = OptionalCu(UNSATISFIED_STRATEGY);
        }
        else if (!std::invoke(criteria.isViableControl, criteria, (picked = res.fromLeft(picked)), tyIndex)) {
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

TychasticFirstNonGuaranteedPickStrategyBitSet::TychasticFirstNonGuaranteedPickStrategyBitSet(ViabiBitSetTrajectoryHelper *viabiHelper) :
    viabiHelper(viabiHelper) {}

TychasticHeavyPickStrategy::TychasticHeavyPickStrategy(int initCu) :
    initCu(initCu) {}


TychasticBubbleBorderPickStrategy::TychasticBubbleBorderPickStrategy(int strategyIndex, const TrajectoryParametersManager *params, const Grid *grid) :
    bubble(params, grid, strategyIndex) {}

TychasticSmoothPickStrategy::TychasticSmoothPickStrategy(double maxAngleRadians) :
    maxAngleRadians(maxAngleRadians) {}

TychasticTemporalControlPickStrategy::TychasticTemporalControlPickStrategy(const TrajectoryParametersManager *params) :
    getTemporalControl(params->getTrajectoryParameters()->TEMPORAL_CONTROL)
{}

void TychasticUserPickStrategy::setName(const std::string &newName) {
    name = newName;
}

const std::string &TychasticUserPickStrategy::getName() const {
    return name;
}

const std::string TychasticFirstGuaranteedPickStrategy::name = toString(FIRST);
const std::string &TychasticFirstGuaranteedPickStrategy::getName() const {
    return name;
}

const std::string TychasticHeavyPickStrategy::name = toString(HEAVY);
const std::string &TychasticHeavyPickStrategy::getName() const {
    return name;
}

const std::string TychasticFirstNonGuaranteedPickStrategyBitSet::name = toString(FIRST);
const std::string &TychasticFirstNonGuaranteedPickStrategyBitSet::getName() const {
    return name;
}

const std::string TychasticFirstOnlyPickStrategy::name = toString(FIRST_ONLY);
const std::string &TychasticFirstOnlyPickStrategy::getName() const {
    return name;
}

const std::string TychasticShuffleStrategy::name = toString(SHUFFLE);
const std::string &TychasticShuffleStrategy::getName() const {
    return name;
}

const std::string TychasticSortStrategy::name = toString(SORT);
const std::string &TychasticSortStrategy::getName() const {
    return name;
}

const std::string TychasticResetOrderStrategy::name = toString(RESET_ORDER);
const std::string &TychasticResetOrderStrategy::getName() const {
    return name;
}

const std::string TychasticClosestPickStrategy::name = toString(CLOSEST);
const std::string &TychasticClosestPickStrategy::getName() const {
    return name;
}

const std::string TychasticPreferedPickStrategy::name = toString(PREFERED);
const std::string &TychasticPreferedPickStrategy::getName() const {
    return name;
}

const std::string TychasticSmoothPickStrategy::name = toString(SMOOTH);
const std::string &TychasticSmoothPickStrategy::getName() const {
    return name;
}

const std::string TychasticTemporalControlPickStrategy::name = toString(TEMPORAL_CONTROL);
const std::string &TychasticTemporalControlPickStrategy::getName() const {
    return name;
}

const std::string TychasticBubbleBorderPickStrategy::name = toString(BUBBLE_BORDER);
const std::string &TychasticBubbleBorderPickStrategy::getName() const {
    return name;
}

