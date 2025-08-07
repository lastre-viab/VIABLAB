#ifndef TYCHASTICCONTROLPICKSTRATEGIES_H
#define TYCHASTICCONTROLPICKSTRATEGIES_H

#include "ViabiBitSetTrajectoryHelper.h"
#include "TychasticControlPickStrategy.h"
#include "Bubble.h"

// Ces stratégies sont dans un fichier différent de ControlPickStrategy pour permettre à
// l'utilisateur de n'avoir à inclure que le fichier TychasticControlPickStrategy.h

class TychasticFirstGuaranteedPickStrategy final : public TychasticControlPickStrategy {
public:
    TychasticFirstGuaranteedPickStrategy() = default;
    OptionalCu pickControl(const OptionalCu &previousStrategyResult,
                           TychasticControlPickCriteria &criteria) override final;
    const std::string &getName() const override final;
    virtual ~TychasticFirstGuaranteedPickStrategy() = default;
private:

    OptionalCu findViableDiscreteControl(const double *xCoordsDouble, TychasticControlPickCriteria &criteria);
    
    static const std::string name;
};

/*!
  Choisit le même que l'itération précédente (le contrôle renvoyé par getLastControlIndex de traj) en modifiant possiblement le pas de temps rho.
  Si celui-ci n'est pas viable, renvoie UNSATISFIED_STRATEGY
*/
class TychasticHeavyPickStrategy final : public TychasticControlPickStrategy {
public:
    TychasticHeavyPickStrategy(int initCu);
    OptionalCu pickControl(const OptionalCu &previousStrategyResult,
                           TychasticControlPickCriteria &criteria) override final;
    const std::string &getName() const override final;
    virtual ~TychasticHeavyPickStrategy() = default;
private:
    static const std::string name;
    int initCu;
};

/*!
  Choisit le premier contrôle qui, pour le pas de temps donné en paramètre,
  rapproche le plus possible le contrôle réel du premier contrôle de grille valide.

  L'ordre déterminant le premier contrôle de grille valide est celui de getPreferedControlIndexes de sysDyn.

  Si le contrôle réel n'est pas viable, renvoie celui de grille.
  Si le contrôle de grille n'est pas viable, renvoie UNSATISFIED_STRATEGY
 */
class TychasticFirstNonGuaranteedPickStrategyBitSet final : public TychasticControlPickStrategy {
public:
    TychasticFirstNonGuaranteedPickStrategyBitSet(ViabiBitSetTrajectoryHelper *);
    OptionalCu pickControl(const OptionalCu &previousStrategyResult,
                           TychasticControlPickCriteria &criteria) override final;
    const std::string &getName() const override final;
    virtual ~TychasticFirstNonGuaranteedPickStrategyBitSet() = default;
private:
    static const std::string name;
    ViabiBitSetTrajectoryHelper *viabiHelper;
};

/*!
  Renvoie le premier contrôle viable réel.

  L'ordre déterminant le premier contrôle de grille valide est celui des preferedControlIndexes de sysDyn.
  
  Si le premier n'est pas viable, renvoie UNSATISFIED_STRATEGY
 */
class TychasticFirstOnlyPickStrategy final : public TychasticControlPickStrategy {
public:
    TychasticFirstOnlyPickStrategy() = default;
    OptionalCu pickControl(const OptionalCu &previousStrategyResult,
                           TychasticControlPickCriteria &criteria) override final;
    const std::string &getName() const override final;
    virtual ~TychasticFirstOnlyPickStrategy() = default;
private:
    static const std::string name;    
};

/*!
  Réordonne les preferedControlIndexes de manière uniformément aléatoire

  Renvoie toujours UNSATISFIED_STRATEGY
 */
class TychasticShuffleStrategy final : public TychasticControlPickStrategy {
public:
    TychasticShuffleStrategy() = default;
    OptionalCu pickControl(const OptionalCu &previousStrategyResult,
                           TychasticControlPickCriteria &criteria) override final;
    const std::string &getName() const override final;
    virtual ~TychasticShuffleStrategy() = default;
private:
    static const std::string name;
};

/*!
  Réordonne les preferedControlIndexes selon l'ordre dicté par
  la fonction utilisateur controlWeight par ordre décroissant des poids.
 */
class TychasticSortStrategy final : public TychasticControlPickStrategy {
public:
    TychasticSortStrategy() = default;
    OptionalCu pickControl(const OptionalCu &previousStrategyResult,
                           TychasticControlPickCriteria &criteria) override final;
    const std::string &getName() const override final;
    virtual ~TychasticSortStrategy() = default;
private:
    static const std::string name;
};

/*!
  Réapplique l'ordre par défaut des contrôles.
 */
class TychasticResetOrderStrategy final : public TychasticControlPickStrategy {
public:
    TychasticResetOrderStrategy() = default;
    OptionalCu pickControl(const OptionalCu &previousStrategyResult,
                           TychasticControlPickCriteria &criteria) override final;
    const std::string &getName() const override final;
    virtual ~TychasticResetOrderStrategy() = default;
private:
    static const std::string name;
};

/*!
  Choisit le contrôle le plus proche de previousStrategyResult viable.
  Ce qui est entendu par "le plus proche" est ici la différence d'indice la plus faible
  dans le tableau preferedControlCoords de sysDyn.
  
  Si previousStrategyResult n'est pas un contrôle, renvoie previousStrategyResult.
  Sinon, renvoie le contrôle le plus proche de previousStrategyResult viable pour
  la trajectoire réelle.
 */
class TychasticClosestPickStrategy final : public TychasticControlPickStrategy {
public:
    TychasticClosestPickStrategy() = default;
    OptionalCu pickControl(const OptionalCu &previousStrategyResult,
                           TychasticControlPickCriteria &criteria) override final;
    const std::string &getName() const override final;
    virtual ~TychasticClosestPickStrategy() = default;
private:
    static const std::string name;
};

/*!
  Choisit le contrôle le plus proche de previousStrategyResult viable
  préféré par apport à previousStrategyResult. Si aucun des préférés ne fonctionne,
  retourne simplement le plus proche viable.
  
  Ce qui est entendu par "le plus proche" est ici la différence d'indice la plus faible
  dans le tableau preferedControlCoords de sysDyn.
  
  Si previousStrategyResult n'est pas un contrôle, renvoie previousStrategyResult.
  Sinon, renvoie le contrôle vérifiant les conditions ci-dessus.
 */
class TychasticPreferedPickStrategy final : public TychasticControlPickStrategy {
public:
    TychasticPreferedPickStrategy() = default;
    OptionalCu pickControl(const OptionalCu &previousStrategyResult,
                           TychasticControlPickCriteria &criteria) override final;
    const std::string &getName() const override final;
    virtual ~TychasticPreferedPickStrategy() = default;
private:
    static const std::string name;
};

/*!
  Choisit le contrôle le plus proche du contrôle passé en paramètre respectant
  la contrainte de l'angle maximale donnée au constructeur de la stratégie.

  La contrainte de l'angle signifie que l'angle maximal entre le segment formé
  par les 2 derniers points de la trajectoire et le segment formé par le dernier
  point et le point d'arrivée du contrôle ne doit pas dépasser la valeur
  maxAngleRadians

  Si le contrôle choisi respecte cette contrainte, il est renvoyé tel quel
  Si previousStrategyResult n'est pas un contrôle, il est renvoyé tel quel
  Sinon, on renvoie le contrôle tel que le point d'arrivé respectant la
  contrainte le plus proche du point d'arrivée du contrôle previousStrategyResult

  Si aucun contrôle ne satisfait ce critère, le contrôle initial est renvoyé tel
  quel et un message de warning est affiché
 */
class TychasticSmoothPickStrategy final : public TychasticControlPickStrategy {
public:
    TychasticSmoothPickStrategy(double maxAngleRadians);

    OptionalCu pickControl(const OptionalCu &previousStrategyResult,
                           TychasticControlPickCriteria &criteria) override final;
    const std::string &getName() const override final;
    virtual ~TychasticSmoothPickStrategy() = default;
private:
    double angleBetween(const double *originA, const double *destinationA, const double *originB, const double *destinationB, int dim) const;
    double squaredDistanceBetween(const double *point1, const double *point2, int dim) const;

    OptionalCu smoothControlClosestTo(pickedControl &p,
                                      const double *previousToLastPosition, const double *lastPosition,
                                      TychasticControlPickCriteria &criteria, double *imageVect);
    
    double maxAngleRadians;
    static const std::string name;
};

/*!
  Choisit l'indice du contrôle le plus proche du contrôle renvoyé par la fonction
  temporalControl de l'utilisateur.
 */
class TychasticTemporalControlPickStrategy final : public TychasticControlPickStrategy {
public:
    TychasticTemporalControlPickStrategy(const TrajectoryParametersManager *);
    OptionalCu pickControl(const OptionalCu &previousStrategyResult,
                           TychasticControlPickCriteria &criteria) override final;
    const std::string &getName() const override final;

    int getClosestControlTo(double *userControl, double **controlCoords, unsigned long long int nbTotalC, int dimC) const;
    double squaredDistanceBetween(const double *point1, const double *point2, int dim) const;
    
    virtual ~TychasticTemporalControlPickStrategy() = default;
private:
    temporalControl_t getTemporalControl;
    static const std::string name;
};

/*!
  Renvoie le contrôle choisi par la suite de la liste de stratégies
  (par rapport à celle-ci) au point le plus proche (au sens de la
  distance euclidienne) du centre de la bulle situé sur le bord du
  noyau de viabilité.
 */
class TychasticBubbleBorderPickStrategy final : public TychasticControlPickStrategy {
public:
    TychasticBubbleBorderPickStrategy(int strategyIndex, const TrajectoryParametersManager *, const Grid *);
    OptionalCu pickControl(const OptionalCu &previousStrategyResult,
                           TychasticControlPickCriteria &criteria) override final;
    const std::string &getName() const override final;
    virtual ~TychasticBubbleBorderPickStrategy() = default;
private:
    static const std::string name;
    Bubble bubble;
};

#endif /* TYCHASTICCONTROLPICKSTRATEGIES_H */
