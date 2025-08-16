#ifndef CONTROLPICKSTRATEGIES_H
#define CONTROLPICKSTRATEGIES_H

#include "ViabiBitSetTrajectoryHelper.h"
#include "Bubble.h"
#include "ControlPickStrategy.h"

/*!
  Choisit le même que l'itération précédente (le contrôle renvoyé par getLastControlIndex de traj) en modifiant possiblement le pas de temps rho.
  Si celui-ci n'est pas viable, renvoie UNSATISFIED_STRATEGY
*/
class VIABLAB_LIBRARY_EXPORT HeavyPickStrategy final : public ControlPickStrategy {
public:
    HeavyPickStrategy(int initCu);
    OptionalCu pickControl(const OptionalCu &previousStrategyResult,
                           ControlPickCriteria &criteria) override final;
    const std::string &getName() const override final;
    virtual ~HeavyPickStrategy() = default;
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
class VIABLAB_LIBRARY_EXPORT FirstPickStrategyBitSet final : public ControlPickStrategy {
public:
    FirstPickStrategyBitSet(ViabiBitSetTrajectoryHelper *);
    OptionalCu pickControl(const OptionalCu &previousStrategyResult,
                           ControlPickCriteria &criteria) override final;
    const std::string &getName() const override final;
    virtual ~FirstPickStrategyBitSet() = default;
private:
    static const std::string name;
    ViabiBitSetTrajectoryHelper *viabiHelper;
};

/*!
  Renvoie le premier contrôle viable réel.

  L'ordre déterminant le premier contrôle de grille valide est celui des preferedControlIndexes de sysDyn.
  
  Si le premier n'est pas viable, renvoie UNSATISFIED_STRATEGY
 */
class VIABLAB_LIBRARY_EXPORT FirstOnlyPickStrategy final : public ControlPickStrategy {
public:
    FirstOnlyPickStrategy() = default;
    OptionalCu pickControl(const OptionalCu &previousStrategyResult,
                           ControlPickCriteria &criteria) override final;
    const std::string &getName() const override final;
    virtual ~FirstOnlyPickStrategy() = default;
private:
    static const std::string name;
};

/*!
  Réordonne les preferedControlIndexes de manière uniformément aléatoire

  Renvoie toujours UNSATISFIED_STRATEGY
 */
class VIABLAB_LIBRARY_EXPORT ShuffleStrategy final : public ControlPickStrategy {
public:
    ShuffleStrategy() = default;
    OptionalCu pickControl(const OptionalCu &previousStrategyResult,
                           ControlPickCriteria &criteria) override final;
    const std::string &getName() const override final;
    virtual ~ShuffleStrategy() = default;
private:
    static const std::string name;
};

/*!
  Réordonne les preferedControlIndexes selon l'ordre dicté par
  la fonction utilisateur controlWeight par ordre décroissant des poids.
 */
class VIABLAB_LIBRARY_EXPORT SortStrategy final : public ControlPickStrategy {
public:
    SortStrategy() = default;
    OptionalCu pickControl(const OptionalCu &previousStrategyResult,
                           ControlPickCriteria &criteria) override final;
    const std::string &getName() const override final;
    virtual ~SortStrategy() = default;
private:
    static const std::string name;
};

/*!
  Réapplique l'ordre par défaut des contrôles.
 */
class VIABLAB_LIBRARY_EXPORT ResetOrderStrategy final : public ControlPickStrategy {
public:
    ResetOrderStrategy() = default;
    OptionalCu pickControl(const OptionalCu &previousStrategyResult,
                           ControlPickCriteria &criteria) override final;
    const std::string &getName() const override final;
    virtual ~ResetOrderStrategy() = default;
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
class VIABLAB_LIBRARY_EXPORT ClosestPickStrategy final : public ControlPickStrategy {
public:
    ClosestPickStrategy() = default;
    OptionalCu pickControl(const OptionalCu &previousStrategyResult,
                           ControlPickCriteria &criteria) override final;
    const std::string &getName() const override final;
    virtual ~ClosestPickStrategy() = default;
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
class VIABLAB_LIBRARY_EXPORT PreferedPickStrategy final : public ControlPickStrategy {
public:
    PreferedPickStrategy() = default;
    OptionalCu pickControl(const OptionalCu &previousStrategyResult,
                           ControlPickCriteria &criteria) override final;
    const std::string &getName() const override final;
    virtual ~PreferedPickStrategy() = default;
private:
    static const std::string name;
};

/*!
  Choisit le contrôle obtenu au bord du noyau de viabilité par simulation
  de la trajectoire du ControlPicker passés en paramètre. Cette simulation
  n'est lancée que si la bulle créé dans le constructeur du picker touche
  elle-même un bord.
  
  Renvoie UNSATISFIED_STRATEGY si la bulle ne touche pas de bord ou si la
  simulation n'atteint pas de bord (ce qui peut être le cas si on dépasse
  NB_MAX_TRAJ_SIMULATIONS itérations, la simulation sort de la bulle sans
  jamais toucher de bord ou si les ControlPicker n'ont pas renvoyé de contrôle
  pour un point de la trajectoire simulé)

  Renvoie le contrôle choisi par le ControlPicker
 */
class VIABLAB_LIBRARY_EXPORT BubblePickStrategy final : public ControlPickStrategy {
public:
    BubblePickStrategy(int strategyIndex, const TrajectoryParametersManager *, const Grid *);
    OptionalCu pickControl(const OptionalCu &previousStrategyResult,
                           ControlPickCriteria &criteria) override final;
    const std::string &getName() const override final;
    virtual ~BubblePickStrategy() = default;
private:
    OptionalCu findBorderControl(ControlPickCriteria &criteria,
                                 unsigned long long int currentPos,
                                 unsigned long long int *currentPosIntCoords, double *doubleCoordsOnDiscreteTraj,
                                 double *imageVect);
    bool isValidControl(double *xCoordsDouble, OptionalCu opt, SysDyn &sysDyn, double *imageVect);
    
    static const std::string name;
    Bubble bubble;
    int realTimeStepsPerDiscreteStep;
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
class VIABLAB_LIBRARY_EXPORT SmoothPickStrategy final : public ControlPickStrategy {
public:
    SmoothPickStrategy(double maxAngleRadians);

    OptionalCu pickControl(const OptionalCu &previousStrategyResult,
                           ControlPickCriteria &criteria) override final;
    const std::string &getName() const override final;
    virtual ~SmoothPickStrategy() = default;
private:
    double angleBetween(const double *originA, const double *destinationA, const double *originB, const double *destinationB, int dim) const;
    double squaredDistanceBetween(const double *point1, const double *point2, int dim) const;

    OptionalCu smoothControlClosestTo(pickedControl &p,
                                      const double *previousToLastPosition, const double *lastPosition,
                                      ControlPickCriteria &criteria, double *imageVect);
    
    double maxAngleRadians;
    static const std::string name;
};

/*!
  Choisit l'indice du contrôle le plus proche du contrôle renvoyé par la fonction
  temporalControl de l'utilisateur.
 */
class VIABLAB_LIBRARY_EXPORT TemporalControlPickStrategy final : public ControlPickStrategy {
public:
    TemporalControlPickStrategy(const TrajectoryParametersManager *);
    OptionalCu pickControl(const OptionalCu &previousStrategyResult,
                           ControlPickCriteria &criteria) override final;
    const std::string &getName() const override final;

    int getClosestControlTo(double *userControl, double **controlCoords, unsigned long long int nbTotalC, int dimC) const;
    double squaredDistanceBetween(const double *point1, const double *point2, int dim) const;
    
    virtual ~TemporalControlPickStrategy() = default;
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
class VIABLAB_LIBRARY_EXPORT BubbleBorderPickStrategy final : public ControlPickStrategy {
public:
    BubbleBorderPickStrategy(int strategyIndex, const TrajectoryParametersManager *, const Grid *);
    OptionalCu pickControl(const OptionalCu &previousStrategyResult,
                           ControlPickCriteria &criteria) override final;
    const std::string &getName() const override final;
    virtual ~BubbleBorderPickStrategy() = default;
private:
    static const std::string name;
    Bubble bubble;
};
#endif /* CONTROLPICKSTRATEGIES_H */
