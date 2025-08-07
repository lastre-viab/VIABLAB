/*
 * SysDyn.h
 *
 *    VIABLAB : a numerical library for Mathematical Viability Computations
 *    Copyright (C) <2020>  <Anna DESILLES, LASTRE>
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU Affero General Public License as
 *   published by the Free Software Foundation, either version 3 of the
 *   License, or (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU Affero General Public License for more details.
 *
 *   You should have received a copy of the GNU Affero General Public License
 *   along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
 *  Created on: reated on: 27 sept. 2013
 *      Author: Anna DESILLES
 */

#ifndef SYSDYN_H_
#define SYSDYN_H_

#include "defs.h"
#include "Params.h"
#include "Grid.h"

/*!
 * \class SysDyn
 * \brief La classe qui regroupe tous les ?����l?����ments d'un syst?����me dynamique
 */
class SysDyn
    {
public:
    /*!
     * Constructeur vide
     */
    SysDyn();
    /*!
     * Destructeur
     */
    virtual ~SysDyn();

    SysDyn(const SysDyn &) = delete;
    SysDyn(const SysDyn &&) = delete;
    SysDyn& operator=(const SysDyn &) = delete;
    SysDyn& operator=(SysDyn&& data);

    /*!
     * \brief Fonction  d?����finissant la dynamique
     *
     * Cette  fonction est d?����finie ici comme pointeur  pour permetre
     *  de la d?����finir  de mani?����re ext?����rieure ?���� la  classe
     *  Ainsi la dynamique est un param?����tre de la classe
     *
     * @param[in] x : la variable d'?����tat
     * @param[in] u : la variable  de contr?����le
     * @param[out] res : le r?����sultat
     */
    void (*dynamics)(const double*, const double*, double*);
    void (*dynamics_tych)(const double*, const double*, const double*, double*);
    void (*dynamics_fd)(const unsigned long long int*, const unsigned long long int*, unsigned long long int*);
    void (*dynamics_tych_fd)(const unsigned long long int*, const unsigned long long int*, const unsigned long long int*, unsigned long long int*);

    /*!
     * \brief Pointeur sur l'une des fonctions membres d?����finissant une m?����thode de discr?����tisation
     *
     * la valeur de ce pointeur est initialis?����e dans le constructeur  en fonction
     *  des param?����tres d'utilisateur d?����finissant ses choix par rapport ?���� la discr?����tisation du syst?����me
     *
     * @param[in] x : la variable d'?����tat
     * @param[in] u : la variable  de contr?����le
     * @param[out] res : le r?����sultat
     */
    void (SysDyn::*discretDynamics)(const double*, const double*, double*, double rho) const;
    void (SysDyn::*discretDynamics_tych)(const double*, const double*, const double*, double*, double rho) const;

    /*!
     * \brief Pointeur sur l'une des fonctions membres d?����finissant
     * une m?����thode de calcul de la constante de Lipshitz
     *
     * la valeur de ce pointeur est initialis?����e dans le constructeur  en fonction
     *  des param?����tres d'utilisateur d?����finissant ses choix par rapport au calcul de
     *  la constante de Lipshitz
     *
     * @param[in] x : la variable d'?����tat
     *
     *\see LC, computeLC, jacobian
     */
    double (SysDyn::*calcul_L)(const double *x) const;
    /*!
     * \brief Pointeur sur l'une des fonctions membres d?����finissant
     * une m?����thode de calcul de la constante M de brne  de la dynamique
     *
     * la valeur de ce pointeur est initialis?����e dans le constructeur  en fonction
     *  des param?����tres d'utilisateur d?����finissant ses choix par rapport au calcul de
     *  la constante M
     *
     * @param[in] x : la variable d'?����tat
     *
     *\see M, computeM, localDynBounds
     */
    double (SysDyn::*calcul_M)(const double *x) const;

    /*!
     * \brief Cette fonction d?����finit les contraintes sur le contr?����le, elle correspond ?���� la d?����finition  de U(x)
     *
     * Attention!  Cette fonction est d?����finie par l'utilisateur  dans une fichier "dataMonModele.h"
     * Ici dans le code c'est un param?����tre dont ne conna?����t que la "signature" : les types des arguments  et de la sortie
     *
     * @param[in] x l'?����tat x
     * @param[in] u le contr?����le
     * @return valeur r?����elle si (x,u)  est dans le graphe de U(x), \f$ +\infty\f$ sinon
     */
    double (*constraintsXU)(const double *x, const double *u);
    double (*constraintsXV_tych)(const double *x, const double *v);
    double (*constraintsXU_fd)(const unsigned long long int *x, const unsigned long long int *u);
    double (*controlEligibilityForTraj_fd)(const unsigned long long int *x, const unsigned long long int *u, const unsigned long long int *previousU);
    /*!
     * \brief  Cette fonction d?����finit les contraintes sur l'?����tat, elle correspond a la d?����finition k(x)
     *
     * Attention!  Cette fonction est d?����finie par l'utilisateur  dans une fichier "dataMonModele.h"
     * Ici dans le code c'est un param?����tre dont ne conna?����t que la "signature" : les types des arguments  et de la sortie
     *
     * @param[in] x l'?����tat x
     * @return valeur r?����elle si\f$x\in K\f$, \f$ +\infty\f$ sinon
     */
    double (*constraintsX)(const double*);
    double (*constraintsX_fd)(const unsigned long long int*);

    /*!
     * \brief  Cette fonction d?����finit les contraintes sur l'?����tat, elle correspond a la d?����finition k(x)
     *
     * Attention!  Cette fonction est d?����finie par l'utilisateur  dans une fichier "dataMonModele.h"
     * Ici dans le code c'est un param?����tre dont ne conna?����t que la "signature" : les types des arguments  et de la sortie
     *
     * @param[in] x l'?����tat x
     * @return valeur r?����elle si\f$x\in K\f$, \f$ +\infty\f$ sinon
     */
    double (*dynConstraintsForTraj)(const double*, double*);

    /*!
     * \brief  Cette fonction d?����finit la cible, elle correspond a la d?����finition de c(x)
     *
     * Attention!  Cette fonction est d?����finie par l'utilisateur  dans une fichier "dataMonModele.h"
     * Ici dans le code c'est un param?����tre dont ne conna?����t que la "signature" : les types des arguments  et de la sortie
     *
     * @param[in] x l'?����tat x
     * @return valeur r?����elle si\f$x\in C\f$, \f$ +\infty\f$ sinon
     */
    double (*target)(const double*);
    double (*target_fd)(const unsigned long long int*);

    /*!
     * \brief  Cette fonction correspond ?���� la d?����finition de l(x,u) dans le cas d'un probl?����me d'optimisation
     * avec la fonction objectif de type \f$ J(x,u)=\int exp(m(x(\tau),u(\tau)))l(x(\tau),,u(\tau)) d\tau\f$
     *
     * Attention!  Cette fonction est d?����finie par l'utilisateur  dans une fichier "dataMonModele.h"
     * Ici dans le code c'est un param?����tre dont ne conna?����t que la "signature" : les types des arguments  et de la sortie
     *
     * @param[in] x l'?����tat x
     * @param[in] u contr?����le
     * @return valeur r?����elle si\f$(x,u)\in Dom(l)\f$, \f$ +\infty\f$ sinon
     */
    double (*lFunc)(const double *x, const double *u);
    double (*lFunc_tych)(const double *x, const double *u, const double *v);
    double (*lFunc_fd)(const unsigned long long int *x, const unsigned long long int *u);
    double (*lFunc_tych_fd)(const unsigned long long int *x,
                            const unsigned long long int *u,
                            const unsigned long long int *v);

    /*!
     * \brief  Cette fonction correspond a la d?����finition de mu(x,u)  qui definit les contraintes
     * sur le controle sous forme epigraphique
     * Attention!  Cette fonction est d?����finie par l'utilisateur  dans une fichier "dataMonModele.h"
     * Ici dans le code c'est un param?����tre dont ne conna?����t que la "signature" : les types des arguments  et de la sortie
     *
     * @param[in] x l'?����tat x
     * @param[in] u contr?����le
     * @return valeur r?����elle si\f$(x,u)\in Dom(l)\f$, \f$ +\infty\f$ sinon
     */

    double (*muFunc_fd)(const unsigned long long int *x, const unsigned long long int *u);

    /*!
     * \brief  Cette fonction correspond ?���� la d?����finition de m(x,u) dans le cas d'un probl?����me d'optimisation
     * avec la fonction objectif de type \f$ J(x,u)=\int exp(m(x(\tau),u(\tau)))l(x(\tau),,u(\tau)) d\tau\f$
     *
     * Attention!  Cette fonction est d?����finie par l'utilisateur  dans une fichier "dataMonModele.h"
     * Ici dans le code c'est un param?����tre dont ne conna?����t que la "signature" : les types des arguments  et de la sortie
     *
     * @param[in] x l'?����tat x
     * @param[in] u contr?����le
     * @return valeur r?����elle si\f$(x,u)\in Dom(l)\f$, \f$ +\infty\f$ sinon
     */
    double (*mFunc)(const double *x, const double *u);
    double (*mFunc_tych)(const double *x, const double *u, const double *v);
    /*!
     * \brief Constructeur principal
     *
     * Ce constructeur est appel?����  dans le main pour passer  ?���� la classe les diff?����rents param?����tres
     *  du syst?����me dynamique ainsi  que  des choix de l'utilisateur concernant
     *  la discr?����tisation et le calcul des constantes
     *
     * @param SP param?����tres du syst?����me, regroup?����s dans une structure
     * @param ds : dimension d'?����tat
     * @param cp : param?����tres de la variable de controle
     * @param refGrid : r?����f?����rence ?���� une grille choisie pour l'?����tude num?����rique du syst?����me
     *
     * \see systemParams, controlParams, Grid
     */
    SysDyn(const systemParams &SP, int, const controlParams &cp, Grid *refGrid);

    /*
     * Dans la conception actuelle la classe sysDyn est la seule ?���� poss?����der l'information sur les controles: la dimension, le
     * nombre de points de dicr?����tisation,
     * la nature g?����om?����trique etc
     *
     * Elle poss?����de donc les m?����thodes d'acc?����s correspondantes pour  fournir aux autres classes qui en auront besoin
     * ces param?����tres
     */

    /*!
     * M?����thode d'acc?����s: renvoie les bornes inf du pav?���� contenant les  contr?����les
     * @return limInfC
     */
    double* getLimInfC();

    /*!
     * M?����thode d'acc?����s: renvoie les bornes SUP du pav?���� contenant les  contr?����les
     * @return limSupC
     */

    double* getLimSupC();

    /*!
     * M?����thode d'acc?����s: renvoie le pas de discr?����tisation du pav?���� contenant les  contr?����les
     * @return stepC
     */

    double* getStepC();
    /*!
     * M?����thode d'acc?����s: renvoie la dimension de l'espace des  contr?����les
     * @return dimC
     */
    unsigned long long int getDimC() const;
    /*!
     * M?����thode d'acc?����s: renvoie les nb de points par axe des  contr?����les
     * @return NbPointsC
     */
    unsigned long long int* getNbPointsC();
        
    /*!
     * M?����thode d'acc?����s: renvoie le pointeur sur le tableau qui stocke les coordonn?����es r?����elles des contr?����les
     * @return  ControlCoords
     */
    double **getControlCoords() const;

    /*!
     * M?����thode d'acc?����s: renvoie le pointeur sur le tableau qui stocke les coordonn?����es r?����elles des contr?����les
     * @return  ControlCoords
     */
        
    unsigned long long int** getControlIntCoords();
    /*!
     * M?����thode d'acc?����s: renvoie le nombre total de points de discr?����tisation de contr?����les
     * @return  ControlCoords
     */
    unsigned long long int getTotalNbPointsC() const;

    /*!
     * M?����thode d'acc?����s: renvoie les bornes inf du pav?���� contenant les  contr?����les
     * @return limInfC
     */
    double* getLimInfTy();

    /*!
     * M?����thode d'acc?����s: renvoie les bornes SUP du pav?���� contenant les  contr?����les
     * @return limSupC
     */

    double* getLimSupTy();

    /*!
     * M?����thode d'acc?����s: renvoie le pas de discr?����tisation du pav?���� contenant les  contr?����les
     * @return stepC
     */

    double* getStepTy();
    /*!
     * M?����thode d'acc?����s: renvoie la dimension de l'espace des  contr?����les
     * @return dimC
     */
    unsigned long long int getDimTy() const;
    /*!
     * M?����thode d'acc?����s: renvoie les nb de points par axe des  contr?����les
     * @return NbPointsC
     */
    unsigned long long int* getNbPointsTy();
    /*!
     * M?����thode d'acc?����s: renvoie le pointeur sur le tableau qui stocke les coordonn?����es r?����elles des contr?����les
     * @return  ControlCoords
     */
    double** getTychCoords() const;
    /*!
     * M?����thode d'acc?����s: renvoie le pointeur sur le tableau qui stocke les coordonn?����es r?����elles des contr?����les
     * @return  ControlCoords
     */
    unsigned long long int** getTychIntCoords();
    /*!
     * M?����thode d'acc?����s: renvoie le nombre total de points de discr?����tisation de contr?����les
     * @return  ControlCoords
     */
    unsigned long long int getTotalNbPointsTy() const;

    /*   ***************************************************************
     *  Fin des m?����thodes d'acc?����s aux param?����tres de controles
     *
     ******************************************************************** */

    /*!
     * \brief M?����thode qui calcul le pas de temps local qunad l'option correspondante a ?����t?���� choisie
     *
     * @param x variable d'?����tat
     * @param u variable de contr?����le
     */
    double calculRho_local(const double *x) const;

    /*!
     * m?����thode d'acc?����s au param?����tre qui d?����finit  si le pas de temps
     *  des algorithmes doit ?����tre local ou global (par rapport ?���� la grille d'?����tat)
     */
    bool isTimeStepGlobal();

    int getFDDynType();
    string getRetroFileName();

    /*!
     *  \brief Méthode  qui retourne le type de dynamique
     * @return dynType
     */
    DynType getDynType();

    void setDynamicsForward();
    void setDynamicsBackward();

    enum PointStatus : unsigned char {
        VALID_TRAJECTORY_POINT,
        OUTSIDE_DOMAIN,
        OUTSIDE_CONSTRAINTS,
        OUTSIDE_GRID,
    };
    PointStatus checkKernelRelation(double *point) const;

    bool isViableControl(const double *x, const double *u, double *image, double rho) const;
    bool isViableControl_tych(const double *x, const double *u, const double *v, double *image, double rho) const;
    bool isViableGuaranteedControl(const double *x, const double *u, double rho) const;

    void getTychasticImage(const double *x, const double *u, const double *v, double *imageVect, double rho) const;
        
    const Grid *getGrid() const;
    int getDim() const;
private:

    void initializeMethods(const systemParams &SP);
        
    double calculL_local_num(const double *x) const;
    double calculMF_local_num(const double *x) const;
    double calculL_local_ana(const double *x) const;
    double calculMF_local_ana(const double *x) const;

    double calculL_local_num_tych(const double *x) const;
    double calculMF_local_num_tych(const double *x) const;
    double calculL_local_ana_tych(const double *x) const;

    double returnL_local_ana(const double *x) const;

    double returnMF_local_ana(const double *x) const;

    /*!
     * Méthode  qui sert  dans le cas où la dynamique est déjà discrète en temps
     * mais continue en espace
     */
    void FDiscret(const double *x, const double *u, double *res, double rho) const;
    /*!
     * Méthodes  qui seront utilisées dans le cas où
     * la dynamique est continue en temps et en espace
     */
    void FDiscretEuler(const double *x, const double *u, double *res, double rho) const;
    void FDiscretRK2(const double *x, const double *u, double *res, double rho) const;
    void FDiscretRK4(const double *x, const double *u, double *res, double rho) const;

    void FDiscret_tych(const double *x, const double *u, const double *v, double *res, double rho) const;
    /*!
     * Méthodes  qui seront utilisées dans le cas où
     * la dynamique est continue en temps et en espace
     */
    void FDiscretEuler_tych(const double *x, const double *u, const double *v, double *res, double rho) const;
    void FDiscretRK2_tych(const double *x, const double *u, const double *v, double *res, double rho) const;
    void FDiscretRK4_tych(const double *x, const double *u, const double *v, double *res, double rho) const;


    int fd_dyn_type;
    double timeHorizon;
    string retroFileName;
    int discretisation;
    DynType dynType;
    int dimS;
    int dimC;
    double *limInfC;
    double *limSupC;
    double *stepC;
    unsigned long long int *nbPointsC;
    double **controlCoords;
    unsigned long long int **controlIntCoords;
    unsigned long long int totalNbPointsC;

    unsigned long long int *nbPointsTy;
    int dimTy;
    double *limInfTy;
    double *limSupTy;
    double *stepTy;
    double **tychCoords;
    unsigned long long int **tychIntCoords;
    unsigned long long int totalNbPointsTych;
    bool isTychastic;

    double L;  // constante de Lipschitz
    double MF;  // majoration de la norme de la dynamique

    double lfunc_L;  // constante de Lipschitz
    double lfunc_MF;  // majoration de la norme de la dynamique

    double *image, *FXmoinsH, *xTemp, *FXplusH;
    Grid *grid;
    bool globalTimeStep;
    ComputeMethod computeMF;
    ComputeMethod computeLC;
    double **jacob;
    void (*localDynBounds)(const double *x, double *res);
    void (*jacobian)(const double *x, const double *u, double **jacob);
    void (*jacobian_tych)(const double *x, const double *u, const double *v, double **jacob);
    double dynSignFactor;
};

#endif /* SYSDYN_H_ */
