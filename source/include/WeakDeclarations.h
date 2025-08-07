/*
 * WeakDeclarations.h
 *
 *  Created on: 27 juil. 2021
 *      Author: adesi
 */

#ifndef DATA_WEAKDECLARATIONS_H_
#define DATA_WEAKDECLARATIONS_H_

#include "defs.h"

class ParametersManager;

void loadModelData(const ParametersManager *PM);

double constraintsXU_fd(const unsigned long long int *x,
                        const unsigned long long int *u);
void dynamics_tych_fd(const unsigned long long int *x, const unsigned long long int *u,
                      const unsigned long long int *v, unsigned long long int *image);
double constraintsX_fd(const unsigned long long int *x);
void dynamics_fd(const unsigned long long int *x, const unsigned long long int *u,
                 unsigned long long int *image);
double dynConstraintsForTraj(const double *x, double *image);
double constraintsXUY_fd(const unsigned long long int *x, const unsigned long long int *u);
double controlEligibilityForTraj_fd(const unsigned long long int *x,
                                    const unsigned long long int *u, const unsigned long long int *previousU);
double target_fd(const unsigned long long int *x);
double l_fd(const unsigned long long int *x, const unsigned long long int *u);
double l_tych_fd(const unsigned long long int *x, const unsigned long long int *u,
                 const unsigned long long int *v);

void dynamics(const double *x, const double *u, double *image);
void dynamics_tych(const double *x, const double *u, const double *v, double *image);
void jacobian(const double *x, const double *u, double **jacob);
void jacobian_tych(const double *x, const double *u, const double *v, double **jacob);
void localDynBounds(const double *x, double *res);
double constraintsXU(const double *x, const double *u);
double constraintsX(const double *x);
double target(const double *x);
double l(const double *x, const double *u);
double m(const double *x, const double *u);
double l_tych(const double *x, const double *u, const double *v);
double m_tych(const double *x, const double *u, const double *v);
void dynamics_hybrid(const double *x, const double *u, double *image);
void postProcess(const ParametersManager *PM);
double constraintsXV_tych(const double *x, const double *v);

#endif /* DATA_WEAKDECLARATIONS_H_ */
