/*
 * WeakDeclarations.h
 *
 *  Created on: 27 juil. 2021
 *      Author: adesi
 */

#ifndef DATA_WEAKDECLARATIONS_H_
#define DATA_WEAKDECLARATIONS_H_

#include "defs.h"
extern double level;
extern int sortieOKinf[];
extern int sortieOKsup[];
extern double l_Lip;
extern double l_max;
extern int saveSubLevel;
unsigned long long int *initPoints_fd;
double *initValues_fd;
double *initValues;

void loadModelData(ParametersManager *PM) __attribute__((weak));

__attribute__((weak)) double constraintsXU_fd(unsigned long long int *x,
	unsigned long long int *u);
void dynamics_tych_fd(unsigned long long int *x, unsigned long long int *u,
	unsigned long long int *v, unsigned long long int *image)
		__attribute__((weak));
double __attribute__((weak)) constraintsX_fd(unsigned long long int *x);
void dynamics_fd(unsigned long long int *x, unsigned long long int *u,
	unsigned long long int *image) __attribute__((weak));
double dynConstraintsForTraj(double *x, double *image) __attribute__((weak));
double constraintsXUY_fd(unsigned long long int *x, unsigned long long int *u)
	__attribute__((weak));
double controlEligibilityForTraj_fd(unsigned long long int *x,
	unsigned long long int *u, unsigned long long int *previousU)
		__attribute__((weak));
double target_fd(unsigned long long int *x) __attribute__((weak));
double l_fd(unsigned long long int *x, unsigned long long int *u)
	__attribute__((weak));
double l_tych_fd(unsigned long long int *x, unsigned long long int *u,
	unsigned long long int *v) __attribute__((weak));

double dynConstraintsForTraj_default(double *x, double *image);

void dynamics(double *x, double *u, double *image) __attribute__((weak));
void dynamics_tych(double *x, double *u, double *v, double *image) __attribute__((weak));
void jacobian(double *x, double *u, double **jacob) __attribute__((weak));
void jacobian_tych(double *x, double *u, double *v, double **jacob) __attribute__((weak));
void localDynBounds(double *x, double *res) __attribute__((weak));
double constraintsXU(double *x, double *u) __attribute__((weak));
double constraintsX(double *x) __attribute__((weak));
double target(double *x) __attribute__((weak));
double l(double *x, double *u) __attribute__((weak));
double m(double *x, double *u) __attribute__((weak));
double l_tych(double *x, double *u, double *v) __attribute__((weak));
double m_tych(double *x, double *u, double *v) __attribute__((weak));
void dynamics_hybrid(double *x, double *u, double *image) __attribute__((weak));
void postProcess(ParametersManager *PM) __attribute__((weak));

inline double l_default(double *x, double *u)
    {
    return 1.0;
    }
inline double m_default(double *x, double *u)
    {
    return 0.0;
    }
inline double l_tych_default(double *x, double *u, double *v)
    {
    return 1.0;
    }
inline double m_tych_default(double *x, double *u, double *v)
    {
    return 0.0;
    }

void dynamics_hybrid_default(double *x, double *u, double *image)
    {

    }
void postProcess_default(ParametersManager *PM)
    {

    }

inline double target_default(double *x)
    {
    return PLUS_INF;
    }

inline double constraintsX_default(double *x)
    {
    return 0.0;
    }

inline double constraintsXU_default(double *x, double *u)
    {
    return 0.0;
    }
inline void localDynBounds_default(double *x, double *res)
    {

    }

inline double constraintsXU_fd_default(unsigned long long int *x,
	unsigned long long int *u)
    {
    return 0.0;
    }
inline void jacobian_default(double *x, double *u, double **jacob)
    {

    }
inline void jacobian_tych_default(double *x, double *u, double *v, double **jacob)
    {

    }
inline void dynamics_default(double *x, double *u, double *image)
    {

    }
inline void dynamics_tych_default(double *x, double *u, double *v, double *image)
    {

    }

double dynConstraintsForTraj_default(double *x, double *image)
    {
    return 1.0;
    }

void loadModelData_void()
    {

    }
inline double target_fd_default(unsigned long long int *x)
    {
    return 0.0;
    }

inline double l_fd_default(unsigned long long int *x, unsigned long long int *u)
    {
    return 0.0;
    }

inline double l_fd_tych_default(unsigned long long int *x,
	unsigned long long int *u, unsigned long long int *v)
    {
    return 0.0;
    }

inline double controlEligibilityForTraj_fd_default(unsigned long long int *x,
	unsigned long long int *u, unsigned long long int *previousU)
    {
    return 1.0;
    }
inline double constraintsXUY_fd_default(unsigned long long int *x,
	unsigned long long int *u)
    {
    return -PLUS_INF;
    }

#endif /* DATA_WEAKDECLARATIONS_H_ */
