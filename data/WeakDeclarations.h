/*
 * WeakDeclarations.h
 *
 *  Created on: 27 juil. 2021
 *      Author: adesi
 */

#ifndef DATA_WEAKDECLARATIONS_H_
#define DATA_WEAKDECLARATIONS_H_

extern double level;
extern int sortieOKinf[];
extern int sortieOKsup[];
extern double l_Lip;
extern double l_max;
extern int saveSubLevel;

double  constraintsXU_fd( unsigned long long int * x, unsigned long long int * u ) __attribute__((weak));
void dynamics_tych_fd(unsigned long long int  * x, unsigned long long int *u,unsigned long long int *v, unsigned long long int * image) __attribute__((weak));
double  constraintsX_fd( unsigned long long int * x ) __attribute__((weak));

void dynamics_fd(unsigned long long int  * x, unsigned long long int *u, unsigned long long int * image) __attribute__((weak));
double dynConstraintsForTraj(double * x, double * image) __attribute__((weak));

double dynConstraintsForTraj_default(double * x, double * image);

double dynConstraintsForTraj_default(double * x, double * image)
{
	return 1.0;
}

#endif /* DATA_WEAKDECLARATIONS_H_ */
