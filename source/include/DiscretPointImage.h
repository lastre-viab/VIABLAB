/*
 * DiscretPointImage.h
 *
 *  Created on: 2 nov. 2024
 *      Author: adesi
 */

#ifndef DISCRETPOINTIMAGE_H_
#define DISCRETPOINTIMAGE_H_
#include "SysDyn.h"
#include "Grid.h"

class DiscretPointImage
    {
public:
    DiscretPointImage();
    DiscretPointImage(unsigned long long int posX, int dim, SysDyn *sd, Grid *gd);
    void Build();
    void setPointNumAndRebuild(unsigned long long int pointNum);
    bool testConstraintesForCell(unsigned long long int numCell);
    bool IsBuilt();
    std::map<unsigned long long int, list<unsigned long long int> >* GetImageCellsWithControls();
    virtual
    ~DiscretPointImage();

private:
    unsigned long long int posX;
    bool _isBuilt;
    double *pointDoubleCoords;
    unsigned long long int *pointIntCoords;
    int dim;
    SysDyn *dynsys;
    Grid *grid;

    std::map<unsigned long long int, list<unsigned long long int> > *imageCellsWithControls;
    };

#endif /* DISCRETPOINTIMAGE_H_ */
