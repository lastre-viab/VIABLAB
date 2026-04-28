/*
 * defs.h
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
 *  Created on: reated on: 27 jul. 2013
 *      Author: Anna DESILLES
 */

#ifndef DEFS_H_

#define DEFS_H_

// #include "CropData.h"
// #include "CycleData.h"

#include <list>
#include <sstream>
#include <stdlib.h>
#include <valarray>
#include <vector>

// #include <boost/thread.hpp>
#include "boost/dynamic_bitset/dynamic_bitset.hpp"
#include <fstream>
#include <iostream>

#include "boost/foreach.hpp"
#include "boost/property_tree/json_parser.hpp"
#include "boost/property_tree/ptree.hpp"
using namespace boost::property_tree;

#include "spdlog/sinks/basic_file_sink.h"
#include "spdlog/sinks/stdout_color_sinks.h" // or "../stdout_sinks.h" if no colors needed
#include "spdlog/spdlog.h"
 #include <functional>

#include <ctime>
#include <map>
#include <numeric> // pour std::iota
#include <omp.h>
#include <sys/time.h>
#include <sys/types.h>

#define NB_MAX_TRAJ_ITER 100000000
#define NB_MAX_TRAJ_SIMULATIONS 50

#define DEV_PRINT 1
#include "Enums.h"
#include "utilities.h"

using namespace std;

/*!
 * \struct discretImageSet
 * \brief  Structure pour le stockage d'une image discrete d'un point;
 *
 * Decrit sous forme compacte l'ensemble image d'un point de grille par la
 * dynamique discrete
 *  \var nbImageCells repr�sente le nombre de mailles de la grille, hors
 * doublons, contenues dans l'image
 *  \var tabImageCells est une tableau d'entiers  de taille nbImageCells
 * repr�sentant les num�ros de mailles contenues dans l'image
 *  \var tabImageControls repr�sente une permutation de l'ensemble
 *  d'entiers num�ros des controles {0,... nbTotalControls-1} dans l'ordre
 * d'association avec les mailles images
 *  \var tabCellsEntrees donne le debut dans tabImageControls de la  liste des
 * controles associes  aux  cellules
 * \see viabiHJB
 */
struct discretImageSet
  {
  unsigned long long int nbImageCells;
  unsigned long long int *tabImageCells;
  unsigned long long int *tabImageControls;
  unsigned long long int *tabCellEntrees;
  };

struct discretImageSet_DD
  {
  unsigned long long int nbImagePoints;
  unsigned long long int *tabImagePoints;
  unsigned long long int *tabImageControls;
  unsigned long long int *tabPointsEntrees;
  };

struct discretImageSet_simple
  {
  unsigned long long int nbImageCells;
  unsigned long long int *tabImageCells;
  };

/*!
 * \struct triple
 *
 * \brief  d�finit un triplet pour un ant�c�dant de cellule
 * \see ViabiHJB
 */
struct triple
  {
  unsigned long long int posX; /*!< num�ro de l'�tat ant�c�dant */
  unsigned long long int
      posU; /*!< num�roe de contr�le tel que x-rho*f(x,u) est dans la cellule*/
  double rho; /*!< le pas de temps */
  };

inline bool equal_triple(triple t1, triple t2)
  {
  return ((t1.posX == t2.posX) & (t1.posU == t2.posU) & (t1.rho == t2.rho));
  }
inline bool compare_triple_first(triple t1, triple t2)
  {
  return ((t1.posX < t2.posX));
  }
inline bool compare_triple_second(triple t1, triple t2)
  {
  return ((t1.posU < t2.posU));
  }
/*!
 * \struct imageCell
 *
 * \brief  d�finit une structure pour une cellule image
 */
struct imageCell
  {
  unsigned long long int cellNum; /*!< numero de cellule */
  double
      minVal; /*!< la valeur optimale de la fonction valeur sur cette cellule*/
  list<triple>
      R_optimal; /*!< la liste des ant�c�dents minimisant la fonction valeur*/
  list<triple> R_viable; /*!< la liste des ant�c�dents viables arrivant dans
                            cette  cellule */
  };

/*!
 * \struct imageCellsList structure utilis�e pour stocker de fa�on compacte
 *  l'image  par la dynamique discrete  d'une couche en cours.
 *
 *  Les valeurs minNum et maxNum servent � gagner du temps pour les algorithmes
 * d'insetion.
 */
struct imageCellsList
  {
  list<imageCell> cellsList;
  int minNum;
  int maxNum;
  };

struct imageCellsIntList
  {
  list<unsigned long long int> cellsList;
  unsigned long long int minNum;
  unsigned long long int maxNum;
  };

/*!
 * \struct imageCell
 *
 * \brief  d�finit une structure pour une cellule image
 */
struct imagePoint
  {
  unsigned long long int PointNum; /*!< numero de cellule */
  double
      minVal; /*!< la valeur optimale de la fonction valeur sur cette cellule*/
  list<triple>
      R_optimal; /*!< la liste des ant�c�dents minimisant la fonction valeur*/
  list<triple> R_viable; /*!< la liste des ant�c�dents viables arrivant dans
                            cette  cellule */
  };
struct imageTempPoint
  {
  unsigned long long int PointNum; /*!< numero de cellule */
  double
      minVal; /*!< la valeur optimale de la fonction valeur sur cette cellule*/
  list<triple>
      *R_optimal; /*!< la liste des ant�c�dents minimisant la fonction valeur*/
  list<triple> *R_viable; /*!< la liste des ant�c�dents viables arrivant dans
                             cette  cellule */
  };

/*!
 * \struct imageCellsList structure utilis�e pour stocker de fa�on compacte
 *  l'image  par la dynamique discrete  d'une couche en cours.
 *
 *  Les valeurs minNum et maxNum servent � gagner du temps pour les algorithmes
 * d'insetion.
 */
struct imagePointsList
  {
  std::list<imagePoint> *pointsList;
  unsigned long long int minNum;
  unsigned long long int maxNum;
  };

struct imagePointsIntList
  {
  list<unsigned long long int> pointsList;
  unsigned long long int minNum;
  unsigned long long int maxNum;
  };

/*!
 * op�rateur de comparaison des ellules images  utilis� pour
 * trier les listes  de cellules  images par le num�ro de cellule (le premier
 * entier)
 *
 */
inline bool imageCellCompare(imageCell c1, imageCell c2)
  {
  return (c1.cellNum <= c2.cellNum);
  }

/***********************************************************************************************
 *  Fonctions abstraites g�n�riques utiles � toutes les classes
 ***********************************************************************************************/

/*! \brief Transformation  de num�rotation inverse
 *  � partir du num�ro du point dans la num�rotation alpha-num�rique
 * on obient le vecteur des coordonn�s enti�res
 */
inline void numToIntCoords_gen(unsigned long long int num,
                               unsigned long long int dim,
                               const unsigned long long int *nbPoints,
                               unsigned long long int *res)
  {
  int temp = num;

  for (int d = dim - 1; d >= 0; d--)
    {
    const unsigned long long int resd =
        temp % nbPoints[d]; // coordonn�es enti�res du point
    temp /= nbPoints[d];
    res[d] = resd;
    }
  }

/*! \brief Transformation  de num�rotation alpha-num�rique
 *  *  � partir du num�ro des coordonn�es enti�res du points dans la grille
 *  sont num�ro
 */
inline unsigned long long int
intCoordsToNum_gen(unsigned long long int dim,
                   const unsigned long long int *nbPoints,
                   unsigned long long int *coords)
  {
  unsigned long long int res = coords[0];

  for (unsigned long long int i = 0; i < dim - 1; i++)
    {
    res = res * nbPoints[i + 1] + coords[i + 1];
    }
  return res;
  }

/*!
 *  le type intPair est d�fini pour le calcul des images
 *   discretes d'un point
 */
typedef pair<int, int> intPair;

/*!
 * op�rateur de comparaison des paires d'entiers utilis� pour
 * trier les listes  de cellules par le num�ro de cellule (le premier entier)
 * le second entier est le num�ro de contr�le  correspondant
 */
inline bool pairCompare(intPair p1, intPair p2)
  {
  return (p1.first <= p2.first);
  }

typedef pair<unsigned long long int, unsigned long long int> uintPair;

/*!
 * op�rateur de comparaison des paires d'entiers utilis� pour
 * trier les listes  de cellules par le num�ro de cellule (le premier entier)
 * le second entier est le num�ro de contr�le  correspondant
 */
inline bool upairCompare(uintPair p1, uintPair p2)
  {
  return (p1.first <= p2.first);
  }

template <typename T>
void printVector(const T *vect, unsigned long long int dim)
  {
  for (unsigned long long int k = 0; k < dim; k++)
    {
    cout << " " << vect[k];
    }
  cout << endl;
  }

template <typename T>
void logVector(const string &msg, const T *vect, unsigned long long int dim)
  {
  ostringstream os;
  os << msg;
  for (unsigned long long int k = 0; k < dim; k++)
    {
    os << " " << vect[k];
    }
  spdlog::debug(os.str().c_str());
  }

inline void logDynBit(const string &msg, const boost::dynamic_bitset<> *masque)
  {
  ostringstream os;
  os << msg;
  os << " " << *masque;
  spdlog::debug(os.str().c_str());
  }

#endif /* DEFS_H_ */
