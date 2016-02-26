// written by Szymon Winczewski

#ifndef cuda_central_force_h
#define cuda_central_force_h

#include <iostream>
#include "defines.h"

typedef struct
{
    bool first_second;
    int atom_j;
    double force[3];
    double r_ij;
    double r_ij_dir[3];
} central_force;


void allocateCentralForces(int number_of_atoms, int max_number_of_central_forces,
                           int *&number_of_central_forces, central_force **&central_forces);

void deallocateCentralForces(int number_of_atoms,
                             int *&number_of_central_forces, central_force **&central_forces);

#endif
