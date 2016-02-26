// written by Szymon Winczewski

#ifndef central_force_h
#define central_force_h

#include <iostream>
#include "cuda_central_force.h"

void allocateCentralForces(int number_of_atoms, int max_number_of_central_forces,
                           int *&number_of_central_forces, central_force **&central_forces);

void deallocateCentralForces(int number_of_atoms,
                             int *&number_of_central_forces, central_force **&central_forces);

#endif
