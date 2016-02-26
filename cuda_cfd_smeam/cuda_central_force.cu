// written by Szymon Winczewski

#include "cuda_central_force.h"

void allocateCentralForces(int number_of_atoms,
                           int max_number_of_central_forces,
                           int *&number_of_central_forces,
                           central_force **&central_forces)
{
    number_of_central_forces = new int [number_of_atoms];
    central_forces = new central_force *[number_of_atoms];
    for (int i = 0; i < number_of_atoms; i++)
        central_forces[i] = new central_force [max_number_of_central_forces];

    #ifdef CFD_SMEAM_DEBUG
    std::cout << "memory usage:" << std::endl;
    std::cout << "  number_of_atoms              = " << number_of_atoms << std::endl;
    std::cout << "  max_number_of_central_forces = " << max_number_of_central_forces << std::endl;
    std::cout << "  number_of_central_forces     = " << number_of_atoms * sizeof(int) << " Bytes" << std::endl;
    std::cout << "  central_forces               = " << number_of_atoms * max_number_of_central_forces * sizeof(central_force) << " Bytes" << std::endl;
    #endif
}


void deallocateCentralForces(int number_of_atoms,
                             int *&number_of_central_forces,
                             central_force **&central_forces)
{
    delete [] number_of_central_forces;
    for (int i = 0; i < number_of_atoms; i++)
        delete [] central_forces[i];
    delete [] central_forces;
}
