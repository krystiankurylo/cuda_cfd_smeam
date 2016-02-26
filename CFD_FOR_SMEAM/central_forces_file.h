// written by Szymon Winczewski

#ifndef central_forces_file_h
#define central_forces_file_h

#include <fstream>

#include "central_force.h"
#include "errors.h"


class CentralForcesFile
{
private:
    int mode_;     // 0 - not set, 1 - write, 2 - read, 3 - closed
    std::string file_name_;
    std::ifstream ifile_;
    std::ofstream ofile_;

    bool memory_allocated_;
    int number_of_atoms_;
    int max_number_of_central_forces_;
    int *number_of_central_forces_;
    central_force **central_forces_;


    void allocateMemory();
    void deallocateMemory();

    void checkAvailability();


public:
    CentralForcesFile();
    ~CentralForcesFile();

    void openOutFile(std::string file_name);
    void closeOutFile();
    void writeToOutFile(int number_of_atoms,
                        int max_number_of_central_forces,
                        int *number_of_central_forces,
                        central_force **central_forces);

    void openInFile(std::string file_name);
    void closeInFile();
    void readFromInFile();

    int getNumberOfAtoms();
    int getMaxNumberOfCentralForces();
    int *getNumberOfCentralForces();
    central_force **getCentralForces();
};

#endif
