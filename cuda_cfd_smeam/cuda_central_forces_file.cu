// written by Szymon Winczewski

#include "cuda_central_forces_file.h"

using namespace std;


CudaCentralForcesFile::CudaCentralForcesFile()
{
    mode_ = 0;

    memory_allocated_ = 0;
    number_of_atoms_ = 0;
    max_number_of_central_forces_ = 0;
    number_of_central_forces_ = NULL;
    central_forces_ = NULL;
}


CudaCentralForcesFile::~CudaCentralForcesFile()
{
    deallocateMemory();
}


void CudaCentralForcesFile::checkAvailability()
{
    if ( memory_allocated_ == 0 )
        raiseError(ERR_CENTRAL_FORCES_FILE, 1, "checkAvailability", "the requested operation is forbidden");
}


void CudaCentralForcesFile::allocateMemory()
{
    if ( memory_allocated_ == 0 )
    {
        allocateCentralForces(number_of_atoms_,
                              max_number_of_central_forces_,
                              number_of_central_forces_,
                              central_forces_);
        memory_allocated_ = 1;
    }
}


void CudaCentralForcesFile::deallocateMemory()
{
    if ( memory_allocated_ == 1 )
    {
        deallocateCentralForces(number_of_atoms_,
                                number_of_central_forces_,
                                central_forces_);
        memory_allocated_ = 0;
    }
}


void CudaCentralForcesFile::openOutFile(std::string file_name)
{
    if ( mode_ != 0 )
        raiseError(ERR_CENTRAL_FORCES_FILE, 2, "openOutFile", "the requested operation is forbidden");
    file_name_ = file_name;
    ofile_.open(file_name_.c_str(), std::ios::binary | std::ios::out);
    openOutputFileError(ERR_CENTRAL_FORCES_FILE, 3, "openOutFile", file_name_, ofile_);
    mode_ = 1;
}


void CudaCentralForcesFile::closeOutFile()
{
    if ( mode_ != 1 )
        raiseError(ERR_CENTRAL_FORCES_FILE, 4, "closeOutFile", "the requested operation is forbidden");
    ofile_.close();
    ofile_.clear();
    mode_ = 3;
}


void CudaCentralForcesFile::writeToOutFile(int number_of_atoms,
                                       int max_number_of_central_forces,
                                       int *number_of_central_forces,
                                       central_force **central_forces)
{
    if ( mode_ != 1 )
        raiseError(ERR_CENTRAL_FORCES_FILE, 5, "writeToOutFile", "the requested operation is forbidden");

    int i, j;
    central_force *this_central_force;

    ofile_.write((char*)&number_of_atoms, sizeof(int));
    ofile_.write((char*)&max_number_of_central_forces, sizeof(int));
    for (i = 0; i < number_of_atoms; i++)
    {
        ofile_.write((char*)&number_of_central_forces[i], sizeof(int));
        for (j = 0; j < number_of_central_forces[i]; j++)
        {
            this_central_force = &central_forces[i][j];
            ofile_.write((char*)&this_central_force->first_second, sizeof(bool));
            ofile_.write((char*)&this_central_force->atom_j, sizeof(int));
            ofile_.write((char*)&this_central_force->force[0], sizeof(double));
            ofile_.write((char*)&this_central_force->force[1], sizeof(double));
            ofile_.write((char*)&this_central_force->force[2], sizeof(double));
            ofile_.write((char*)&this_central_force->r_ij, sizeof(double));
            ofile_.write((char*)&this_central_force->r_ij_dir[0], sizeof(double));
            ofile_.write((char*)&this_central_force->r_ij_dir[1], sizeof(double));
            ofile_.write((char*)&this_central_force->r_ij_dir[2], sizeof(double));
        }
    }
}


void CudaCentralForcesFile::openInFile(std::string file_name)
{
    if ( mode_ != 0 )
        raiseError(ERR_CENTRAL_FORCES_FILE, 6, "openInFile", "the requested operation is forbidden");
    file_name_ = file_name;
    ifile_.open(file_name_.c_str(), std::ios::binary | std::ios::in);
    openInputFileError(ERR_CENTRAL_FORCES_FILE, 7, "openInFile", file_name_, ifile_);
    mode_ = 2;
}


void CudaCentralForcesFile::closeInFile()
{
    if ( mode_ != 2 )
        raiseError(ERR_CENTRAL_FORCES_FILE, 8, "closeInFile", "the requested operation is forbidden");
    ifile_.close();
    ifile_.clear();
    mode_ = 3;
}


void CudaCentralForcesFile::readFromInFile()
{
    if ( mode_ != 2 )
        raiseError(ERR_CENTRAL_FORCES_FILE, 9, "readFromInFile", "the requested operation is forbidden");

    int i, j;
    central_force *this_central_force;
    int tmp_number_of_atoms;
    int tmp_max_number_of_central_forces;

    ifile_.read((char*)&tmp_number_of_atoms, sizeof(int));
    if ( ifile_.eof() == 1 )
        raiseError(ERR_CENTRAL_FORCES_FILE, 10, "readFromInFile", "end of file reached");

    ifile_.read((char*)&tmp_max_number_of_central_forces, sizeof(int));
    if ( ifile_.eof() == 1 )
        raiseError(ERR_CENTRAL_FORCES_FILE, 10, "readFromInFile", "end of file reached");

    if ( memory_allocated_ == 0 )
    {
        number_of_atoms_ = tmp_number_of_atoms;
        max_number_of_central_forces_ = tmp_max_number_of_central_forces;
        allocateMemory();
    }
    else
    {
        if ( tmp_number_of_atoms != number_of_atoms_ )
            raiseError(ERR_CENTRAL_FORCES_FILE, 11, "readFromInFile", "number of atoms has changed");
        if ( tmp_max_number_of_central_forces != max_number_of_central_forces_ )
            raiseError(ERR_CENTRAL_FORCES_FILE, 12, "readFromInFile", "max number of central forces has changed");
        number_of_atoms_ = tmp_number_of_atoms;
        max_number_of_central_forces_ = tmp_max_number_of_central_forces;
    }

    for (i = 0; i < number_of_atoms_; i++)
    {
        ifile_.read((char*)&number_of_central_forces_[i], sizeof(int));
        if ( ifile_.eof() == 1 )
            raiseError(ERR_CENTRAL_FORCES_FILE, 10, "readFromInFile", "end of file reached");

        for (j = 0; j < number_of_central_forces_[i]; j++)
        {
            this_central_force = &central_forces_[i][j];
            ifile_.read((char*)&this_central_force->first_second, sizeof(bool));
            ifile_.read((char*)&this_central_force->atom_j, sizeof(int));
            ifile_.read((char*)&this_central_force->force[0], sizeof(double));
            ifile_.read((char*)&this_central_force->force[1], sizeof(double));
            ifile_.read((char*)&this_central_force->force[2], sizeof(double));
            ifile_.read((char*)&this_central_force->r_ij, sizeof(double));
            ifile_.read((char*)&this_central_force->r_ij_dir[0], sizeof(double));
            ifile_.read((char*)&this_central_force->r_ij_dir[1], sizeof(double));
            ifile_.read((char*)&this_central_force->r_ij_dir[2], sizeof(double));
        }
    }
}


int CudaCentralForcesFile::getNumberOfAtoms()
{
    checkAvailability();
    return number_of_atoms_;
}


int CudaCentralForcesFile::getMaxNumberOfCentralForces()
{
    checkAvailability();
    return max_number_of_central_forces_;
}


int *CudaCentralForcesFile::getNumberOfCentralForces()
{
    checkAvailability();
    return number_of_central_forces_;
}


central_force **CudaCentralForcesFile::getCentralForces()
{
    checkAvailability();
    return central_forces_;
}
