// written by Szymon Winczewski

#include "cfd_smeam.h"

using namespace std;


CFD_sMEAM::CFD_sMEAM(int max_number_of_central_forces)
{
    #ifdef CFD_sMEAM_DEBUG
    std::cout << "*** CFD_sMEAM::CFD_sMEAM() called! ***" << std::endl;
    #endif

    max_number_of_central_forces_ = max_number_of_central_forces;

    #ifdef CFD_SMEAM_TIMINGS
    timer_ = new Timer("CFD_sMEAM", 20);
    timer_->addRoutine("readPotentialFile",        1);
    timer_->addRoutine("readNeighboursList",       2);
    timer_->addRoutine("compute_central_forces",   3);
    timer_->addRoutine("write_central_foces",      4);
    timer_->addRoutine("allocateMemory",           5);
    timer_->addRoutine("deallocateMemory",         6);
    #endif

    memory_allocated_ = 0;

    number_of_atoms_ = 0;

    max_number_of_n_neighbours_ = 0;
    n_num_ = NULL;
    n_list_ = NULL;
    n_bonds_ = NULL;

    max_number_of_s_neighbours_ = 0;
    s_num_ = NULL;
    s_list_ = NULL;
    s_bonds_ = NULL;

    number_of_central_forces_ = NULL;
    central_forces_ = NULL;

    #ifdef CFD_sMEAM_DEBUG
    std::cout << "*** CFD_sMEAM::CFD_sMEAM() done! ***" << std::endl;
    #endif
}


CFD_sMEAM::~CFD_sMEAM()
{
    #ifdef CFD_sMEAM_DEBUG
    std::cout << "*** CFD_sMEAM::~CFD_sMEAM() called! ***" << std::endl;
    #endif

    deallocateMemory();

    if ( phi_spline_ != NULL )
        delete phi_spline_;
    if ( rho_spline_ != NULL )
        delete rho_spline_;
    if ( U_spline_ != NULL )
        delete U_spline_;
    if ( f_spline_ != NULL )
        delete f_spline_;
    if ( g_spline_ != NULL )
        delete g_spline_;

    timer_->printTimings();

    delete timer_;

    #ifdef CFD_sMEAM_DEBUG
    std::cout << "*** CFD_sMEAM::~CFD_sMEAM() done! ***" << std::endl;
    #endif
}


void CFD_sMEAM::allocateMemory()
{
    #ifdef CFD_SMEAM_TIMINGS
    timer_->routineStarted(5);
    #endif

    if ( memory_allocated_ == 0 )
    {
        allocateCentralForces(number_of_atoms_, max_number_of_central_forces_,
                              number_of_central_forces_, central_forces_);

        n_num_ = new int [number_of_atoms_];
        n_list_ = new int *[number_of_atoms_];
        n_bonds_ = new vec3d *[number_of_atoms_];
        for (int i = 0; i < number_of_atoms_; i++)
        {
            n_list_[i] = new int [max_number_of_n_neighbours_];
            n_bonds_[i] = new vec3d [max_number_of_n_neighbours_];
        }

        s_num_ = new int [number_of_atoms_];
        s_list_ = new int *[number_of_atoms_];
        s_bonds_ = new vec3d *[number_of_atoms_];
        for (int i = 0; i < number_of_atoms_; i++)
        {
            s_list_[i] = new int [max_number_of_s_neighbours_];
            s_bonds_[i] = new vec3d [max_number_of_s_neighbours_];
        }

        Uprime_ = new double [number_of_atoms_];
        bonds_list_ = new meam_bond [max_number_of_n_neighbours_];

        number_of_bonds_central_ = new int [number_of_atoms_];
        bonds_list_central_ = new meam_bond *[number_of_atoms_];
        for (int i = 0; i < number_of_atoms_; i++)
            bonds_list_central_[i] = new meam_bond[max_number_of_n_neighbours_];
        smart_neighbours_list_ = new int *[number_of_atoms_];
        for (int i = 0; i < number_of_atoms_; i++)
            smart_neighbours_list_[i] = new int [2];

        #ifdef CFD_SMEAM_DEBUG
        std::cout << "memory usage:" << std::endl;
        std::cout << "  number_of_atoms_             = " << number_of_atoms_ << std::endl;
        std::cout << "  max_number_of_n_neighbours_  = " << max_number_of_n_neighbours_ << std::endl;
        std::cout << "  number_of_bonds_central_     = " << number_of_atoms_ * sizeof(int) << " Bytes" << std::endl;
        std::cout << "  bonds_list_central_          = " << number_of_atoms_ * max_number_of_n_neighbours_ * sizeof(meam_bond) << " Bytes" << std::endl;
        std::cout << "  smart_neighbours_list_       = " << number_of_atoms_ * 2 * sizeof(int) << " Bytes" << std::endl;
        #endif

        memory_allocated_ = 1;
    }

    #ifdef CFD_SMEAM_TIMINGS
    timer_->routineStopped(5);
    #endif
}


void CFD_sMEAM::deallocateMemory()
{
    #ifdef CFD_SMEAM_TIMINGS
    timer_->routineStarted(6);
    #endif

    if ( memory_allocated_ == 1 )
    {
        deallocateCentralForces(number_of_atoms_,
                                number_of_central_forces_, central_forces_);

        delete [] n_num_;
        for (int i = 0; i < number_of_atoms_; i++)
        {
            delete [] n_list_[i];
            delete [] n_bonds_[i];
        }
        delete [] n_list_;
        delete [] n_bonds_;

        delete [] s_num_;
        for (int i = 0; i < number_of_atoms_; i++)
        {
            delete [] s_list_[i];
            delete [] s_bonds_[i];
        }
        delete [] s_list_;
        delete [] s_bonds_;

        delete [] Uprime_;
        delete [] bonds_list_;

        delete [] number_of_bonds_central_;
        for (int i = 0; i < number_of_atoms_; i++)
            delete [] bonds_list_central_[i];
        delete [] bonds_list_central_;
        for (int i = 0; i < number_of_atoms_; i++)
            delete [] smart_neighbours_list_[i];
        delete [] smart_neighbours_list_;

        memory_allocated_ = 0;
    }

    #ifdef CFD_SMEAM_TIMINGS
    timer_->routineStopped(6);
    #endif
}


void CFD_sMEAM::readPotentialFile(std::string file_name)
{
    #ifdef CFD_SMEAM_TIMINGS
    timer_->routineStarted(1);
    #endif

    std::ifstream input;
    bool end_keyword_found;
    std::string keyword;

    input.open(file_name.c_str());
    openInputFileError(ERR_CFD_SMEAM, 1, "readPotentialFile", file_name, input);

    end_keyword_found = 0;
    input >> keyword;
    do
    {
        if ( keyword == "func_phi" )
        {
            #ifdef CFD_SMEAM_DEBUG
            std::cout << "reading phi..." << std::endl;
            #endif
            phi_spline_ = new Spline();
            phi_spline_->readFromASCIIFile(input);
            #ifdef CFD_SMEAM_DEBUG
            phi_spline_->print();
            std::cout << "done!" << std::endl;
            std::cout << std::endl;
            #endif
        }
        else if ( keyword == "func_rho" )
        {
            #ifdef CFD_SMEAM_DEBUG
            std::cout << "reading rho..." << std::endl;
            #endif
            rho_spline_ = new Spline();
            rho_spline_->readFromASCIIFile(input);
            #ifdef CFD_SMEAM_DEBUG
            rho_spline_->print();
            std::cout << "done!" << std::endl;
            std::cout << std::endl;
            #endif
        }
        else if ( keyword == "func_U" )
        {
            #ifdef CFD_SMEAM_DEBUG
            std::cout << "reading U..." << std::endl;
            #endif
            U_spline_ = new Spline();
            U_spline_->readFromASCIIFile(input);
            #ifdef CFD_SMEAM_DEBUG
            U_spline_->print();
            std::cout << "done!" << std::endl;
            std::cout << std::endl;
            #endif
        }
        else if ( keyword == "func_f" )
        {
            #ifdef CFD_SMEAM_DEBUG
            std::cout << "reading f..." << std::endl;
            #endif
            f_spline_ = new Spline();
            f_spline_->readFromASCIIFile(input);
            #ifdef CFD_SMEAM_DEBUG
            f_spline_->print();
            std::cout << "done!" << std::endl;
            std::cout << std::endl;
            #endif
        }
        else if ( keyword == "func_g" )
        {
            #ifdef CFD_SMEAM_DEBUG
            std::cout << "reading g..." << std::endl;
            #endif
            g_spline_ = new Spline();
            g_spline_->readFromASCIIFile(input);
            #ifdef CFD_SMEAM_DEBUG
            g_spline_->print();
            std::cout << "done!" << std::endl;
            std::cout << std::endl;
            #endif
        }
        else if ( keyword == "end" )
            end_keyword_found = 1;
        else
            std::cout << "unrecognized keyword \"" << keyword << "\" during parsing the config file!" << std::endl;
        input >> keyword;
    }
    while ( ( end_keyword_found == 0 ) && ( input.eof() == 0 ) );

    input.close();

    zero_atom_energy_ = U_spline_->evaluate(0.0);
    #ifdef CFD_SMEAM_DEBUG
    std::cout << "zero atom energy = " << zero_atom_energy_ << std::endl;
    #endif

    cutoff_ = 0.0;
    if ( phi_spline_->getCutoff() > cutoff_ )
        cutoff_ = phi_spline_->getCutoff();
    if ( rho_spline_->getCutoff() > cutoff_ )
        cutoff_ = rho_spline_->getCutoff();
    if ( f_spline_->getCutoff() > cutoff_ )
        cutoff_ = f_spline_->getCutoff();
    if ( cutoff_ <= 0.0 )
        raiseError(ERR_CFD_SMEAM, 2, "readPotentialFile", "incorrect cutoff radius");

    #ifdef CFD_SMEAM_DEBUG
    std::cout << "cutoffs: " << std::endl;
    std::cout << "   phi      " << phi_spline_->getCutoff() << std::endl;
    std::cout << "   rho      " << rho_spline_->getCutoff() << std::endl;
    std::cout << "   f        " << f_spline_->getCutoff() << std::endl;
    std::cout << "   global   " << cutoff_ << std::endl;
    #endif

    #ifdef CFD_SMEAM_TIMINGS
    timer_->routineStopped(1);
    #endif
}


void CFD_sMEAM::readNeighboursList(std::string file_name)
{
    #ifdef CFD_SMEAM_TIMINGS
    timer_->routineStarted(2);
    #endif

    std::ifstream in;
    int i, j;

    in.open(file_name.c_str(), std::ios::binary | std::ios::in);
    openInputFileError(ERR_CFD_SMEAM, 1, "readNeighboursList", file_name, in);

    in.read((char*)&number_of_atoms_, sizeof(int));
    in.read((char*)&max_number_of_n_neighbours_, sizeof(int));
    in.read((char*)&max_number_of_s_neighbours_, sizeof(int));

    allocateMemory();

    for (i = 0; i < number_of_atoms_; i++)
    {
        in.read((char*)&i, sizeof(int));

        in.read((char*)&n_num_[i], sizeof(int));
        for (j = 0; j < n_num_[i]; j++)
        {
            in.read((char*)&n_list_[i][j], sizeof(int));
            in.read((char*)&n_bonds_[i][j].vec[0], sizeof(double));
            in.read((char*)&n_bonds_[i][j].vec[1], sizeof(double));
            in.read((char*)&n_bonds_[i][j].vec[2], sizeof(double));
            in.read((char*)&n_bonds_[i][j].r, sizeof(double));
        }

        in.read((char*)&s_num_[i], sizeof(int));
        for (j = 0; j < s_num_[i]; j++)
        {
            in.read((char*)&s_list_[i][j], sizeof(int));
            in.read((char*)&s_bonds_[i][j].vec[0], sizeof(double));
            in.read((char*)&s_bonds_[i][j].vec[1], sizeof(double));
            in.read((char*)&s_bonds_[i][j].vec[2], sizeof(double));
            in.read((char*)&s_bonds_[i][j].r, sizeof(double));
        }
    }

    in.close();

    #ifdef CFD_SMEAM_TIMINGS
    timer_->routineStopped(2);
    #endif
}
