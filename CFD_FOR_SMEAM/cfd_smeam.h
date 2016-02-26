// written by Szymon Winczewski

#ifndef cfd_smeam_h
#define cfd_smeam_h

#include "defines.h"

#include <cmath>
#include <sstream>
#include "cuda_cfd_smeam.h"
#include "central_force.h"
#include "central_forces_file.h"
#include "errors.h"
#include "spline.h"

#ifdef CFD_SMEAM_TIMINGS
#include "timer.h"
#endif



class CFD_sMEAM
{
private:
    #ifdef CFD_SMEAM_TIMINGS
    Timer *timer_;
    #endif

    bool memory_allocated_;

    Spline *phi_spline_, *rho_spline_, *U_spline_, *f_spline_, *g_spline_;
    double cutoff_;
    double zero_atom_energy_;

    int number_of_atoms_;

    int max_number_of_n_neighbours_;
    int *n_num_;
    int **n_list_;
    vec3d **n_bonds_;

    int max_number_of_s_neighbours_;
    int *s_num_;
    int **s_list_;
    vec3d **s_bonds_;

    double *Uprime_;
// pochodna energi wbudowania poszczegolnych atomow

    meam_bond *bonds_list_;

    int *number_of_bonds_central_;
    meam_bond **bonds_list_central_;

    int **smart_neighbours_list_;
// sprytna lista sasiadow, wektor o wymiarach number_of_atoms_ x 2
// jezeli w polu [i][0] wpisana jest wartosc atom_j
// to atom i-ty jest sasiadem atomu j-tego
// pole [i][1] informuje o tym, ktorym z kolei sasiadem atomu j-tego jest atom i-ty
// (tj. ktora z kolei pozycje zajmuje atom j-ty na liscie sasiadow i wiazan atomu i-tego)

    int max_number_of_central_forces_;
    int *number_of_central_forces_;
    central_force **central_forces_;

    void allocateMemory();
    void deallocateMemory();

    int createBondsList(int atom_i, meam_bond *bonds_list);
    double computeElectronDensity(int number_of_bonds, meam_bond *bonds_list);

    void computeCosThetaJIK(double r_ij_dir[3], double r_ik_dir[3], double &cos_theta_jik);
    bool compute_Theta4_kji(int atom_i, int number_of_bonds_j,
                            meam_bond *bonds_list_i, meam_bond *bonds_list_j,
                            double &T4kji);


public:
    CFD_sMEAM(int max_number_of_central_forces);
    ~CFD_sMEAM();

    void readPotentialFile(std::string file_name);
    void readNeighboursList(std::string file_name);

    void compute_central_forces();
    void write_central_foces(std::string file_name);
};


inline void CFD_sMEAM::computeCosThetaJIK(double r_ij_dir[3], double r_ik_dir[3], double &cos_theta_jik)
{
    cos_theta_jik = r_ij_dir[0] * r_ik_dir[0] + r_ij_dir[1] * r_ik_dir[1] + r_ij_dir[2] * r_ik_dir[2];
}


inline bool CFD_sMEAM::compute_Theta4_kji(int atom_i, int number_of_bonds_j,
                                          meam_bond *bonds_list_i, meam_bond *bonds_list_j,
                                          double &T4kji)
{
// nie ma koniecznosci by atom j-ty byl sasiadem atomu i-tego,
// musimy zlokalizowac te atomy k-te, ktore rownoczesnie sa sasiadami
// atomu i-tego oraz atomu j-tego
    bool nonzero_force = 0;
    int k, atom_k;
    double r_jk, inv_of_r_jk, r_jk_dir[3];
    double f_r_jk;
    double r_ik, inv_of_r_ik, r_ik_dir[3];
    double f_r_ik;
    double cos_theta_kji, gprime_cos_theta_kji;

    int jjj;

    T4kji = 0.0;

// petla po atomach j-tych bedacych sasiadami atomy i-tego
    for (k = 0; k < number_of_bonds_j; k++)
    {
        atom_k = bonds_list_j[k].atom_j;
        if ( atom_k == atom_i )
            continue;

        f_r_jk = bonds_list_j[k].f_r_ij;
        if ( f_r_jk == 0.0 )
            continue;

        r_jk = bonds_list_j[k].r_ij;
        inv_of_r_jk = 1.0 / r_jk;
        r_jk_dir[0] = bonds_list_j[k].r_ij_dir[0];
        r_jk_dir[1] = bonds_list_j[k].r_ij_dir[1];
        r_jk_dir[2] = bonds_list_j[k].r_ij_dir[2];

// sprawdzamy czy atom k-ty jest sasiadem atomu i-tego
// wykorzystujemy sprytna liste sasiadow
        if ( smart_neighbours_list_[atom_k][0] == atom_i )
        {
            jjj = smart_neighbours_list_[atom_k][1];

            f_r_ik = bonds_list_i[jjj].f_r_ij;
            if ( f_r_ik == 0.0 )
                continue;

            r_ik = bonds_list_i[jjj].r_ij;
            inv_of_r_ik = 1.0 / r_ik;
            r_ik_dir[0] = bonds_list_i[jjj].r_ij_dir[0];
            r_ik_dir[1] = bonds_list_i[jjj].r_ij_dir[1];
            r_ik_dir[2] = bonds_list_i[jjj].r_ij_dir[2];

            computeCosThetaJIK(r_jk_dir, r_ik_dir, cos_theta_kji);

            gprime_cos_theta_kji = g_spline_->evaluate_deriv(cos_theta_kji);

            T4kji += Uprime_[atom_k] * f_r_jk * inv_of_r_jk * f_r_ik * inv_of_r_ik * gprime_cos_theta_kji;

            nonzero_force = 1;
        }
    }

    return nonzero_force;
}

#endif
