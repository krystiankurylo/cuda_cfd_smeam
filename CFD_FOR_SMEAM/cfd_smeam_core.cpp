// written by Szymon Winczewski

#include "cfd_smeam.h"

using namespace std;


int CFD_sMEAM::createBondsList(int atom_i, meam_bond *bonds_list)
{
// tworzymy liste wiazan atomu i-tego,
// dla kazdego z sasiadow na liscie wiazan wyznaczamy i zapamietujemy:
// id sasiada, r_ij, r_ij_dir, f_r_ij, fprime_r_ij
    int jj;
    vec3d *this_bond;
    double r_ij, inv_of_r_ij;

    int number_of_bonds = n_num_[atom_i];
    for (jj = 0; jj  < number_of_bonds; jj++)
    {
        this_bond = &n_bonds_[atom_i][jj];
        bonds_list[jj].atom_j = n_list_[atom_i][jj];
        r_ij = this_bond->r;
        bonds_list[jj].r_ij = r_ij;
        inv_of_r_ij = 1.0 / r_ij;
        bonds_list[jj].inv_of_r_ij = inv_of_r_ij;
        bonds_list[jj].r_ij_dir[0] = inv_of_r_ij * this_bond->vec[0];
        bonds_list[jj].r_ij_dir[1] = inv_of_r_ij * this_bond->vec[1];
        bonds_list[jj].r_ij_dir[2] = inv_of_r_ij * this_bond->vec[2];
        bonds_list[jj].f_r_ij = f_spline_->evaluate(r_ij, bonds_list[jj].fprime_r_ij);
    }

    return number_of_bonds;
}


double CFD_sMEAM::computeElectronDensity(int number_of_bonds, meam_bond *bonds_list)
{
// wyznaczamy gestosc n_i dana jako:
//    n_i = sum_j rho(r_ij) + 1/2 sum_jk f(r_ij) f(r_ik) g(cos_theta_jik)
    int j, k;
    double n_i_total = 0.0;
    double n_i_three_body;
    double cos_theta_jik;
    meam_bond bond_ij;
    meam_bond bond_ik;

    for (j = 0; j < number_of_bonds; j++)
    {
        memcpy(&bond_ij, &bonds_list[j], sizeof(meam_bond));

        n_i_three_body = 0.0;
        for (k = j + 1; k < number_of_bonds; k++)
        {
            memcpy(&bond_ik, &bonds_list[k], sizeof(meam_bond));

            computeCosThetaJIK(bond_ij.r_ij_dir, bond_ik.r_ij_dir, cos_theta_jik);

            n_i_three_body += bond_ik.f_r_ij * g_spline_->evaluate(cos_theta_jik);
        }

        n_i_total += bond_ij.f_r_ij * n_i_three_body;
        n_i_total += rho_spline_->evaluate(bond_ij.r_ij);
    }

    return n_i_total;
}


void CFD_sMEAM::compute_central_forces()
{
    #ifdef CFD_SMEAM_TIMINGS
    timer_->routineStarted(3);
    #endif

    #ifdef CFD_SMEAM_DEBUG
    std::cout << "*** CFD_sMEAM::compute_central_forces() called! ***" << std::endl;
    #endif

    int i, j;
    int atom_i, atom_j, atom_k;
    int number_of_bonds_i, number_of_bonds_j;
    meam_bond *bonds_list_i, *bonds_list_j;
    double n_i_total;

    double Theta1_jik, Theta2_jik, Theta3_jik;
    double Theta1_ijk, Theta2_ijk, Theta3_ijk;

    double Theta4_kji;
    double r_ij_vec[3];

    bool nonzero_force;

    vec3d *bond_r_ij;
    double r_ij, inv_of_r_ij, r_ij_dir[3];
    double f_r_ij, fprime_r_ij;
    double inv_of_r_ik, r_ik_dir[3];
    double f_r_ik;

    double cos_theta_jik, g_cos_theta_jik, gprime_cos_theta_jik;
    double inv_of_r_jk, r_jk_dir[3];
    double f_r_jk;
    double cos_theta_kji, g_cos_theta_kji, gprime_cos_theta_kji;
    double prefactor;
    double phi_prime_r_ij, rho_prime_r_ij;

    central_force *this_central_force;
    int tmp_number_of_central_forces;
    double tmp_central_force[3];

// sprytna lista sasiadow - zerowanie listy
    int minus_one = -1;
    for (atom_i = 0; atom_i < number_of_atoms_; atom_i++)
        memcpy(&smart_neighbours_list_[atom_i][0], &minus_one, sizeof(int));

	
// dla kazdego z atomow tworzymy liste wiazan oraz wyznaczamy n_i, U(n_i), U'(n_i)
    for (atom_i = 0; atom_i < number_of_atoms_; atom_i++)
    {
        number_of_bonds_central_[atom_i] = createBondsList(atom_i, bonds_list_central_[atom_i]);
        n_i_total = computeElectronDensity(number_of_bonds_central_[atom_i], bonds_list_central_[atom_i]);
        Uprime_[atom_i] = U_spline_->evaluate_deriv(n_i_total);
    }

    int max_number_of_central_forces = 0;

    for (atom_i = 0; atom_i < number_of_atoms_; atom_i++)
    {
        tmp_number_of_central_forces = 0;

        number_of_bonds_i = number_of_bonds_central_[atom_i];
        bonds_list_i = bonds_list_central_[atom_i];

// sprytna lista sasiadow - aktualizacja
        for (j = 0; j < n_num_[atom_i]; j++)
        {
            atom_j = n_list_[atom_i][j];
            memcpy(&smart_neighbours_list_[atom_j][0], &atom_i, sizeof(int));
            memcpy(&smart_neighbours_list_[atom_j][1], &j, sizeof(int));
        }

// start: sily centralne typu NN
        for (i = 0; i < number_of_bonds_i; i++)
        {
            f_r_ij = bonds_list_i[i].f_r_ij;

            if ( f_r_ij != 0.0 )
            {
                atom_j = bonds_list_i[i].atom_j;
                r_ij = bonds_list_i[i].r_ij;
                inv_of_r_ij = bonds_list_i[i].inv_of_r_ij;
                r_ij_dir[0] = bonds_list_i[i].r_ij_dir[0];
                r_ij_dir[1] = bonds_list_i[i].r_ij_dir[1];
                r_ij_dir[2] = bonds_list_i[i].r_ij_dir[2];
                fprime_r_ij = bonds_list_i[i].fprime_r_ij;

// obliczanie Theta1_jik, Theta2_jik, Theta3_jik
                Theta1_jik = 0.0;
                Theta2_jik = 0.0;
                Theta3_jik = 0.0;
                for (j = 0; j < number_of_bonds_i; j++)
                    if ( j != i )
                    {
                        f_r_ik = bonds_list_i[j].f_r_ij;

                        if ( f_r_ik == 0.0 )
                            continue;

                        inv_of_r_ik = bonds_list_i[j].inv_of_r_ij;
                        r_ik_dir[0] = bonds_list_i[j].r_ij_dir[0];
                        r_ik_dir[1] = bonds_list_i[j].r_ij_dir[1];
                        r_ik_dir[2] = bonds_list_i[j].r_ij_dir[2];

                        computeCosThetaJIK(r_ij_dir, r_ik_dir, cos_theta_jik);
                        g_cos_theta_jik = g_spline_->evaluate(cos_theta_jik, gprime_cos_theta_jik);

                        Theta1_jik += ( f_r_ik * g_cos_theta_jik );
                        Theta2_jik += ( inv_of_r_ik * f_r_ik * gprime_cos_theta_jik );
                        Theta3_jik += ( f_r_ik * gprime_cos_theta_jik * cos_theta_jik );
                    }

// obliczanie Theta1_ijk, Theta2_ijk, Theta3_ijk
                Theta1_ijk = 0.0;
                Theta2_ijk = 0.0;
                Theta3_ijk = 0.0;

                number_of_bonds_j = number_of_bonds_central_[atom_j];
                bonds_list_j = bonds_list_central_[atom_j];
                for (j = 0; j < number_of_bonds_j; j++)
                {
                    atom_k = bonds_list_j[j].atom_j;
                    if ( atom_k == atom_i )
                        continue;

                    f_r_jk = bonds_list_j[j].f_r_ij;
                    if ( f_r_jk == 0.0 )
                        continue;

                    inv_of_r_jk = bonds_list_j[j].inv_of_r_ij;
                    r_jk_dir[0] = bonds_list_j[j].r_ij_dir[0];
                    r_jk_dir[1] = bonds_list_j[j].r_ij_dir[1];
                    r_jk_dir[2] = bonds_list_j[j].r_ij_dir[2];

                    computeCosThetaJIK(r_jk_dir, r_ij_dir, cos_theta_kji);
                    cos_theta_kji *= -1.0;
                    g_cos_theta_kji = g_spline_->evaluate(cos_theta_kji, gprime_cos_theta_kji);

                    Theta1_ijk += ( f_r_jk * g_cos_theta_kji );
                    Theta2_ijk += ( inv_of_r_jk * f_r_jk * gprime_cos_theta_kji );
                    Theta3_ijk += ( f_r_jk * gprime_cos_theta_kji * cos_theta_kji );
                }

                prefactor  = Uprime_[atom_i] * ( fprime_r_ij * Theta1_jik + f_r_ij * Theta2_jik - f_r_ij * inv_of_r_ij * Theta3_jik );
                prefactor += Uprime_[atom_j] * ( fprime_r_ij * Theta1_ijk + f_r_ij * Theta2_ijk - f_r_ij * inv_of_r_ij * Theta3_ijk );
            }
            else
            {
                prefactor = 0.0;
                atom_j = bonds_list_i[i].atom_j;
                r_ij = bonds_list_i[i].r_ij;
                r_ij_dir[0] = bonds_list_i[i].r_ij_dir[0];
                r_ij_dir[1] = bonds_list_i[i].r_ij_dir[1];
                r_ij_dir[2] = bonds_list_i[i].r_ij_dir[2];
                number_of_bonds_j = number_of_bonds_central_[atom_j];
                bonds_list_j = bonds_list_central_[atom_j];
            }

            phi_prime_r_ij = phi_spline_->evaluate_deriv(r_ij);
            prefactor += phi_prime_r_ij;

            rho_prime_r_ij = rho_spline_->evaluate_deriv(r_ij);
            prefactor += ( rho_prime_r_ij * ( Uprime_[atom_i] + Uprime_[atom_j] ) );

            tmp_central_force[0] = prefactor * r_ij_dir[0];
            tmp_central_force[1] = prefactor * r_ij_dir[1];
            tmp_central_force[2] = prefactor * r_ij_dir[2];

            compute_Theta4_kji(atom_i, number_of_bonds_j,
                               bonds_list_i, bonds_list_j,
                               Theta4_kji);

            tmp_central_force[0] -= Theta4_kji * r_ij_dir[0] * r_ij;
            tmp_central_force[1] -= Theta4_kji * r_ij_dir[1] * r_ij;
            tmp_central_force[2] -= Theta4_kji * r_ij_dir[2] * r_ij;

            if ( tmp_number_of_central_forces == max_number_of_central_forces_ )
                raiseError(ERR_CFD_SMEAM, 2, "compute_central_forces", "too many central forces");

            this_central_force = &central_forces_[atom_i][tmp_number_of_central_forces];
            this_central_force->first_second = 0;
            this_central_force->atom_j = atom_j;
            this_central_force->force[0] = tmp_central_force[0];
            this_central_force->force[1] = tmp_central_force[1];
            this_central_force->force[2] = tmp_central_force[2];
            this_central_force->r_ij = r_ij;
            this_central_force->r_ij_dir[0] = r_ij_dir[0];
            this_central_force->r_ij_dir[1] = r_ij_dir[1];
            this_central_force->r_ij_dir[2] = r_ij_dir[2];
            tmp_number_of_central_forces++;
        }
// stop: sily centralne typu NN

// uwaga: poprzez czlon Theta4_kji z atomem i-tym oddzialywuja rowniez pozostale atomy ukladu,
//        w powyzszej petli obliczylismy juz sily dla atomow bedacych najblizszymi sasiadami atomu i-tego,
//        wymagane jest jeszcze obliczenie pozostalych sil ,,trojcialowych''
//        pomiedzy atomem i-tym a atomami, ktore nie sa jego bezposrednimi sasiadami,
//        atomy takie charakteryzuja sie tym, iz posiadaja razem z atomem i-tym co najmniej jednego wspolnego sasiada

// start: sily centralne typu non-NN
        for (i = 0; i < s_num_[atom_i]; i++)
        {
            atom_j = s_list_[atom_i][i];
            number_of_bonds_j = number_of_bonds_central_[atom_j];
            bonds_list_j = bonds_list_central_[atom_j];
            nonzero_force = compute_Theta4_kji(atom_i, number_of_bonds_j,
                                               bonds_list_i, bonds_list_j,
                                               Theta4_kji);

            if ( nonzero_force == 1 )
            {
                bond_r_ij = &s_bonds_[atom_i][i];
                r_ij = bond_r_ij->r;
                inv_of_r_ij = 1.0 / r_ij;
                r_ij_dir[0] = bond_r_ij->vec[0] * inv_of_r_ij;
                r_ij_dir[1] = bond_r_ij->vec[1] * inv_of_r_ij;
                r_ij_dir[2] = bond_r_ij->vec[2] * inv_of_r_ij;

                r_ij_vec[0] = bond_r_ij->vec[0];
                r_ij_vec[1] = bond_r_ij->vec[1];
                r_ij_vec[2] = bond_r_ij->vec[2];
                tmp_central_force[0] = - Theta4_kji * r_ij_vec[0];
                tmp_central_force[1] = - Theta4_kji * r_ij_vec[1];
                tmp_central_force[2] = - Theta4_kji * r_ij_vec[2];

                if ( tmp_number_of_central_forces == max_number_of_central_forces_ )
                    raiseError(ERR_CFD_SMEAM, 3, "compute_central_forces", "too many central forces");

                this_central_force = &central_forces_[atom_i][tmp_number_of_central_forces];
                this_central_force->first_second = 1;
                this_central_force->atom_j = atom_j;
                this_central_force->force[0] = tmp_central_force[0];
                this_central_force->force[1] = tmp_central_force[1];
                this_central_force->force[2] = tmp_central_force[2];
                this_central_force->r_ij = r_ij;
                this_central_force->r_ij_dir[0] = r_ij_dir[0];
                this_central_force->r_ij_dir[1] = r_ij_dir[1];
                this_central_force->r_ij_dir[2] = r_ij_dir[2];
                tmp_number_of_central_forces++;
            }
        }
// stop: sily centralne typu non-NN

        number_of_central_forces_[atom_i] = tmp_number_of_central_forces;
        if ( tmp_number_of_central_forces > max_number_of_central_forces )
            max_number_of_central_forces = tmp_number_of_central_forces;
    }

    #ifdef CFD_SMEAM_DEBUG
    std::cout << std::endl;
    std::cout << "*** CFD_sMEAM::compute_central() done! ***" << std::endl;
    std::cout << std::endl;
    #endif

    #ifdef CFD_SMEAM_TIMINGS
    timer_->routineStopped(3);
    #endif
}


void CFD_sMEAM::write_central_foces(std::string file_name)
{
    #ifdef CFD_SMEAM_TIMINGS
    timer_->routineStarted(4);
    #endif

    CentralForcesFile *cforces_file;

    cforces_file = new CentralForcesFile();
    cforces_file->openOutFile(file_name);
    cforces_file->writeToOutFile(number_of_atoms_, max_number_of_central_forces_,
                                 number_of_central_forces_, central_forces_);
    cforces_file->closeOutFile();
    delete cforces_file;

    #ifdef CFD_SMEAM_TIMINGS
    timer_->routineStopped(4);
    #endif
}
