// written by Szymon Winczewski

#include <iostream>

#include "errors.h"
#include "cfd_smeam.h"
#include "cuda_cfd_smeam.h"

using namespace std;


int main(int argc, char *argv[])
{
    initializeErrors();

	CUDA_CFD_sMEAM *my_cuda_cfd_for_sMEAM = new CUDA_CFD_sMEAM(200);
	my_cuda_cfd_for_sMEAM->readNeighboursList("bin/nbh_list.bin");
	my_cuda_cfd_for_sMEAM->readPotentialFile("bin/Mo.MEAM.params");
	my_cuda_cfd_for_sMEAM->readNeighboursList("bin/nbh_list.bin");
	my_cuda_cfd_for_sMEAM->compute_central_forces();
	my_cuda_cfd_for_sMEAM->write_central_foces("bin/cforces.bin");
    /*CFD_sMEAM *my_cfd_for_sMEAM;
    my_cfd_for_sMEAM = new CFD_sMEAM(200);
    my_cfd_for_sMEAM->readPotentialFile("bin/Mo.MEAM.params");
    my_cfd_for_sMEAM->readNeighboursList("bin/nbh_list.bin");
    my_cfd_for_sMEAM->compute_central_forces();
    my_cfd_for_sMEAM->write_central_foces("bin/cforces.bin");*/

    /*delete my_cfd_for_sMEAM;*/
	system("pause");
    destroyErrors();
	
    return 0;
}
