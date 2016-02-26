#include "cudaerror.h"

//w przypadku jakiegokolwiek b��du w operacjach cuda wyrzuca w konsoli b��d i ko�czy dzia�anie programu
void RaiseError(cudaError cudaStatus)
{
	if (cudaStatus != cudaSuccess)
	{
		printf("error, %s", cudaGetErrorString(cudaStatus));
		exit(-1);
	}
}