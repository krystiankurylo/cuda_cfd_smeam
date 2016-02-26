#include "cudaerror.h"

//w przypadku jakiegokolwiek b³êdu w operacjach cuda wyrzuca w konsoli b³¹d i koñczy dzia³anie programu
void RaiseError(cudaError cudaStatus)
{
	if (cudaStatus != cudaSuccess)
	{
		printf("error, %s", cudaGetErrorString(cudaStatus));
		exit(-1);
	}
}