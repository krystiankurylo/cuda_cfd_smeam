#ifndef cudaerrors_h
#define cudaerrors_h

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

#include "cuda_runtime.h"
#include "device_launch_parameters.h"


void RaiseError(cudaError);


#endif