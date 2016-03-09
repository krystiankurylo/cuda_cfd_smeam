#include <stdio.h>


#include "cudaspline.cuh"



#include <windows.h>
#include <time.h>
#include "cudaerror.h"

//kernel obliczania pochodnej
__global__ void kernelEvaluateDeriv(double* result,double* array_x,int size,double xmin_,double deriv0_,double derivN_,double h_,double hsq_,double xmax_shifted_,values_struct* values_vec_)
{
	int tidx = threadIdx.x + blockDim.x * blockIdx.x;

	if (tidx < size)
	{
		double x = array_x[tidx];

		x -= xmin_;
		if (x <= 0.0)
		{
			result[tidx] = deriv0_;
		}
		else if (x >= xmax_shifted_)
		{
			result[tidx] = derivN_;
		}
			
		else
		{
			int k = (int)(x / h_);
			double a = values_vec_[k].Xs_next - x;
			double b = h_ - a;
			double asq = a * a;
			double bsq = b * b;
			result[tidx] = values_vec_[k].Ydelta + ((3.0 * bsq - hsq_) * values_vec_[k].Y2_next - (3.0 * asq - hsq_) * values_vec_[k].Y2);
		}
	}
}


//obliczanie pochodnej
void CudaSpline::EvaluateDerivCUDA(double* result,double *array_x,int size)
{
	cudaError cudaStatus;

	dim3 block(1024, 1024);
	dim3 grid(32, 32);

	double* dev_result;
	double* dev_array_x;
	values_struct* dev_values_vec;

	//alokowanie pamiêci dla wyniku
	cudaStatus = cudaMalloc(&dev_result, size*sizeof(double));

	RaiseError(cudaStatus);

	//alokowanie pamiêci dla dev_array_x
	cudaStatus = cudaMalloc(&dev_array_x, size*sizeof(double));

	RaiseError(cudaStatus);

	//kopiowaine pamiêci dla dev_array_x
	cudaStatus = cudaMemcpy(dev_array_x, array_x, size*sizeof(double), cudaMemcpyHostToDevice);

	RaiseError(cudaStatus);

	//alokowanie pamiêci dla dev_values_vec

	cudaStatus = cudaMalloc(&dev_values_vec,N_*sizeof(values_struct));

	RaiseError(cudaStatus);

	cudaStatus = cudaMemcpy(dev_values_vec, values_vec_, N_*sizeof(values_struct), cudaMemcpyHostToDevice);

	RaiseError(cudaStatus);

	kernelEvaluateDeriv<<<block,grid>>>(dev_result, dev_array_x, size, xmin_, deriv0_, derivN_, h_, hsq_, xmax_shifted_, dev_values_vec);

	//kopiowanie pamiêci dla result
	cudaStatus = cudaMemcpy(result, dev_result, size*sizeof(double), cudaMemcpyDeviceToHost);

	RaiseError(cudaStatus);

	cudaFree(dev_values_vec);
	cudaFree(dev_array_x);
	cudaFree(dev_result);
}

CudaSpline::CudaSpline()
{
	X_ = NULL;
	Xs_ = NULL;
	Y_ = NULL;
	Y2_ = NULL;
	Ydelta_ = NULL;
	N_ = 0;
	value0_ = 0.0;
	deriv0_ = 0.0;
	valueN_ = 0.0;
	derivN_ = 0.0;
	xmin_ = 0.0;
	xmax_ = 0.0;
	cutoff_ = 0.0;
	h_ = 0.0;
	hsq_ = 0.0;
	xmax_shifted_ = 0.0;
#ifdef SPLINE_FAST
	values_vec_ = NULL;
#endif
}
CudaSpline::~CudaSpline()
{
	if (X_ != NULL)
		delete[] X_;
	if (Xs_ != NULL)
		delete[] Xs_;
	if (Y_ != NULL)
		delete[] Y_;
	if (Y2_ != NULL)
		delete[] Y2_;
	if (Ydelta_ != NULL)
		delete[] Ydelta_;
#ifdef SPLINE_FAST
	if (values_vec_ != NULL)
		delete[] values_vec_;
#endif
}

double CudaSpline::getCutoff()
{
	return cutoff_;
}

//stary evaluate - nieu¿ywany
__global__ void evaluateKernel(pair* results, int N, double x, double xmin_, double h_, double xmax_shifted_, double dx, double value0_, double deriv0_, double valueN_, double derivN_, double hsq_, values_struct* values_vec_)
{
	int tid = blockIdx.x;
	if (tid < N)
	{
		x += tid*dx;
		results[tid].key = x;
		x -= xmin_;
		if (x <= 0.0)
		{
			results[tid].value = value0_ + deriv0_ * x;
		}
		else if (x >= xmax_shifted_)
		{
			results[tid].value = valueN_ + derivN_ * (x - xmax_shifted_);
		}
		else
		{
			int k = (int)(x / h_);
			double a = values_vec_[k].Xs_next - x;
			double b = h_ - a;
			results[tid].value = values_vec_[k].Y_next - a * values_vec_[k].Ydelta + ((a * a - hsq_) * a * values_vec_[k].Y2 + (b * b - hsq_) * b * values_vec_[k].Y2_next);
		}
	}
    
}

//kernel evaluate
__global__ void evaluateKernel(vec3d* n_bonds_1D, meam_bond* bonds_list_central_1D,int* n_num_, int number_of_atoms_, int max_number_of_n_neighbours_,double xmin_,double deriv0_,double derivN_,double value0_,double valueN_,double xmax_shifted_,values_struct *values_vec_,double h_,double hsq_)
{

	int tidx = threadIdx.x + blockDim.x*blockIdx.x;
	int tidy = threadIdx.y + blockDim.y*blockIdx.y;



	if (tidx < number_of_atoms_)
	{	
			int number_of_bonds = n_num_[tidx];
			if (tidy < number_of_bonds)
			{
				vec3d *this_bond = &n_bonds_1D[tidx*max_number_of_n_neighbours_ + tidy];
				double r_ij = this_bond->r;
				r_ij -= xmin_;
				if (r_ij <= 0.0)
				{
					bonds_list_central_1D[tidx*max_number_of_n_neighbours_ + tidy].fprime_r_ij = deriv0_;
					bonds_list_central_1D[tidx*max_number_of_n_neighbours_ + tidy].f_r_ij = value0_ + deriv0_ * r_ij;
				}
				else if (r_ij >= xmax_shifted_)
				{
					bonds_list_central_1D[tidx*max_number_of_n_neighbours_ + tidy].fprime_r_ij = derivN_;
					bonds_list_central_1D[tidx*max_number_of_n_neighbours_ + tidy].f_r_ij = valueN_ + derivN_ * (r_ij - xmax_shifted_);
				}
				else
				{
					// Xs_next, Ydelta, Y2_next, Y2, Y_next
					int k = (int)(r_ij / h_);
					double a = values_vec_[k].Xs_next - r_ij;
					double b = h_ - a;
					double asq = a * a;

					double bsq = b * b;

					//	printf("%lf:%lf\t",asq,bsq);
					bonds_list_central_1D[tidx*max_number_of_n_neighbours_ + tidy].fprime_r_ij = values_vec_[k].Ydelta + ((3.0 * bsq - hsq_) * values_vec_[k].Y2_next - (3.0 * asq - hsq_) * values_vec_[k].Y2);
					bonds_list_central_1D[tidx*max_number_of_n_neighbours_ + tidy].f_r_ij = values_vec_[k].Y_next - a * values_vec_[k].Ydelta + ((asq - hsq_) * a * values_vec_[k].Y2 + (bsq - hsq_) * b * values_vec_[k].Y2_next);
				}
			}
	}

}

//evaluate
void CudaSpline::evaluateCUDA(vec3d* n_bonds_1D,meam_bond* bonds_list_central_1D,int* n_num_,int number_of_atoms_,int max_number_of_n_neighbours_)
{
	vec3d *dev_n_bonds_1D;
	meam_bond *dev_bonds_list_central_1D;
	int * dev_n_num_;
	values_struct* dev_values_vec_;

	cudaError_t cudaStatus;

	//alokowanie pamiêci dla values_vec
	cudaStatus = cudaMalloc(&dev_values_vec_, N_*sizeof(values_struct));
	
	RaiseError(cudaStatus);

	//kopiowanie pamiêci dla values_vec
	cudaStatus = cudaMemcpy(dev_values_vec_, values_vec_, N_*sizeof(values_struct), cudaMemcpyHostToDevice);

	RaiseError(cudaStatus);

	//alokowanie pamiêci dla n_bonds_1D
	cudaStatus = cudaMalloc(&dev_n_bonds_1D, number_of_atoms_*max_number_of_n_neighbours_*sizeof(vec3d));

	RaiseError(cudaStatus);

	//kopiowanie pamiêci dla n_bonds_1D
	cudaStatus = cudaMemcpy(dev_n_bonds_1D, n_bonds_1D, number_of_atoms_*max_number_of_n_neighbours_*sizeof(vec3d), cudaMemcpyHostToDevice);

	RaiseError(cudaStatus);

	//alokowanie pamiêci dla bonds_list_central_1D
	cudaStatus = cudaMalloc(&dev_bonds_list_central_1D, max_number_of_n_neighbours_*number_of_atoms_*sizeof(meam_bond));

	RaiseError(cudaStatus);

	cudaStatus = cudaMemcpy(dev_bonds_list_central_1D, bonds_list_central_1D, max_number_of_n_neighbours_*number_of_atoms_*sizeof(meam_bond),cudaMemcpyHostToDevice);

	RaiseError(cudaStatus);

	//alokowanie pamiêci dla n_num
	cudaStatus = cudaMalloc(&dev_n_num_, number_of_atoms_*sizeof(int));
	
	RaiseError(cudaStatus);

	//kopiowanie pamiêci dla n_num
	cudaStatus = cudaMemcpy(dev_n_num_, n_num_, number_of_atoms_*sizeof(int), cudaMemcpyHostToDevice);

	RaiseError(cudaStatus);

	dim3 block(1024, 1024);
	dim3 grid(32,32);

	evaluateKernel << <block, grid >> >(dev_n_bonds_1D, dev_bonds_list_central_1D, dev_n_num_, number_of_atoms_, max_number_of_n_neighbours_, xmin_ , deriv0_, derivN_, value0_, valueN_, xmax_shifted_, dev_values_vec_, h_,hsq_);
	
	//kopiowanie pamiêci dla bonds_list_central_1D
	cudaStatus = cudaMemcpy(bonds_list_central_1D, dev_bonds_list_central_1D, max_number_of_n_neighbours_*number_of_atoms_*sizeof(meam_bond), cudaMemcpyDeviceToHost);

	RaiseError(cudaStatus);

	cudaFree(dev_bonds_list_central_1D);
	cudaFree(dev_n_bonds_1D);
	cudaFree(dev_values_vec_);
	cudaFree(dev_n_num_);
}

//stary evaluate
double CudaSpline::evaluateCUDA(pair* results,double x,double dx,int N)
{
	values_struct* dev_values_vec = 0;
	pair *dev_results = 0;

	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	float elapsedTime = 0;

	cudaEventRecord(start, 0);

	cudaError_t cudaStatus;
	


	cudaStatus = cudaMalloc((void**)&dev_results, N_*sizeof(pair));
	cudaStatus = cudaMalloc((void**)&dev_values_vec, N_*sizeof(values_struct));

	cudaStatus = cudaMemcpy(dev_values_vec, values_vec_, N_*sizeof(values_struct), cudaMemcpyHostToDevice);

	evaluateKernel <<<N_,1>>>(dev_results, N_, x,xmin_, h_, xmax_shifted_, dx, value0_, deriv0_, valueN_, derivN_, hsq_, dev_values_vec);


	cudaStatus = cudaMemcpy(results, dev_results, N_*sizeof(pair), cudaMemcpyDeviceToHost);

	if (cudaStatus != cudaSuccess) 
	{
		fprintf(stderr, "evaluateCUDA failed!");
	}

	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&elapsedTime, start, stop);

	cudaEventDestroy(start);
	cudaEventDestroy(stop);

	cudaFree(dev_results);
	cudaFree(dev_values_vec);

	return (double)elapsedTime;
}

void CudaSpline::initialize(int N, double deriv0, double derivN,
	double xmin, double xmax, double *Y)
{
	int i;

	N_ = N;
	if ((N_ < 3) || (N_ != N_))
		raiseError(ERR_SPLINE, 15, "initialize", "incorrect N");

	deriv0_ = deriv0;
	derivN_ = derivN;

	xmin_ = xmin;
	xmax_ = xmax;
	if (xmin_ >= xmax_)
		raiseError(ERR_SPLINE, 16, "initialize", "incorrect xmin and xmax");

	cutoff_ = xmax_;
	h_ = (xmax_ - xmin_) / (double(N_) - 1.0);
	hsq_ = h_ * h_;

	if (Y == NULL)
		raiseError(ERR_SPLINE, 17, "initialize", "NULL pointer passed as an argument");

	X_ = new double[N_];
	Y_ = new double[N_];
	for (i = 0; i < N_; i++)
	{
		X_[i] = xmin_ + i * h_;
		Y_[i] = Y[i];

		if (Y_[i] != Y_[i])
			raiseError(ERR_SPLINE, 18, "initialize", "incorrect Y[i]");
	}

	init();
}

void CudaSpline::init()
{
	int i, k;
	double sig, p;
	double qn, un;

	xmax_shifted_ = xmax_ - xmin_;

	Xs_ = new double[N_];
	Ydelta_ = new double[N_];
	Y2_ = new double[N_];


	for (i = 0; i < N_; i++)
		Xs_[i] = i * h_;

	double *u = new double[N_];
	Y2_[0] = -0.5;
	u[0] = (3.0 / (X_[1] - X_[0])) * ((Y_[1] - Y_[0]) / (X_[1] - X_[0]) - deriv0_);
	for (i = 1; i < (N_ - 1); i++)
	{
		sig = (X_[i] - X_[i - 1]) / (X_[i + 1] - X_[i - 1]);
		p = sig * Y2_[i - 1] + 2.0;
		Y2_[i] = (sig - 1.0) / p;
		u[i] = (Y_[i + 1] - Y_[i]) / (X_[i + 1] - X_[i]) - (Y_[i] - Y_[i - 1]) / (X_[i] - X_[i - 1]);
		u[i] = (6.0 * u[i] / (X_[i + 1] - X_[i - 1]) - sig * u[i - 1]) / p;
	}

	qn = 0.5;
	un = (3.0 / (X_[N_ - 1] - X_[N_ - 2])) * (derivN_ - (Y_[N_ - 1] - Y_[N_ - 2]) / (X_[N_ - 1] - X_[N_ - 2]));
	Y2_[N_ - 1] = (un - qn *u[N_ - 2]) / (qn * Y2_[N_ - 2] + 1.0);
	for (k = N_ - 2; k >= 0; k--)
		Y2_[k] = Y2_[k] * Y2_[k + 1] + u[k];

	delete[] u;

	for (i = 0; i < N_; i++)
	{
		if (i < (N_ - 1))
			Ydelta_[i] = (Y_[i + 1] - Y_[i]) / h_;
		Y2_[i] /= (h_ * 6.0);
	}

	value0_ = Y_[0];
	valueN_ = Y_[N_ - 1];

#ifdef SPLINE_FAST
	values_vec_ = new values_struct[N_];
	for (i = 0; i < N_ - 1; i++)
	{
		values_vec_[i].Xs_next = Xs_[i + 1];
		values_vec_[i].Y_next = Y_[i + 1];
		values_vec_[i].Ydelta = Ydelta_[i];
		values_vec_[i].Y2 = Y2_[i];
		values_vec_[i].Y2_next = Y2_[i + 1];
	}
#endif

}

void CudaSpline::ShowParams()
{
	std::cout << "N      = " << N_ << std::endl;
	std::cout << "value0 = " << std::setprecision(14) << value0_ << std::endl;
	std::cout << "deriv0 = " << std::setprecision(14) << deriv0_ << std::endl;
	std::cout << "valueN = " << std::setprecision(14) << valueN_ << std::endl;
	std::cout << "derivN = " << std::setprecision(14) << derivN_ << std::endl;
	std::cout << "xmin   = " << std::setprecision(14) << xmin_ << std::endl;
	std::cout << "xmax   = " << std::setprecision(14) << xmax_ << std::endl;
	std::cout << "values: " << std::endl;
	for (int i = 0; i < N_; i++)
		std::cout << std::setprecision(14) << Y_[i] << std::endl;
}

void CudaSpline::readFromASCIIFile(std::ifstream &input)
{
	int i;
	std::string keyword;

	input >> keyword;
	input >> N_;
	if (keyword != "N")
		raiseError(ERR_SPLINE, 1, "readFromASCIIFile", "\"N\" was expected");
	if ((N_ < 3) || (N_ != N_))
		raiseError(ERR_SPLINE, 2, "readFromASCIIFile", "incorrect N");

	input >> keyword;
	input >> deriv0_;
	if (keyword != "deriv0")
		raiseError(ERR_SPLINE, 3, "readFromASCIIFile", "\"deriv0\" was expected");
	if (deriv0_ != deriv0_)
		raiseError(ERR_SPLINE, 4, "readFromASCIIFile", "incorrect deriv0");

	input >> keyword;
	input >> derivN_;
	if (keyword != "derivN")
		raiseError(ERR_SPLINE, 5, "readFromASCIIFile", "\"derivN\" was expected");
	if (derivN_ != derivN_)
		raiseError(ERR_SPLINE, 6, "readFromASCIIFile", "incorrect derivN");

	input >> keyword;
	input >> xmin_;
	if (keyword != "xmin")
		raiseError(ERR_SPLINE, 7, "readFromASCIIFile", "\"xmin\" was expected");
	if (xmin_ != xmin_)
		raiseError(ERR_SPLINE, 8, "readFromASCIIFile", "incorrect xmin");

	input >> keyword;
	input >> xmax_;
	if (keyword != "xmax")
		raiseError(ERR_SPLINE, 9, "readFromASCIIFile", "\"xmax\" was expected");
	if (xmax_ != xmax_)
		raiseError(ERR_SPLINE, 10, "readFromASCIIFile", "incorrect xmax");

	if (xmin_ >= xmax_)
		raiseError(ERR_SPLINE, 11, "readFromASCIIFile", "incorrect xmin and xmax");

	cutoff_ = xmax_;
	h_ = (xmax_ - xmin_) / (double(N_) - 1.0);
	hsq_ = h_ * h_;

	input >> keyword;
	if (keyword != "values")
		raiseError(ERR_SPLINE, 12, "readFromASCIIFile", "\"values\" was expected");

	X_ = new double[N_];
	Y_ = new double[N_];
	for (i = 0; i < N_; i++)
	{
		X_[i] = xmin_ + i * h_;
		input >> Y_[i];

		if (Y_[i] != Y_[i])
			raiseError(ERR_SPLINE, 13, "readFromASCIIFile", "incorrect Y[i]");

		if (input.eof() == 1)
			raiseError(ERR_SPLINE, 14, "readFromASCIIFile", "unexpected end of file");
	}

	init();
}




