#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include "errors.h"
#include "defines.h"

#ifndef cuda_spline_h
#define cuda_spline_h

typedef struct
{
	double vec[3];
	double r;
} vec3d;


typedef struct
{
	double f_r_ij;
	int atom_j;
	double r_ij;
	double inv_of_r_ij;
	double r_ij_dir[3];
	double fprime_r_ij;
} meam_bond;

typedef struct
{
	double Xs_next;
	double Ydelta;
	double Y2_next;
	double Y2;
	double Y_next;
} values_struct;

typedef struct
{
	double key;
	double value;
} pair;

class __declspec(dllexport) CudaSpline 
{
private:
	
	double *X_;              // wezly spline'a
	double *Xs_;             // przesuniete wezly spline'a
	double *Y_;              // wartosci funkcji w wezlach
	double *Y2_;             // wartosci drugiej pochodnej w wezlach
	double *Ydelta_;         // okrelsone jako: Ydelta[i] = (Y[i+1]-Y[i])/h
	int N_;                  // liczba wezlow spline'a
	double deriv0_;          // wartosc pierwszej pochodnej w wezle 0
	double value0_;          // wartosc funkcji w wezle 0
	double derivN_;          // wartosc pierwszej pochodnej w wezle (N-1)
	double valueN_;          // wartosc funkcji pochodnej w wezle (N-1)
	double xmin_;            // poczatek przedzialu na ktorym okreslony jest spline
	double xmax_;            // konie przedzialu na ktorym okreslony jest spline
	double cutoff_;
	double h_;               // odleglosc pomiedzy wezlami
	double hsq_;             // kwadrat odleglosci pomiedzy wezlami
	double xmax_shifted_;    // koniec przedzialu na ktorym okreslony jest spline po przesunieciu jego poczatku do 0

#ifdef SPLINE_FAST
	values_struct *values_vec_;
#endif


	void init();


public:

	
	
	CudaSpline();
	~CudaSpline();
	double getCutoff();
	void readFromASCIIFile(std::ifstream &file);
	void initialize(int, double, double, double, double, double*);
	void ShowParams();
	void evaluateCUDA(vec3d*,meam_bond*,int*,int,int);
	void EvaluateDerivCUDA(double*, double*, int);
	double evaluateCUDA(pair*,double,double,int);
	inline double evaluate_deriv(double x)
	{
		x -= xmin_;
		if (x <= 0.0)
			return deriv0_;
		else if (x >= xmax_shifted_)
			return derivN_;
		else
		{
			int k = (int)(x / h_);
			double a = values_vec_[k].Xs_next - x;
			double b = h_ - a;
			double asq = a * a;
			double bsq = b * b;
			return values_vec_[k].Ydelta + ((3.0 * bsq - hsq_) * values_vec_[k].Y2_next - (3.0 * asq - hsq_) * values_vec_[k].Y2);
		}
	}

#ifndef SPLINE_FAST
	inline double evaluate(double x)
	{
		x -= xmin_;
		if (x <= 0.0)
			return Y_[0] + deriv0_ * x;
		else if (x >= xmax_shifted_)
			return Y_[N_ - 1] + derivN_ * (x - xmax_shifted_);
		else
		{
			int klo = (int)(x / h_);
			int khi = klo + 1;
			double a = Xs_[khi] - x;
			double b = h_ - a;
			return Y_[khi] - a * Ydelta_[klo] + ((a * a - hsq_) * a * Y2_[klo] + (b * b - hsq_) * b * Y2_[khi]);
		}
	}
#else
    inline double evaluate(double x)
	{
		x -= xmin_;
		if (x <= 0.0)
			return value0_ + deriv0_ * x;
		else if (x >= xmax_shifted_)
			return valueN_ + derivN_ * (x - xmax_shifted_);
		else
		{
			int k = (int)(x / h_);
			double a = values_vec_[k].Xs_next - x;
			double b = h_ - a;
			return values_vec_[k].Y_next - a * values_vec_[k].Ydelta + ((a * a - hsq_) * a * values_vec_[k].Y2 + (b * b - hsq_) * b * values_vec_[k].Y2_next);
		}
	}
#endif


#ifndef SPLINE_FAST
	inline double evaluate(double x, double &deriv)
	{
		x -= xmin_;
		if (x <= 0.0)
		{
			deriv = deriv0_;
			return Y_[0] + deriv0_ * x;
		}
		else if (x >= xmax_shifted_)
		{
			deriv = derivN_;
			return Y_[N_ - 1] + derivN_ * (x - xmax_shifted_);
		}
		else
		{
			int klo = (int)(x / h_);
			int khi = klo + 1;
			double a = Xs_[khi] - x;
			double b = h_ - a;
			double asq = a * a;
			double bsq = b * b;
			deriv = Ydelta_[klo] + ((3.0 * bsq - hsq_) * Y2_[khi] - (3.0 * asq - hsq_) * Y2_[klo]);
			return Y_[khi] - a * Ydelta_[klo] + ((asq - hsq_) * a * Y2_[klo] + (bsq - hsq_) * b * Y2_[khi]);
		}
	}
#else
	inline double evaluate(double x, double &deriv)
	{
		x -= xmin_;
		if (x <= 0.0)
		{
			deriv = deriv0_;
			return value0_ + deriv0_ * x;
		}
		else if (x >= xmax_shifted_)
		{
			deriv = derivN_;
			return valueN_ + derivN_ * (x - xmax_shifted_);
		}
		else
		{
			// Xs_next, Ydelta, Y2_next, Y2, Y_next
			int k = (int)(x / h_);
			double a = values_vec_[k].Xs_next - x;
			double b = h_ - a;
			double asq = a * a;
			double bsq = b * b;
			//printf("From Host\n");
			//printf("%lf:%lf\n", asq, bsq);
			deriv = values_vec_[k].Ydelta + ((3.0 * bsq - hsq_) * values_vec_[k].Y2_next - (3.0 * asq - hsq_) * values_vec_[k].Y2);
			return values_vec_[k].Y_next - a * values_vec_[k].Ydelta + ((asq - hsq_) * a * values_vec_[k].Y2 + (bsq - hsq_) * b * values_vec_[k].Y2_next);
		}
	}
#endif
};

#endif
