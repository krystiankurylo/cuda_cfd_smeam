// written by Szymon Winczewski
// partially based on pair_meam_spline library (LAMMPS)

#include "spline.h"

using namespace std;


Spline::Spline()
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


Spline::~Spline()
{
    if ( X_ != NULL )
        delete [] X_;
    if ( Xs_ != NULL )
        delete [] Xs_;
    if ( Y_ != NULL )
        delete [] Y_;
    if ( Y2_ != NULL )
        delete [] Y2_;
    if ( Ydelta_ != NULL )
        delete [] Ydelta_;
    #ifdef SPLINE_FAST
    if ( values_vec_ != NULL )
        delete [] values_vec_;
    #endif
}


void Spline::readFromASCIIFile(std::ifstream &input)
{
    int i;
    std::string keyword;

    input >> keyword;
    input >> N_;
    if ( keyword != "N" )
        raiseError(ERR_SPLINE, 1, "readFromASCIIFile", "\"N\" was expected");
    if ( ( N_ < 3 ) || ( N_ != N_ ) )
        raiseError(ERR_SPLINE, 2, "readFromASCIIFile", "incorrect N");

    input >> keyword;
    input >> deriv0_;
    if ( keyword != "deriv0" )
        raiseError(ERR_SPLINE, 3, "readFromASCIIFile", "\"deriv0\" was expected");
    if ( deriv0_ != deriv0_ )
        raiseError(ERR_SPLINE, 4, "readFromASCIIFile", "incorrect deriv0");

    input >> keyword;
    input >> derivN_;
    if ( keyword != "derivN" )
        raiseError(ERR_SPLINE, 5, "readFromASCIIFile", "\"derivN\" was expected");
    if ( derivN_ != derivN_ )
        raiseError(ERR_SPLINE, 6, "readFromASCIIFile", "incorrect derivN");

    input >> keyword;
    input >> xmin_;
    if ( keyword != "xmin" )
        raiseError(ERR_SPLINE, 7, "readFromASCIIFile", "\"xmin\" was expected");
    if ( xmin_ != xmin_ )
        raiseError(ERR_SPLINE, 8, "readFromASCIIFile", "incorrect xmin");

    input >> keyword;
    input >> xmax_;
    if ( keyword != "xmax" )
        raiseError(ERR_SPLINE, 9, "readFromASCIIFile", "\"xmax\" was expected");
    if ( xmax_ != xmax_ )
        raiseError(ERR_SPLINE, 10, "readFromASCIIFile", "incorrect xmax");

    if ( xmin_ >= xmax_ )
        raiseError(ERR_SPLINE, 11, "readFromASCIIFile", "incorrect xmin and xmax");

    cutoff_ = xmax_;
    h_ = ( xmax_ - xmin_ ) / ( double(N_) - 1.0 );
    hsq_ = h_ * h_;

    input >> keyword;
    if ( keyword != "values" )
        raiseError(ERR_SPLINE, 12, "readFromASCIIFile", "\"values\" was expected");

    X_ = new double [N_];
    Y_ = new double [N_];
    for (i = 0; i < N_; i++)
    {
        X_[i] = xmin_ + i * h_;
        input >> Y_[i];

        if ( Y_[i] != Y_[i] )
            raiseError(ERR_SPLINE, 13, "readFromASCIIFile", "incorrect Y[i]");

        if ( input.eof() == 1 )
            raiseError(ERR_SPLINE, 14, "readFromASCIIFile", "unexpected end of file");
    }

    init();
}


void Spline::initialize(int N, double deriv0, double derivN,
                        double xmin, double xmax, double *Y)
{
    int i;

    N_ = N;
    if ( ( N_ < 3 ) || ( N_ != N_ ) )
        raiseError(ERR_SPLINE, 15, "initialize", "incorrect N");

    deriv0_ = deriv0;
    derivN_ = derivN;

    xmin_ = xmin;
    xmax_ = xmax;
    if ( xmin_ >= xmax_ )
        raiseError(ERR_SPLINE, 16, "initialize", "incorrect xmin and xmax");

    cutoff_ = xmax_;
    h_ = ( xmax_ - xmin_ ) / ( double(N_) - 1.0 );
    hsq_ = h_ * h_;

    if ( Y == NULL )
        raiseError(ERR_SPLINE, 17, "initialize", "NULL pointer passed as an argument");

    X_ = new double [N_];
    Y_ = new double [N_];
    for (i = 0; i < N_; i++)
    {
        X_[i] = xmin_ + i * h_;
        Y_[i] = Y[i];

        if ( Y_[i] != Y_[i] )
            raiseError(ERR_SPLINE, 18, "initialize", "incorrect Y[i]");
    }

    init();
}


void Spline::print()
{
    int i;

    std::cout << "N      = " << N_ << std::endl;
    std::cout << "value0 = " << setprecision(14) << value0_ << std::endl;
    std::cout << "deriv0 = " << setprecision(14) << deriv0_ << std::endl;
    std::cout << "valueN = " << setprecision(14) << valueN_ << std::endl;
    std::cout << "derivN = " << setprecision(14) << derivN_ << std::endl;
    std::cout << "xmin   = " << setprecision(14) << xmin_ << std::endl;
    std::cout << "xmax   = " << setprecision(14) << xmax_ << std::endl;
    std::cout << "values: " << std::endl;
    for (i = 0; i < N_; i++)
        std::cout << setprecision(14) << Y_[i] << std::endl;
}


void Spline::writeToMEAMFile(ofstream &output)
{
    int i;
    double x;

    h_ = ( xmax_ - xmin_ ) / ( double(N_) - 1.0 );
    output << setprecision(0) << N_ << std::endl;
    output << setprecision(24) << deriv0_ << "   " << setprecision(24) << derivN_ << std::endl;
    output << "0 0 0 0" << std::endl;
    for (i = 0; i < N_; i++)
    {
        x = xmin_ + i * h_;
        output << setprecision(24) << x << "   " << setprecision(24) << Y_[i] << "   " << "0.0000" << std::endl;
    }
}


double Spline::getCutoff()
{
    return cutoff_;
}


void Spline::init()
{
    int i, k;
    double sig, p;
    double qn, un;

    xmax_shifted_ = xmax_ - xmin_;

    Xs_ = new double [N_];
    Ydelta_ = new double [N_];
    Y2_ = new double [N_];

    for (i = 0; i < N_; i++)
        Xs_[i] = i * h_;

    double *u = new double [N_];
    Y2_[0] = -0.5;
    u[0] = ( 3.0 / ( X_[1] - X_[0] ) ) * ( ( Y_[1] - Y_[0] ) / ( X_[1] - X_[0] ) - deriv0_ );
    for (i = 1; i < ( N_ - 1 ); i++)
    {
        sig = ( X_[i] - X_[i - 1] ) / ( X_[i + 1] - X_[i - 1] );
        p = sig * Y2_[i - 1] + 2.0;
        Y2_[i] = ( sig - 1.0 ) / p;
        u[i] = ( Y_[i + 1] - Y_[i] ) / ( X_[i + 1] - X_[i] ) - ( Y_[i] - Y_[i - 1] ) / ( X_[i] - X_[i - 1] );
        u[i] = ( 6.0 * u[i] / ( X_[i + 1] - X_[i - 1] ) - sig * u[i - 1] ) / p;
    }

    qn = 0.5;
    un = ( 3.0 / ( X_[N_ - 1] - X_[N_ - 2] ) ) * ( derivN_ - ( Y_[N_ - 1] - Y_[N_ - 2] ) / ( X_[N_ - 1] - X_[N_ - 2] ) );
    Y2_[N_ - 1] = ( un - qn *u[N_ - 2] ) / ( qn * Y2_[N_ - 2] + 1.0 );
    for (k = N_ - 2; k >= 0; k--)
        Y2_[k] = Y2_[k] * Y2_[k + 1] + u[k];

    delete [] u;

    for (i = 0; i < N_; i++)
    {
        if ( i < ( N_ - 1 ) )
            Ydelta_[i] = ( Y_[i + 1] - Y_[i] ) / h_;
        Y2_[i] /= ( h_ * 6.0 );
    }

    value0_ = Y_[0];
    valueN_ = Y_[N_ - 1];

    #ifdef SPLINE_FAST
    values_vec_ = new values_struct [N_];
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


void Spline::writeTestFile(std::string file_name)
{
    ofstream file;
    double x;
    double y, y2;
    double y_alt, y2_alt;
    double xrange = xmax_ - xmin_;
    double dx = xrange / 200.0;

    file.open(file_name.c_str());
    openOutputFileError(91, 15, "writeTestFile", file_name, file);
    for (x = xmin_ - 0.5 * xrange; x < xmax_ + 0.5 * xrange; x += dx)
    {
        y = evaluate(x, y2);
        y_alt = evaluate(x);
        y2_alt = evaluate_deriv(x);

        file << x << "   " << y << "   " << y2 << "   " << y_alt << "   " << y2_alt << std::endl;
    }
    file.close();
}
