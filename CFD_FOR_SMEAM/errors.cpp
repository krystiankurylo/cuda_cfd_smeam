// written by Szymon Winczewski

#include "errors.h"

using namespace std;


std::string *module_name = NULL;


void initializeErrors()
{
    module_name = new std::string [1000];

    module_name[ERR_MAIN]                          = "main.cpp";
    module_name[ERR_TIMER]                         = "timer.cpp";
    module_name[ERR_SPLINE]                        = "spline.cpp";
    module_name[ERR_CFD_SMEAM]                     = "cfd_smeam.cpp";
    module_name[ERR_CENTRAL_FORCES_FILE]           = "central_forces_file.cpp";
}


void destroyErrors()
{
    delete [] module_name;
}


void raiseError(int module_id, int error_code, std::string routine_name, std::string error_message)
{
    int full_error_code;
    full_error_code = module_id * 100 + error_code;
    std::cout << std::endl;
    std::cout << "error " << full_error_code;
    std::cout << ", in " << module_name[module_id] << ", routine " << routine_name << ": ";
    std::cout << error_message << std::endl;
    std::cout << std::endl;
    exit(full_error_code);
}


void raiseError(int module_id, int error_code, std::string routine_name, std::string error_message, int value)
{
    int full_error_code;
    full_error_code = module_id * 100 + error_code;
    std::cout << std::endl;
    std::cout << "error " << full_error_code;
    std::cout << ", in " << module_name[module_id] << ", routine " << routine_name << ": ";
    std::cout << error_message;
    std::cout << " " << value;
    std::cout << std::endl;
    exit(full_error_code);
}


void raiseError(int module_id, int error_code, std::string error_message)
{
    int full_error_code;
    full_error_code = module_id * 100 + error_code;
    std::cout << std::endl;
    std::cout << "error " << full_error_code;
    std::cout << ", in " << module_name[module_id] << ": ";
    std::cout << error_message << std::endl;
    std::cout << std::endl;
    exit(full_error_code);
}


void raiseError(int module_id, int error_code)
{
    int full_error_code;
    full_error_code = module_id * 100 + error_code;
    exit(full_error_code);
}


void raiseWarning(int module_id, std::string routine_name, std::string warning_message)
{
    std::cout << std::endl;
    std::cout << "warning";
    std::cout << ", in " << module_name[module_id] << ", routine " << routine_name << ": ";
    std::cout << warning_message << std::endl;
    std::cout << std::endl;
}


void openInputFileError(int module_id, int error_code, std::string routine_name, std::string file_name, std::ifstream &file)
{
    if ( file.is_open() == 0 )
    {
        int full_error_code;
        full_error_code = module_id * 100 + error_code;
        std::cout << std::endl;
        std::cout << "error " << full_error_code;
        std::cout << ", in " << module_name[module_id] << ", routine " << routine_name << ": ";
        std::cout << "could not open input file \"" << file_name << "\"" << std::endl;
        std::cout << std::endl;
        exit(full_error_code);
    }
}


void openOutputFileError(int module_id, int error_code, std::string routine_name, std::string file_name, ofstream &file)
{
    if ( file.is_open() == 0 )
    {
        int full_error_code;
        full_error_code = module_id * 100 + error_code;
        std::cout << std::endl;
        std::cout << "error " << full_error_code;
        std::cout << ", in " << module_name[module_id] << ", routine " << routine_name << ": ";
        std::cout << "could not create output file \"" << file_name << "\"" << std::endl;
        std::cout << std::endl;
        exit(full_error_code);
    }
}
