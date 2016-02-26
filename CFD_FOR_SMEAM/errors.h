// written by Szymon Winczewski

#ifndef errors_h
#define errors_h

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

extern std::string *module_name;


const int ERR_MAIN                          =  0;
const int ERR_TIMER                         = 10;
const int ERR_SPLINE                        = 20;
const int ERR_CFD_SMEAM                     = 30;
const int ERR_CENTRAL_FORCES_FILE           = 40;


void initializeErrors();
void destroyErrors();
void raiseError(int module_id, int error_code, std::string routine_name, std::string error_message);
void raiseError(int module_id, int error_code, std::string routine_name, std::string error_message, int value);
void raiseError(int module_id, int error_code, std::string error_message);
void raiseError(int module_id, int error_code);
void raiseWarning(int module_id, std::string routine_name, std::string warning_message);

void openInputFileError(int module_id, int error_code, std::string routine_name, std::string file_name, std::ifstream &file);
void openOutputFileError(int module_id, int error_code, std::string routine_name, std::string file_name, std::ofstream &file);

#endif
