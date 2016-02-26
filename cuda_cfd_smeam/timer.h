// written by Szymon Winczewski

#ifndef timer_h
#define timer_h

#include <cstdlib>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <string>
#include <time.h>

#include "errors.h"


class __declspec(dllexport) Timer
{
private:
    int max_number_of_routines_;
    std::string timer_name_;
    bool *routine_defined_;
    std::string *routine_name_;
    long int *routine_executions_;
    double *routine_total_time_;
    double *routine_av_time_;
    clock_t *routine_total_clock_;
    clock_t *routine_started_;
    clock_t *routine_stopped_;

public:
    Timer(std::string timer_name, int max_number_of_routines);
    ~Timer();

    void resetTimings();
    void addRoutine(std::string routine_name, int routine_id);
    void routineStarted(int routine_id);
    void routineStopped(int routine_id);

    void printTimings();
};

#endif
