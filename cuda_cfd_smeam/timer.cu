// written by Szymon Winczewski

#include "timer.h"

using namespace std;


Timer::Timer(std::string timer_name, int max_number_of_routines)
{
    if ( timer_name == "" )
        raiseError(ERR_TIMER, 1, "Timer", "wrong timer name specified");
    timer_name_ = timer_name;

    if ( max_number_of_routines < 0 )
        raiseError(ERR_TIMER, 2, "Timer", "incorrect maximum number of routines parameter");
    max_number_of_routines_ = max_number_of_routines;

    routine_defined_     = new bool [max_number_of_routines_];
    for (int i = 0; i < max_number_of_routines_; i++)
        routine_defined_[i] = 0;
    routine_name_        = new std::string[max_number_of_routines_];
    routine_executions_  = new long int [max_number_of_routines_];
    routine_total_time_  = new double [max_number_of_routines_];
    routine_av_time_     = new double [max_number_of_routines_];
    routine_total_clock_ = new clock_t [max_number_of_routines_];
    routine_started_     = new clock_t [max_number_of_routines_];
    routine_stopped_     = new clock_t [max_number_of_routines_];
}


Timer::~Timer()
{
    delete [] routine_defined_;
    delete [] routine_name_;
    delete [] routine_executions_;
    delete [] routine_total_time_;
    delete [] routine_av_time_;
    delete [] routine_total_clock_;
    delete [] routine_started_;
    delete [] routine_stopped_;
}


void Timer::resetTimings()
{
    for (int i = 0; i < max_number_of_routines_; i++)
    {
        routine_executions_[i] = 0;
        routine_total_time_[i] = 0.0;
        routine_av_time_[i] = 0.0;
        routine_total_clock_[i] = 0;
    }
}


void Timer::addRoutine(std::string routine_name, int routine_id)
{
    if ( ( routine_id < 0 ) || ( routine_id > max_number_of_routines_ ) )
        raiseError(ERR_TIMER, 3, "addRoutine", "wrong routine id");
    if ( routine_defined_[routine_id] == 1 )
        raiseError(ERR_TIMER, 4, "addRoutine", "routine already defined");

    routine_defined_[routine_id] = 1;
    routine_name_[routine_id] = routine_name;
    routine_executions_[routine_id] = 0;
    routine_total_time_[routine_id] = 0.0;
    routine_av_time_[routine_id] = 0.0;
    routine_total_clock_[routine_id] = 0;
}


void Timer::routineStarted(int routine_id)
{
    routine_started_[routine_id] = clock();
}


void Timer::routineStopped(int routine_id)
{
    routine_stopped_[routine_id] = clock();
    routine_total_clock_[routine_id] += ( routine_stopped_[routine_id] - routine_started_[routine_id] );
    routine_executions_[routine_id]++;
}


void Timer::printTimings()
{
    std::cout << std::endl;
    std::cout << "Timer: " << timer_name_ << std::endl;
    for (int i = 0; i < max_number_of_routines_; i++)
    {
        if ( routine_defined_[i] == 1 )
        {
            routine_total_time_[i] = double(routine_total_clock_[i]) / double(CLOCKS_PER_SEC);
            if ( routine_executions_[i] != 0 )
                routine_av_time_[i] = double(routine_total_time_[i]) / double(routine_executions_[i]);
            else
                routine_av_time_[i] = 0.0;

            std::cout << "   " << setw(50) << routine_name_[i] << " / "
                          << setw(10) << routine_executions_[i] << " / "
                          << setw(15)  << setprecision(9) << routine_total_time_[i] << " / "
                          << setw(15) << setprecision(9) << routine_av_time_[i] << " /" << std::endl;
        }
    }
    std::cout << std::endl;
    std::cout << std::endl;
}
