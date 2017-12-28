#ifndef _PROFILER_H
#define _PROFILER_H

#include <boost/thread.hpp>
#include <boost/timer/timer.hpp>
#include <vector>
#include <iterator>
#include <cstdio>
#include <sys/time.h>
#include <sys/resource.h>

/*!
@struct  	task_timer_data
@author		Sandeep Chatterjee
@details	Stores information related to time-profiling of code sections. Each
            instance is should be identifiable by a unique key, preferably
            describing the purpose of the profiled code segment.
*/
struct task_timer_data
{
	/*! Can be a function name, or may describe a part of a function doing a
	    specific task. This string uniquely identifies the task.*/
	std::string task_description;
	
	/*! Counts the number of times this specific task has been called. */
	size_t task_count;

	/*! A struct provided by boost that stores the wall clock time, system time
	    and user time spent in doing a task. All times are in ns (1e-9 secs). The 
	    definition is as follows:
    \code{.cpp} typedef boost::int_least64_t nanosecond_type;
    struct cpu_times
    {
      nanosecond_type wall;
      nanosecond_type user;
      nanosecond_type system;

      void clear();
    };
    \endcode
	*/
	boost::timer::cpu_times elapsed;
	
	/*! The constructor
	@param[in] _desc The name (that acts as the key) to identify the time-profiler
	records
	*/
	task_timer_data( std::string _desc );
};

/*!
@struct  	profiler
@author		Sandeep Chatterjee
@details	This utility class provides methods that help in time-profiling a code.
*/
class profiler
{
protected:
	task_timer_data* new_task_timer_data;
	boost::timer::cpu_timer timer;
	
	// static vector to hold records
	static std::vector< task_timer_data* > task_timer_records;
	
	profiler(){}

public:	
	profiler( std::string desc );
	
	~profiler();
	
	static const task_timer_data* get_task_timer_data_by_desc( std::string desc );
	static void print_all_task_timer_records( FILE * out );
	static void print_all_task_timer_records_terminal();
	static void clear_all_task_timer_records();
	static void get_resource_usage( char *buffer, int who );
};

class mt_profiler:protected profiler
{
	// additional class static mutex to serialize access to 
	// mt_profiler::task_timer_records, Could be a shared lock in future so that
	// mutiple threads can access time_data by description.
	static boost::mutex ttr_mutex;
	
public:	
	mt_profiler( std::string desc );
	
	~mt_profiler();
	
	static const task_timer_data* mt_get_task_timer_data_by_desc( std::string desc );
	static void mt_print_all_task_timer_records( FILE * out );
	static void mt_clear_all_task_timer_records();
};
#endif
