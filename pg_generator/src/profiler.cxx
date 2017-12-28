#include "profiler.h"
#include <sstream>
#include <fstream>

using namespace std;

void line2arr (char* str, vector<string>* arr, char *tokenizer)
{
	string ts;
	char* tok;
	(*arr).clear();
	tok = strtok(str,tokenizer);
	while ( tok != NULL )
	{
		//printf("%s", tok);
		ts.assign(tok);
		(*arr).push_back(ts);
		tok = strtok(NULL,tokenizer);
	}
}
task_timer_data::task_timer_data( string _desc )
{
	task_count = 1;
	task_description = _desc;
	elapsed.wall   = 0;
	elapsed.user   = 0; 
	elapsed.system = 0;
}

/*-------------------------------------------------------------------------
		Method and static variable definitions for profiler class
  -------------------------------------------------------------------------*/

vector< task_timer_data*> profiler::task_timer_records;

profiler::profiler( string desc )
{
	int idx = -1;
	for (unsigned int i = 0; i < task_timer_records.size(); i++ )
	{
		// hopefully size checking is not an extra step
		if ( desc.size() == task_timer_records[i]->task_description.size() )
		{
			if ( desc == task_timer_records[i]->task_description )
			{
				idx = i;
				break;
			}
		}
	}

	if ( idx == -1 )
	{
		new_task_timer_data = new task_timer_data( desc );
		task_timer_records.push_back( new_task_timer_data );
		timer.start();
	}
	else
	{
		new_task_timer_data = task_timer_records[idx];
		new_task_timer_data->task_count++;
		timer.start();
	}
}

profiler::~profiler()
{
	boost::timer::cpu_times elasped_present = timer.elapsed();
	new_task_timer_data->elapsed.wall += elasped_present.wall;
	new_task_timer_data->elapsed.user += elasped_present.user;
	new_task_timer_data->elapsed.system += elasped_present.system;
}

void profiler::print_all_task_timer_records( FILE * out )
{
	fprintf(out, "\nPRINTING TIME PROFILING DETAILS (TIMER RECORDS)\
	              \n---------------------------------------------\n");
	fprintf(out, "%4s %40s %10s %10s %10s %10s | %7s %10s\n", "S.no", "Task Description",
	 "Wall", "User", "System", "Percent", "#Calls", "Avg CPUt");
	fprintf(out, "%4s %40s %10s %10s %10s %10s | %7s %10s\n", "----", "----------------",
	 "----", "----", "------", "-------", "------", "--------");
	
	for ( unsigned int i = 0; i < task_timer_records.size(); i++ )
	{
		double time_wall    = task_timer_records[i]->elapsed.wall/1.0e9,
		       time_user    = task_timer_records[i]->elapsed.user/1.0e09,
		       time_system  = task_timer_records[i]->elapsed.system/1.0e09;
		
		double percent_use = 100.0*(time_user + time_system)/time_wall;
		
		// print the first two columns
		fprintf(out, "%4d %40s ", i+1, 
		    (char*)task_timer_records[i]->task_description.c_str() );
		    
		// the third, fourth and fifth column are printed based on their values
		if ( time_wall > 3600.0 )    fprintf(out, "%8.4f h ", time_wall/3600.0);
        else if ( time_wall > 60.0 ) fprintf(out, "%8.4f m ", time_wall/60.0);
        else                         fprintf(out, "%8.4f s ", time_wall);
        
        if ( time_user > 3600.0 )    fprintf(out, "%8.4f h ", time_user/3600.0);
        else if ( time_user > 60.0 ) fprintf(out, "%8.4f m ", time_user/60.0);
        else                         fprintf(out, "%8.4f s ", time_user);
        
        if ( time_system > 3600.0 )    fprintf(out, "%8.4f h ", time_system/3600.0);
        else if ( time_system > 60.0 ) fprintf(out, "%8.4f m ", time_system/60.0);
        else                           fprintf(out, "%8.4f s ", time_system);
        
        // finally print percentage use, num_calls and average time taken per call
		fprintf(out, "%10.2f | %7d ", percent_use, (int)task_timer_records[i]->task_count);

		double avg_time_cpu = (time_user+time_system)/(1.0*task_timer_records[i]->task_count);
		if ( avg_time_cpu > 3600.0 )    fprintf(out, "%8.4f h\n", avg_time_cpu/3600.0);
        else if ( avg_time_cpu > 60.0 ) fprintf(out, "%8.4f m\n", avg_time_cpu/60.0);
        else if ( avg_time_cpu > 1 )    fprintf(out, "%8.4f s\n", avg_time_cpu);
        else if ( avg_time_cpu > 1e-3 ) fprintf(out, "%8.4f ms\n", avg_time_cpu/1e-3);
        else if ( avg_time_cpu > 1e-6 ) fprintf(out, "%8.4f us\n", avg_time_cpu/1e-6);
        else                            fprintf(out, "%8.4f ns\n", avg_time_cpu/1e-9);
	}
	fprintf(out, "\n");
}

void profiler::print_all_task_timer_records_terminal()
{
	print_all_task_timer_records( stdout );
	/*printf("\nPRINTING TIME PROFILING DETAILS (TIMER RECORDS)\
	              \n-----------------------------------------------\n");
	printf("%4s %40s %11s %11s %11s %11s\n", "S.no", "Task Description", "Wall",
	"User", "System", "Percent");
	printf("%4s %40s %11s %11s %11s %11s\n", "----", "----------------", "----",
	"----", "------", "-------");
	
	for ( unsigned int i = 0; i < task_timer_records.size(); i++ )
	{
		double time_wall    = task_timer_records[i]->elapsed.wall/1.0e9,
		       time_user    = task_timer_records[i]->elapsed.user/1.0e09,
		       time_system  = task_timer_records[i]->elapsed.system/1.0e09;
		
		double percent_use = 100.0*(time_user + time_system)/time_wall;
		
		// print the first two columns
		printf("%4d %40s ", i+1, 
		    (char*)task_timer_records[i]->task_description.c_str() );
		    
		// the third, fourth and fifth column are printed based on their values
		if ( time_wall > 3600.0 )    printf("%9.4f h ", time_wall/3600.0);
        else if ( time_wall > 60.0 ) printf("%9.4f m ", time_wall/60.0);
        else                         printf("%9.4f s ", time_wall);
        
        if ( time_user > 3600.0 )    printf("%9.4f h ", time_user/3600.0);
        else if ( time_user > 60.0 ) printf("%9.4f m ", time_user/60.0);
        else                         printf("%9.4f s ", time_user);
        
        if ( time_system > 3600.0 )    printf("%9.4f h ", time_system/3600.0);
        else if ( time_system > 60.0 ) printf("%9.4f m ", time_system/60.0);
        else                           printf("%9.4f s ", time_system);
        
        // finally print percentage use
		printf("%11.2f\n", percent_use);
	}
	printf("\n");*/
}


void profiler::clear_all_task_timer_records()
{
	for ( unsigned int i = 0; i < task_timer_records.size(); i++)
	{
		delete task_timer_records[i];
	}
	task_timer_records.clear();
}

/*  Ideally, we would like to have a map here, but since timer records are typically
	less than 5000, a linear search performance is equivalent to a map. We also
	prevent lockup times for updating the map.
*/
const task_timer_data* profiler::get_task_timer_data_by_desc( string desc )
{
	for (unsigned int i = 0; i < task_timer_records.size(); i++ )
	{
		// hopefully size checking is not an extra step
		if ( desc.size() == task_timer_records[i]->task_description.size() )
		{
			if ( desc == task_timer_records[i]->task_description )
			{
				return task_timer_records[i];
				break;
			}
		}
	}
	
	return NULL;
}

/*-------------------------------------------------------------------------
	Method and static variable definitions for mt_profiler class
  -------------------------------------------------------------------------*/

boost::mutex mt_profiler::ttr_mutex;

mt_profiler::mt_profiler( string desc )
{
	new_task_timer_data = new task_timer_data( desc );
	{
		boost::lock_guard<boost::mutex> lock(ttr_mutex);
		task_timer_records.push_back( new_task_timer_data );
	}
	timer.start();
}

mt_profiler::~mt_profiler()
{
	new_task_timer_data->elapsed = timer.elapsed();
}

void mt_profiler::mt_print_all_task_timer_records( FILE * out )
{
	boost::lock_guard<boost::mutex> lock(ttr_mutex);
	print_all_task_timer_records( out );
}

void mt_profiler::mt_clear_all_task_timer_records()
{
	boost::lock_guard<boost::mutex> lock(ttr_mutex);
	clear_all_task_timer_records();
}

const task_timer_data* mt_profiler::mt_get_task_timer_data_by_desc( string desc )
{
	boost::lock_guard<boost::mutex> lock(ttr_mutex);
	return get_task_timer_data_by_desc( desc );
}

void profiler::get_resource_usage( char *buffer, int who )
{
	int VmPeak = 0, VmSize = 0, VmLck = 0, VmHWM = 0, VmRSS = 0, VmData = 0, VmStk = 0, VmExe = 0, VmLib = 0,
	    VmPTE = 0, Threads = 0;
	    
	struct rusage usage;
	int result = getrusage(who, &usage);
	vector<string> arr;

	bool status_available = true;
	ifstream status;
	status.open("/proc/self/status");
	if (!status)
	{
		printf("%s\n", "Could not open file for reading status\n");
		status_available = false;
	}
	if( status_available )
	{
		char str[255];
		while (!status.eof())
		{
			status.getline(str, 255);
			//printf("%s\n", str);
			arr.clear();
			line2arr(str, &arr, (char*)"\t: ");
			//printf("arr[0] = %s, arr[1] = %s\n", arr[0].c_str(), arr[1].c_str());
			if ( arr.empty() )
				continue;
				
			if ( arr[0] == "VmPeak" )
				VmPeak = atoi((char*)arr[1].c_str());
			else if ( arr[0] == "VmSize" )
				VmSize = atoi((char*)arr[1].c_str());
			else if ( arr[0] == "VmLck" )
				VmLck = atoi((char*)arr[1].c_str());
			else if ( arr[0] == "VmHWM" )
				VmHWM = atoi((char*)arr[1].c_str());
			else if ( arr[0] == "VmRSS" )
				VmRSS = atoi((char*)arr[1].c_str());
			else if ( arr[0] == "VmData" )
				VmData = atoi((char*)arr[1].c_str());
			else if ( arr[0] == "VmStk" )
				VmStk = atoi((char*)arr[1].c_str());
			else if ( arr[0] == "VmExe" )
				VmExe = atoi((char*)arr[1].c_str());
			else if ( arr[0] == "VmLib" )
				VmLib = atoi((char*)arr[1].c_str());
			else if ( arr[0] == "VmPTE" )
				VmPTE = atoi((char*)arr[1].c_str());
			else if ( arr[0] == "Threads" )
				Threads = atoi((char*)arr[1].c_str());
		}
	}
	status.close();
	
	double total_user_cpu_time = usage.ru_utime.tv_sec + usage.ru_utime.tv_usec*pow(10,-6),
	       total_syst_cpu_time = usage.ru_stime.tv_sec + usage.ru_stime.tv_usec*pow(10,-6);	
	
	double KiB_to_MiB = 0.000976562,
	       KiB_to_GiB = 0.000976562*0.000976562;
	if ( result != -1 )
	{		
		sprintf(buffer,"\
%40s : %9.4f sec\n\
%40s : %9.4f sec\n\
%40s : %9.2f MiB (%2.4f GiB)\n\
%40s : %9.2f MiB (%2.4f GiB)\n\
%40s : %9.2f MiB (%2.4f GiB)\n\
%40s : %9.2f MiB (%2.4f GiB)\n\
%40s : %9.2f MiB (%2.4f GiB)\n\
%40s : %9.2f MiB (%2.4f GiB)\n\
%40s : %9.2f MiB (%2.4f GiB)\n\
%40s : %9.2f MiB (%2.4f GiB)\n\
%40s : %9.2f MiB (%2.4f GiB)\n\
%40s : %9.2f MiB (%2.4f GiB)\n\
%40s : %9d\n\
%40s : %9ld\n\
%40s : %9ld\n\
%40s : %9ld\n\
%40s : %9ld\n\
%40s : %9ld\n\
%40s : %9ld\n\
%40s : %9ld\n\
%40s : %9ld\n\
%40s : %9ld\n",
		"Total user CPU-time used",                total_user_cpu_time, 
		"Total system CPU-time used",              total_syst_cpu_time,
		"Peak virtual memory size",                KiB_to_MiB*VmPeak,  KiB_to_GiB*VmPeak,
		"Current virtual memory size",             KiB_to_MiB*VmSize,  KiB_to_GiB*VmSize,
		"Current locked memory size",              KiB_to_MiB*VmLck,   KiB_to_GiB*VmLck,
		"Peak resident set size",                  KiB_to_MiB*VmHWM,   KiB_to_GiB*VmHWM,
		"Current resident set size",               KiB_to_MiB*VmRSS,   KiB_to_GiB*VmRSS,
		"Size of Data",                            KiB_to_MiB*VmData,  KiB_to_GiB*VmData,
		"Size of stack",                           KiB_to_MiB*VmStk,   KiB_to_GiB*VmStk,
		"Size of text segments",                   KiB_to_MiB*VmExe,   KiB_to_GiB*VmExe,
		"Shared library code size",                KiB_to_MiB*VmLib,   KiB_to_GiB*VmLib,
		"Page table entry size (af Linux 2.6.10)", KiB_to_MiB*VmPTE,   KiB_to_GiB*VmPTE,
		"Number of Threads",                       Threads,
		"page reclaims (soft page faults) no I/O", usage.ru_minflt,   
		"hard page faults that required I/O",      usage.ru_majflt,
		"Number of Swaps",                         usage.ru_nswap,
		"Number of times the file system inputs",  usage.ru_inblock,
		"Number of times the file system outputs", usage.ru_oublock,
		"Number of IPC messages sent",             usage.ru_msgsnd,
		"Number of IPC messages received",         usage.ru_msgrcv,
		"Number of voluntary context switches",    usage.ru_nvcsw,
		"Number of involuntary context switches",  usage.ru_nivcsw );
	}
	else
	{
	 	printf("There was an error retreiving the stats!!!!!\n");
	 	buffer = NULL;
	}

	//printf("Sorry, but reporting resource usage is not yet supported on this platform.");
	//buffer = NULL;
}
