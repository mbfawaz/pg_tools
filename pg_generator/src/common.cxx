#include "common.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <iostream>
using namespace std;

void common::line2arr (char* str, vector<string>* arr, char *tokenizer)
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

ostream& common::operator<<(ostream& os, color_code code)
{
    return os << "\033[" << static_cast<int>(code) << "m";
}

//! This method overloads the output operator<<
//! @param os [in]: ostream
//! @param vector<double> [in]: vector of doubles
ostream& common::operator<<(ostream& os, vector<double> vec)
{
    for(unsigned int i = 0 ; i < vec.size() ; i++)
    {
        cout<<vec[i]<<" ";
    }
    cout<<endl;
    return os;
}

//! This method overloads the output operator<<
//! @param os [in]: ostream
//! @param vector<int> [in]: vector of integers
ostream& common::operator<<(ostream& os, vector<int> vec)
{
    for(unsigned int i = 0 ; i < vec.size() ; i++)
    {
        cout<<vec[i]<<" ";
    }
    cout<<endl;
    return os;
}

//! This method overloads the output operator<<
//! @param os [in]: ostream
//! @param vector<int> [in]: vector of booleans
ostream& common::operator<<(ostream& os, vector<bool> vec)
{
    for(unsigned int i = 0 ; i < vec.size() ; i++)
    {
        cout<<vec[i]<<" ";
    }
    cout<<endl;
    return os;
}


char* common::get_current_local_time_string()
{
    time_t rawtime;
  	struct tm * timeinfo;
  	time ( &rawtime );
  	timeinfo = localtime ( &rawtime );
  	char* time_str = new char[100];
  	strftime ( time_str, 100, "%c", timeinfo);
  	return time_str;
}

void common::print_error( FILE* out, char* errmsg  )
{
    char* time_str = get_current_local_time_string();
    fprintf(out, "\n\033[1;%dm---------- %s --------------\n", FG_RED, time_str);
    fprintf(out, "ERROR!!! %s \033[%dm\n\n", errmsg, FG_RESET);
    delete[] time_str;
} 

void common::print_warning( FILE* out, char* warnmsg )
{
    char* time_str = get_current_local_time_string();
    fprintf(out, "\n\033[1;%dm---------- %s --------------\n", FG_MAGENTA, time_str);
    fprintf(out, "WARNING: %s \033[%dm\n", warnmsg, FG_RESET);
    delete[] time_str;
}

void common::draw_progress_bar_in_terminal(int len, double fraction)
{
    char *progress = new char[len];
	int progress_mark = static_cast<int>( len*fraction ) ;
	if ( progress_mark > len )	progress_mark = len;

	for (int i = 0; i < len; i++)
	{
		if (i < progress_mark)
		{
			progress[i] = '=';
		}
		else 
			progress[i] = ' ';
	}
	fprintf(stderr, "\033[1;%dm[%s\b]\033[0m (%3.2f%%)\r",
	FG_LIGHT_BLUE, progress, 100.0*fraction);
	delete[] progress;
}

void common::mt_line2arr (char* str, vector<string>* arr, char *tokenizer)
{
	string ts;
	char* tok;
	char* saveptr;
	(*arr).clear();
	tok = strtok_r(str,tokenizer, &saveptr);
	while ( tok != NULL )
	{
		//printf("%s", tok);
		ts.assign(tok);
		(*arr).push_back(ts);
		tok = strtok_r(NULL, tokenizer, &saveptr);
	}
}

void common::mt_fprintf(FILE* out, char* buffer)
{
    common::ScopedFileLock sfl(out);
    fprintf(out, "%s", buffer);
}

void common::mt_print_error( FILE* out, char* errmsg  )
{
    common::ScopedFileLock sfl(out);
    common::print_error( out, errmsg );
}

void common::mt_print_warning( FILE* out, char* warnmsg )
{
    common::ScopedFileLock sfl(out);
    common::print_warning( out, warnmsg );
}

//! This method generates a random vector of booleans with specific number of ones
//! @param size [in]: size of vector
//! @param num_of_ones [in]: number of 1s in the vector
//! \return vector<bool>
vector<bool> common::random_vector(int size, int num_of_ones)
{
    vector<bool> vec;
    vec.resize(size,0);
    for(int i = 0; i < num_of_ones ; i++)
    {
        vec[i] = 1;
    }
    random_shuffle(vec.begin(),vec.end());
    
    return vec;
}
