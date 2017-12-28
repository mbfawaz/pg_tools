#ifndef _COMMON_H_
#define _COMMON_H_

#include <string>
#include <vector>
#include <ostream>
#include <time.h>
#include <cstdio>
#include "sparse_matrix.h"

namespace common
{
	void line2arr (char* str, std::vector<std::string>* arr, char *tokenizer);
    
    void print_error  ( FILE* out, char* errmsg  );
    void print_warning( FILE* out, char* warnmsg );
    char* get_current_local_time_string();
    
    void draw_progress_bar_in_terminal(int len, double fraction);
    
    /* thread-safe versions for multithreaded applications */
	void mt_fprintf(FILE* out, char* buffer);
    void mt_line2arr (char* str, std::vector<std::string>* arr, char *tokenizer);
    void mt_print_error  ( FILE* out, char* errmsg  );
    void mt_print_warning( FILE* out, char* warnmsg );
    
    vector<bool> random_vector(int size, int num_of_ones);
    
    enum color_code 
    {
    	FG_RESET    = 0,
    	FG_BOLD     = 1,
    	FG_BLACK    = 30,
        FG_RED      = 31,
        FG_GREEN    = 32,
        FG_YELLOW   = 33,
        FG_BLUE     = 34,
        FG_MAGENTA  = 35,
        FG_CYAN     = 36,
        FG_LIGHT_GRAY    = 37,
        FG_DEFAULT       = 39,
        FG_LIGHT_RED     = 91,
        FG_LIGHT_GREEN   = 92,
        FG_LIGHT_YELLOW  = 93,
        FG_LIGHT_BLUE    = 94,
        FG_LIGHT_MAGENTA = 95,
        FG_LIGHT_CYAN    = 96,
        FG_WHITE    = 97, 
        
        // use background colors rarely.
        BG_RED      = 41,
        BG_GREEN    = 42,
        BG_BLUE     = 44,
        BG_DEFAULT  = 49
    };
    
    std::ostream& operator<<(std::ostream& os, color_code code);
    std::ostream& operator<<(ostream& os, vector<double> vec);
    std::ostream& operator<<(ostream& os, vector<int> vec);
    std::ostream& operator<<(ostream& os, vector<bool> vec);
   
    class Str2int
    {
    public:
        Str2int(): value_() {}

        inline int value() const { return value_; }
        inline bool operator() (const char* str, size_t len)
        {
            value_ = 0;
            /*int sign = 1;
            if (str[0] == '-') { // handle negative
                sign = -1;
                ++str;
                --len;
            }*/

            switch (len) { // handle up to 10 digits, assume we're 32-bit
                case 10:    value_ += (str[len-10] - '0') * 1000000000;
                case  9:    value_ += (str[len- 9] - '0') * 100000000;
                case  8:    value_ += (str[len- 8] - '0') * 10000000;
                case  7:    value_ += (str[len- 7] - '0') * 1000000;
                case  6:    value_ += (str[len- 6] - '0') * 100000;
                case  5:    value_ += (str[len- 5] - '0') * 10000;
                case  4:    value_ += (str[len- 4] - '0') * 1000;
                case  3:    value_ += (str[len- 3] - '0') * 100;
                case  2:    value_ += (str[len- 2] - '0') * 10;
                case  1:    value_ += (str[len- 1] - '0');
                    //value_ *= sign;
                    return value_ > 0;
                default:
                    return false;
            }
        }
    private:
        int value_;
    };
    
    class ScopedFileLock
    {
        FILE* fp;
        
        public:
        ScopedFileLock(FILE* out):fp(out)
        {
            flockfile(fp);
        }
        
        ~ScopedFileLock()
        {
            funlockfile(fp);
        }

    };
}

#endif
