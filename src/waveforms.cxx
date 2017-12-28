#include "waveforms.h"
#include "math.h"
#include <algorithm>
#include <iostream>

using namespace std;

// Default constructor
waveforms::waveforms()
{
    cholmod_start(&c);
}

// Constructor
waveforms::waveforms(int _SOURCES, int _NUM_TIME_POINTS)
{
    cholmod_start(&c);
    SOURCES          = _SOURCES;
    NUM_TIME_POINTS  = _NUM_TIME_POINTS;
    
}
        
// Copy constructor
waveforms::waveforms(waveforms & to_copy)
{
    cholmod_start(&c);
}
    
// Overload the equality operator 
waveforms & waveforms::operator= (waveforms & to_copy)
{
    if (this == &to_copy)
        return (*this);

    return *this;
}

// Generate the waveforms
void waveforms::generate (int    _SOURCES,
                          int    _NUM_TIME_POINTS,
                          double _MIN_TIME_STEP,
                          double _MAX_TIME_STEP,
                          double _MIN_I,
                          double _MAX_I,
                          int    _RANDOM_SEED)
{
    // Set up dimensions and parameters
    SOURCES = _SOURCES;
    NUM_TIME_POINTS = _NUM_TIME_POINTS;
    MIN_TIME_STEP   = _MIN_TIME_STEP;
    MAX_TIME_STEP   = _MAX_TIME_STEP;
    MIN_I           = _MIN_I;
    MAX_I           = _MAX_I;
    
    // Parameters
    int m  = SOURCES;
    int N  = NUM_TIME_POINTS;

    // Seed for the RNG
    srand(_RANDOM_SEED);

    // Allocating memory
    current_vectors.resize(N);
    for (int i = 0; i<N; i++)
        current_vectors[i].resize(m);
    time_points.resize(N); 

    // Time points
    time_points[0] = 0;
    for (int i = 1; i<N; i++)
        time_points[i] = time_points[i-1] + _MIN_TIME_STEP + (_MAX_TIME_STEP - _MIN_TIME_STEP)*1.0*rand()/RAND_MAX;

    // Two different approaches here:

    // Current values (APPROACH 1)
//     for (int k = 0; k<N; k++)
//         for (int j = 0; j<m; j++)
//             current_vectors[k][j] = MIN_I + (_MAX_I - MIN_I)*1.0*rand()/RAND_MAX;
   
    // Current values (APPROACH 2)
    for (int k = 0; k<N/4; k++)
        for (int j = 0; j<m; j++)
            current_vectors[k][j] = _MIN_I + (_MAX_I - _MIN_I)*1.0*rand()/RAND_MAX;
    _MAX_I *= 1.2;
    for (int k = N/4; k<N/2; k++)
        for (int j = 0; j<m; j++)
            current_vectors[k][j] = _MIN_I + (_MAX_I - _MIN_I)*1.0*rand()/RAND_MAX;
    _MAX_I *= 1.1;
    for (int k = N/2; k<3*N/4; k++)
        for (int j = 0; j<m; j++)
            current_vectors[k][j] = _MIN_I + (_MAX_I - _MIN_I)*1.0*rand()/RAND_MAX;
    _MAX_I *= 0.9;
    for (int k = 3*N/4; k<N; k++)
        for (int j = 0; j<m; j++)
            current_vectors[k][j] = _MIN_I + (_MAX_I - _MIN_I)*1.0*rand()/RAND_MAX;
}

void waveforms::cholmod_currents(vector<cholmod_dense* > &currents, int THREADS, int NODES, vector<int> & csi)
{   
    int TASK_LOAD = NUM_TIME_POINTS/THREADS;
    for (int tid = 0; tid<THREADS; tid++)
    {
        cholmod_dense * new_current_vectors = cholmod_zeros(NODES, TASK_LOAD, CHOLMOD_REAL, &c);
        double * cvx = (double*)new_current_vectors->x;
        for (int j = 0; j<TASK_LOAD; j++)
        {
            for (int  q = 0; q<SOURCES; q++)
            {
                cvx[j*NODES + csi[q]] = current_vectors[tid*TASK_LOAD + j][q];
            }
        }
        currents.push_back(new_current_vectors);
    }
}


// Destroyer
waveforms::~waveforms()
{
    ;
}












