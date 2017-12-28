#include "container.h"
#include "math.h"
#include <algorithm>
#include <iostream>

using namespace std;

// Default constructor
container::container()
{
    cholmod_start(&c);
}

// Constructor
container::container(int _DIMENSION, int _CONSTRAINTS)
{
    cholmod_start(&c);
    DIMENSION          = _DIMENSION;
    CONSTRAINTS        = _CONSTRAINTS;
    lb.resize(DIMENSION);
    gb.resize(DIMENSION);
    
}
        
// Copy constructor
container::container(container & to_copy)
{
    cholmod_start(&c);
    DIMENSION           = to_copy.DIMENSION;
    CONSTRAINTS         = to_copy.CONSTRAINTS;
    lb                  = to_copy.lb;
    gb                  = to_copy.gb;
    U                   = to_copy.U;
}
    
// Overload the equality operator 
container & container::operator= (container & to_copy)
{
    if (this == &to_copy)
        return (*this);

    DIMENSION           = to_copy.DIMENSION;
    CONSTRAINTS         = to_copy.CONSTRAINTS;
    lb                  = to_copy.lb;
    gb                  = to_copy.gb;
    U                   = to_copy.U;
    return *this;
}

// Generate the container
void container::generate(int    _DIMENSION,
                         int    _CONSTRAINTS,
                         double _MAX_IL,
                         int    _RANDOM_SEED)
{
    // Set up dimensions
    DIMENSION = _DIMENSION;
    CONSTRAINTS = _CONSTRAINTS;
    lb.resize(DIMENSION);
    gb.resize(CONSTRAINTS);

    // We are going to build the transpose of the global constraints matrix first
    // this is way more efficient
 
    // RANDOM SEED
    srand(_RANDOM_SEED);

    // Define the local constraints randomly
    for (int i =0; i<DIMENSION; i++)
        lb[i] = _MAX_IL*((double)rand()/RAND_MAX);

    // Define the global constraint vector
    // First, fill in the first column with all ones (i.e. all current sources are involved)
    vector<int> UTpv, UTiv;
    vector<double> UTxv;
    int nz = DIMENSION; // number of nz in the global constraints matrix
    UTpv.push_back(0);
    for (int i = 0; i<DIMENSION; i++)
    {
        UTxv.push_back(1);
        UTiv.push_back(i);
    }
    UTpv.push_back(DIMENSION);
    
    // Fill in the rest of the rows
    for (int i = 1; i<CONSTRAINTS; i++)
    {
        // Get the number of sources to add in the ith global constraint
        int num_cs = floor(2+0.3*((double)rand()/RAND_MAX)*DIMENSION);

        // Get the required indices
        vector<int> idx;
        for (int j =0; j<num_cs; j++)
            idx.push_back((int)(DIMENSION*((double)rand()/RAND_MAX)));

        // Remove duplicates
        sort(idx.begin(), idx.end());
        idx.erase(unique(idx.begin(), idx.end()), idx.end());
        
        // Fill in the columns of the global constraints matrix
        for (unsigned int j = 0; j<idx.size(); j++)
        {
            UTxv.push_back(1);
            UTiv.push_back(idx[j]);
        }
        UTpv.push_back(nz + idx.size());
        nz += idx.size();
    }

    // Fill in the actual matrix vectors of U
    cholmod_sparse * UT = cholmod_allocate_sparse(DIMENSION, CONSTRAINTS, nz, 0, 1, 0, CHOLMOD_REAL, &c); 
    int    * UTp = (int*)   UT->p;
    int    * UTi = (int*)   UT->i;
    double * UTx = (double*)UT->x;

    for (unsigned int i = 0; i<UTpv.size(); i++)
        UTp[i] = UTpv[i];
    for (unsigned int i = 0; i<UTiv.size(); i++)
    {   
        UTi[i] = UTiv[i];
        UTx[i] = UTxv[i];
    }
       
    // Find the global constraints vectors
    for (int i = 0; i<CONSTRAINTS; i++)
    {
        double sum_il = 0; 
        int start_index = UTp[i];
        int end_index   = UTp[i+1];
        for (int index = start_index; index < end_index; index++)
            sum_il += lb[UTi[index]];
        gb[i] = max( (0.5 + 0.2*((double)rand()/RAND_MAX)) * sum_il, 
                     _MAX_IL + (sum_il - _MAX_IL)*((double)rand()/RAND_MAX));
    }
    U = cholmod_allocate_sparse(UT->ncol, UT->nrow, UT->nzmax, 0, 1, 0, CHOLMOD_REAL, &c);
    cholmod_transpose_unsym(UT, 1, NULL, NULL, 0, U, &c);
 
    cholmod_free_sparse(&UT, &c);
}

// Destroyer
container::~container()
{
    ;
}












