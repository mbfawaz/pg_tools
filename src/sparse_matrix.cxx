#include "sparse_matrix.h"
#include "umfpack.h"
#include "stdio.h"
#include "stdlib.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include "cholmod.h"
#include "amd.h"
#include "assert.h"
#include <cmath>

using namespace std;
using namespace matrix;

int sparse_matrix::_global_sparse_matrix_iterator = 0;

/* Namespace matrix functions */

//! This method faztorizes a sparse matrix A using umfpack, and stores the result in the void pointer Numeric
//! @param A [in]: Matrix to factorize
//! \return bool
bool matrix::lu_factorize(sparse_matrix &A)
{
    // Get matrix CCS vectors and params
    int*    Ai  = &A.get_Ai().at(0);   
    int*    Ap  = &A.get_Ap().at(0);   
    double* Ax  = &A.get_Ax().at(0);
    int nrows = A.get_number_of_rows();
    int ncols = A.get_number_of_columns();

    // Symbolic factorization
    void* Symbolic;
    int symbolic_status = umfpack_di_symbolic (nrows, ncols, Ap, Ai, Ax, &Symbolic, NULL, NULL);
    switch (symbolic_status)
    {
        case UMFPACK_OK:                     /*cerr<<"Symbolic factorization successful"<<endl;*/              break;
        case UMFPACK_ERROR_n_nonpositive:    cerr<<"Symbolic factorization error: n nonpositive"<<endl;    exit(1);
        case UMFPACK_ERROR_invalid_matrix:   cerr<<"Symbolic factorization error: invalid matrix"<<endl;   exit(1);
        case UMFPACK_ERROR_out_of_memory:    cerr<<"Symbolic factorization error: out of memory"<<endl;    exit(1);
        case UMFPACK_ERROR_argument_missing: cerr<<"Symbolic factorization error: argument missing"<<endl; exit(1);
        case UMFPACK_ERROR_internal_error:   cerr<<"Symbolic factorization: UMFPACK BUG"<<endl;            exit(1);
        default:                             cerr << "Unknown symbolic factorization error"<<endl;         exit(1);
    }

    // Numeric factorization
    void * Numeric;
    int numeric_status = umfpack_di_numeric (Ap, Ai, Ax, Symbolic, &Numeric, NULL, NULL);
    switch (numeric_status)
    {
      case UMFPACK_OK:                            /*cerr<<"Numeric factorization successful"<<endl;*/
                                                  A.set_lu_factorized(true); break;
      case UMFPACK_WARNING_singular_matrix:       cerr<<"Numeric factorization error: singular matrix"<<endl;                   exit(1);
      case UMFPACK_ERROR_out_of_memory:           cerr<<"Numeric factorization error: out of memory"<<endl;                     exit(1);
      case UMFPACK_ERROR_argument_missing:        cerr<<"Numeric factorization error: argument missing"<<endl;                  exit(1);
      case UMFPACK_ERROR_invalid_Symbolic_object: cerr<<"Numeric factorization error: invalid symbolic object"<<endl;           exit(1);
      case UMFPACK_ERROR_different_pattern:       cerr<<"Numeric factorization error: change in symbolic object pattern"<<endl; exit(1);
      default:                                    cerr<<"Unknown numeric factorization error"<<endl;                            exit(1);
    }    
    A.set_lu_numeric(Numeric);

    // Free memory
    umfpack_di_free_symbolic(&Symbolic);

    return true;
}

//! This method solves a system Ax = b using backward/forward substitution
//! @param A [in]: System representaion 
//! @param b [in]: right hand side of the system
//! \return vector<double>
vector<double> matrix::backward_forward_subtitution(sparse_matrix &A, vector<double> &b)
{
    // Check if lu_numeric pointer is assigned properly
    if (A.get_lu_numeric() == NULL)
    {
        cout<<"Error: Numeric factorization not available."<<endl;
        exit(0);
    }

    // Get CCS vectors and parameters
    int*    Ai = A.get_Ai().empty()? (int*)   malloc(0): &A.get_Ai().at(0); 
    int*    Ap = A.get_Ap().empty()? (int*)   malloc(0): &A.get_Ap().at(0); 
    double* Ax = A.get_Ai().empty()? (double*)malloc(0): &A.get_Ax().at(0);
    double* bb = &b.at(0); 
   
    // Pointer to the solution
    double* xx  = (double*) malloc (b.size()*sizeof(double));
    
    // UMFPACK backward/forward solve
    int solve_status = umfpack_di_solve(UMFPACK_A, Ap, Ai, Ax, xx, bb, A.get_lu_numeric(), (double*) NULL, (double*) NULL);
    if (solve_status != UMFPACK_OK) 
    { 
        cerr<<"Error in matrix solution";
        exit(1);
    }

    // Vector solution
    vector<double> solution(b.size());
    for (int unsigned i = 0; i<b.size(); i++)
        solution[i] = xx[i];

    free(xx);

    // Return the solution
    return solution;
}

//! This method faztorizes a sparse matrix A using cholmod, and stores the result in the pointer _chol_L
//! @param A [in]: Matrix to factorize
void matrix::chol_factorize(sparse_matrix &A)
{
    // Get the cholmod representation of A
    cholmod_common c;
    cholmod_start (&c);
    cholmod_factor * chol_L;

    cholmod_sparse * chol_A;
    vector<int>      Ai;
    vector<int>      Ap;
    vector<double>   Ax;

    chol_A = cholmod_allocate_sparse(A.get_number_of_rows(), A.get_number_of_columns(), A.get_nz(), 0, 1, 0, CHOLMOD_REAL, &c);
    int    *chol_Ai = (int*)   chol_A->i;
    int    *chol_Ap = (int*)   chol_A->p;
    double *chol_Ax = (double*)chol_A->x;
    chol_A->stype = 1; // A matrix is symetric positive definite SPD

    Ai = A.get_Ai(); Ap = A.get_Ap(); Ax = A.get_Ax();
    for(int jj=0;jj<A.get_nz();jj++)
    {
        chol_Ai[jj] = Ai[jj];
        chol_Ax[jj] = Ax[jj];
    }
    for(int jj=0;jj<A.get_number_of_columns()+1;jj++)
        chol_Ap[jj] = Ap[jj];

    chol_L = cholmod_analyze (chol_A, &c);
    cholmod_factorize (chol_A, chol_L, &c);
    A.set_chol_L(chol_L);

    cholmod_free_sparse (&chol_A, &c) ;
    cholmod_finish (&c) ;
}

//! This method solves a system Ax = b using cholmod
//! @param A [in]: System representaion 
//! @param b [in]: right hand side of the system
//! \return vector<double>
vector<double> matrix::chol_solve(sparse_matrix &A, vector<double> &b)
{
    //return vector
    vector<double> rv(A.get_number_of_rows());

    // Number of equations
    int nrow = A.get_number_of_rows();

    // Start cholmod environment
    cholmod_common c;
    cholmod_start (&c);
    
    // will contain the solution
    cholmod_dense *x;

    // contains the RHS
    cholmod_dense* chol_b = cholmod_ones(nrow, 1, CHOLMOD_REAL, &c);
    double * chol_bx = (double*) chol_b->x;
    for (unsigned int i = 0; i<b.size(); i++)
        chol_bx[i] = b[i];
    
    x = cholmod_solve (CHOLMOD_A, A.get_chol_L(), chol_b, &c) ;

    double * chol_xx = (double*) x->x;
    for (int i = 0; i<A.get_number_of_rows(); i++)
        rv[i] = chol_xx[i];

    cholmod_free_dense (&x, &c) ;
    cholmod_free_dense (&chol_b, &c) ;
    cholmod_finish (&c);

    return rv;
}

/* sparse matrix addition routine */
/* C = alpha*A + beta*B */
//! This method performs addition of two sparse matrices such that C = alpha*A + beta*B
//! @param alpha [in]: Scalar to be multiplied by the first matrix
//! @param A     [in]: First matrix
//! @param beta  [in]: Scalar to be multiplied by the second matrix
//! @param alpha [in]: Second matrix
//! @param C     [out]: Resulting sparse matrix
void matrix::sparse_matrix_add(double alpha, sparse_matrix & A, 
                               double beta , sparse_matrix & B,
                                             sparse_matrix & C)
{
    /* check if the matrices have the right dimensions */
    if(   ( A.get_number_of_columns() != B.get_number_of_columns())  
       && ( A.get_number_of_rows()    != B.get_number_of_rows()   ))
	    printf("Matrix dimensions do not agree.\n");
    
    else
    {
        /* Set the dimensions of C */
        C.set_number_of_rows   (A.get_number_of_rows());
        C.set_number_of_columns(A.get_number_of_columns());

        /* Cholmod variables */
        cholmod_common Common, *cc;

        /* Matrix A */
        cholmod_sparse * chol_A;
        int        * chol_Ai, * chol_Ap;
        double         * chol_Ax;
        vector<int>      Ai;
        vector<int>      Ap;
        vector<double>   Ax;

        /* Matrix B */
        cholmod_sparse * chol_B;
        int        * chol_Bi, * chol_Bp;
        double         * chol_Bx;
        vector<int>      Bi;
        vector<int>      Bp;
        vector<double>   Bx;

        /* the result matrix C */
        cholmod_sparse *chol_C;
        int        * chol_Ci, * chol_Cp;
        double         * chol_Cx;
        vector<int>      Ci;
        vector<int>      Cp;
        vector<double>   Cx;

        /* start cholmod */
        cc = &Common; cholmod_l_start(cc);

        /* create the cholmod matrix A */
        chol_A = cholmod_l_allocate_sparse(A.get_number_of_rows(), A.get_number_of_columns(),A.get_nz(), 0,1,0,CHOLMOD_REAL,cc);
        chol_Ai = (int*)chol_A->i;
        chol_Ap = (int*)chol_A->p;
        chol_Ax = (double*) chol_A->x;
        chol_A->nzmax = A.get_nz();

        Ai = A.get_Ai(); Ap = A.get_Ap(); Ax = A.get_Ax();
        for(int jj=0;jj<A.get_nz();jj++)
	    {
	        chol_Ai[jj] = Ai[jj];
	        chol_Ax[jj] = Ax[jj];
	    }
	    for(int jj=0;jj<A.get_number_of_columns()+1;jj++)
	        chol_Ap[jj] = Ap[jj];

	    /* create the cholmod matrix B */
        chol_B = cholmod_l_allocate_sparse(B.get_number_of_rows(), B.get_number_of_columns(),B.get_nz(), 0,1,0,CHOLMOD_REAL,cc);
        chol_Bi = (int*)chol_B->i;
        chol_Bp = (int*)chol_B->p;
        chol_Bx = (double*) chol_B->x;
        chol_B->nzmax = B.get_nz();

        Bi = B.get_Ai(); Bp = B.get_Ap(); Bx = B.get_Ax();
        for(int jj=0;jj<B.get_nz();jj++)
        {
            chol_Bi[jj] = Bi[jj];
            chol_Bx[jj] = Bx[jj];
        }

        for(int jj=0;jj<B.get_number_of_columns()+1;jj++)
            chol_Bp[jj] = Bp[jj];

        /* Scaling factors */
        double Alpha[2] = {0, 0}; Alpha[0] = alpha;
        double Beta[2]  = {0, 0}; Beta[0] = beta;

        /* perform the addition */
        chol_C = cholmod_l_add(chol_A,chol_B,Alpha,Beta,1,1,cc); 

        /* Copy the result to our matrix class */
        chol_Ci = (int*)chol_C->i;
        chol_Cp = (int*)chol_C->p;
        chol_Cx = (double*) chol_C->x;

        /* copy the compressed column format */
        Ci.resize(chol_C->nzmax);
        Cx.resize(chol_C->nzmax);
        Cp.resize(C.get_number_of_columns()+1);
        for(unsigned int jj=0; jj<chol_C->nzmax; jj++)
        {
            Ci[jj] = chol_Ci[jj];
            Cx[jj] = chol_Cx[jj];
        }

        for(int jj=0;jj<C.get_number_of_columns()+1;jj++)
	        Cp[jj] = chol_Cp[jj];

        C.set_Ai(Ci);
        C.set_Ap(Cp);
        C.set_Ax(Cx);
        C.set_nz(chol_C->nzmax);

        /* free the cholmod structures */
        cholmod_l_free_sparse(&chol_A,cc);
        cholmod_l_free_sparse(&chol_B,cc);
        cholmod_l_free_sparse(&chol_C,cc);
        cholmod_l_finish(cc);
    }
}

/* sparse matrix multiplication routine */
/* C = A*B */
//! This method performs multiplication of two sparse matrices 
//! @param A [in]: First matrix
//! @param B [in]: Second matrix
//! @param C [out]: Resulting sparse matrix
void matrix::sparse_matrix_mult(sparse_matrix & A, 
                                sparse_matrix & B,
                                sparse_matrix & C)
{
    /* check if the matrices have the right dimensions */
    if(A.get_number_of_columns() != B.get_number_of_rows())  
	    printf("Matrix dimensions do not agree.\n");

    else
    {
        /* Set the dimensions of C */
        C.set_number_of_rows   (A.get_number_of_rows());
        C.set_number_of_columns(B.get_number_of_columns());

        /* Cholmod variables */
        cholmod_common Common, *cc;

        /* Matrix A */
        cholmod_sparse * chol_A;
        int        * chol_Ai, * chol_Ap;
        double         * chol_Ax;
        vector<int>      Ai;
        vector<int>      Ap;
        vector<double>   Ax;

        /* Matrix B */
        cholmod_sparse * chol_B;
        int        * chol_Bi, * chol_Bp;
        double         * chol_Bx;
        vector<int>      Bi;
        vector<int>      Bp;
        vector<double>   Bx;

        /* the result matrix C */
        cholmod_sparse * chol_C;
        int        * chol_Ci, * chol_Cp;
        double         * chol_Cx;
        vector<int>      Ci;
        vector<int>      Cp;
        vector<double>   Cx;

        /* start cholmod */
        cc = &Common; cholmod_l_start(cc);

        /* create the cholmod matrix A */
        chol_A = cholmod_l_allocate_sparse(A.get_number_of_rows(), A.get_number_of_columns(),A.get_nz(), 0,1,0,CHOLMOD_REAL,cc);
        chol_Ai = (int*)chol_A->i;
        chol_Ap = (int*)chol_A->p;
        chol_Ax = (double*) chol_A->x;
        chol_A->nzmax = A.get_nz();

        Ai = A.get_Ai(); Ap = A.get_Ap(); Ax = A.get_Ax();
        for(int jj=0;jj<A.get_nz();jj++)
	    {
	        chol_Ai[jj] = Ai[jj];
	        chol_Ax[jj] = Ax[jj];
	    }
	    for(int jj=0;jj<A.get_number_of_columns()+1;jj++)
	        chol_Ap[jj] = Ap[jj];

	    /* create the cholmod matrix B */
        chol_B = cholmod_l_allocate_sparse(B.get_number_of_rows(), B.get_number_of_columns(),B.get_nz(), 0,1,0,CHOLMOD_REAL,cc);
        chol_Bi = (int*)chol_B->i;
        chol_Bp = (int*)chol_B->p;
        chol_Bx = (double*) chol_B->x;
        chol_B->nzmax = B.get_nz();

        Bi = B.get_Ai(); Bp = B.get_Ap(); Bx = B.get_Ax();
        for(int jj=0;jj<B.get_nz();jj++)
        {
            chol_Bi[jj] = Bi[jj];
            chol_Bx[jj] = Bx[jj];
        }

        for(int jj=0;jj<B.get_number_of_columns()+1;jj++)
            chol_Bp[jj] = Bp[jj];

        /* perform the multiplication */
        chol_C = cholmod_l_ssmult(chol_A,chol_B,0,1,1,cc); 

        /* Copy the result to our matrix class */
        /* [ There is a bug in cholmod when the result is 0. In this case, it sees
         *   that nz = 1 and the entry is zero and the column/row index in invalid ]*/
        chol_Ci = (int*)chol_C->i;
        chol_Cp = (int*)chol_C->p;
        chol_Cx = (double*) chol_C->x;
       
        if ((chol_Cx[0] == 0) && (chol_C->nzmax == 1) )
        {
            // Result is zero. Just set nz = 0
	        C.set_nz(0);
	    }
	    else
	    {
            /* copy the compressed column format as usual */
            Ci.resize(chol_C->nzmax);
            Cx.resize(chol_C->nzmax);
            Cp.resize(C.get_number_of_columns()+1);
            for(unsigned int jj=0;jj<chol_C->nzmax;jj++)
            {
                Ci[jj] = chol_Ci[jj];
                Cx[jj] = chol_Cx[jj];
            }

            for(int jj=0;jj<C.get_number_of_columns()+1;jj++)
                Cp[jj] = chol_Cp[jj];

            C.set_Ai(Ci);
            C.set_Ap(Cp);
            C.set_Ax(Cx);
            C.set_nz(chol_C->nzmax);
        }

        /* free the cholmod structures */
        cholmod_l_free_sparse(&chol_A,cc);
        cholmod_l_free_sparse(&chol_B,cc);
        cholmod_l_free_sparse(&chol_C,cc);
        cholmod_l_finish(cc);
    }
}


/* sparse matrix-vecotr multiplication routines */
//! This method performs sparse matrix vector multiplication
//! @param A [in]: Sparse matrix to multiply 
//! @param b [in]: Vector to multiply
//! \return vector<double>
vector<double> matrix::sparse_matrix_vector_mult(sparse_matrix &A, vector<double> &b)
{
    /* Make sure the number of columns of A is equal to the size of b */
    assert ((unsigned)A.get_number_of_columns() == b.size());

    /* start the matrix-vector multiplication */
    vector<double> c(A.get_number_of_rows(),0);
    vector<int> Ap, Ai;
    vector<double> Ax;

    Ap = A.get_Ap(); Ai = A.get_Ai(); Ax = A.get_Ax();

    for(int j = 0; j < A.get_number_of_columns(); j++)
        for(int p = Ap[j]; p < Ap[j+1]; p++)
            c[Ai[p]] += Ax[p]*b[j];

    return c;
}

/* sparse vecotr-matrix multiplication routines */
//! This method performs vector sparse matrix multiplication
//! @param b [in]: Vector to multiply
//! @param A [in]: Sparse matrix to multiply 
//! \return vector<double>
vector<double> matrix::sparse_vector_matrix_mult(vector<double> &b, sparse_matrix &A)
{
    /* Make sure the number of rows of A is equal to the size of b */
    assert ((unsigned)A.get_number_of_rows() == b.size());

    /* start the vector-matrix multiplication */
    vector<double> c(A.get_number_of_columns(), 0);
    vector<int> Ap, Ai;
    vector<double> Ax;

    Ap = A.get_Ap(); Ai = A.get_Ai(); Ax = A.get_Ax();

    for (unsigned int i = 0; i<Ap.size()-1; i++)
        for (int j = Ap[i]; j<Ap[i+1]; j++)
            c[i] += b[Ai[j]] * Ax[j];     
    
    return c;
}

/* drop either the positive or the negative elements */
/* 1: drops the positive entries, -1: drops the negative entries */
//! This method drops either the positive of the negative elements of a matrix
//! @param A    [in]:  Sparse matrix to drop entries from
//! @param C    [out]: Resulting sparse matrix
//! @param sign [in]:  If 1, then positive. If -1, then negative.
void matrix::sparse_drop_entries(sparse_matrix &A, int sign, sparse_matrix &C)
{
    // Get the triplet format to make things easier
    bool status = A.get_triplet_format();
    
    // Number of non-zero entries in C
    int c_nz = 0;

    // Drop unwanted entries
    if((sign == 1 || sign == -1) && status)
    {
        // Get all the required entries
        for (unsigned int i = 0; i<A.get_Ti().size(); i++)
        {
            if (sign * (A.get_Tx()[i]) < 0)
            {
                C.get_Ti().push_back(A.get_Ti()[i]);
                C.get_Tj().push_back(A.get_Tj()[i]);
                C.get_Tx().push_back(A.get_Tx()[i]);
                c_nz++;
            }
        }

        // Update the parameters of C
        C.set_nz(c_nz);
        C.set_number_of_rows   (A.get_number_of_rows   ());
        C.set_number_of_columns(A.get_number_of_columns());
       
        // Get the column format of C
        C.get_column_format();

        // Delete the triplet format
        C.get_Ti().clear();
        C.get_Tj().clear();
        C.get_Tx().clear();
    }
    else
    	cerr<<"Failed to drop entries."<<endl;
}

/* transpose matrix A and put the result in C */
//! This method performs a transpose of a sparse matrix
//! @param A [in]: Sparse matrix to transpose
//! @param B [in]: Resulting sparse matrix
void matrix::transpose(sparse_matrix & A, 
                       sparse_matrix & C)
{
    /* Set the dimensions of C */
    C.set_number_of_rows   (A.get_number_of_columns());
    C.set_number_of_columns(A.get_number_of_rows()   );

    /* Cholmod variables */
    cholmod_common Common, *cc;

    /* Matrix A */
    cholmod_sparse * chol_A;
    int        * chol_Ai, * chol_Ap;
    double         * chol_Ax;
    vector<int>      Ai;
    vector<int>      Ap;
    vector<double>   Ax;

    /* the result matrix C */
    cholmod_sparse * chol_C;
    int        * chol_Ci, * chol_Cp;
    double         * chol_Cx;
    vector<int>      Ci;
    vector<int>      Cp;
    vector<double>   Cx;

    /* start cholmod */
    cc = &Common; cholmod_l_start(cc);

    /* create the cholmod matrix A */
    chol_A = cholmod_l_allocate_sparse(A.get_number_of_rows(), A.get_number_of_columns(),A.get_nz(), 0,1,0,CHOLMOD_REAL,cc);
    chol_Ai = (int*)chol_A->i;
    chol_Ap = (int*)chol_A->p;
    chol_Ax = (double*) chol_A->x;
    chol_A->nzmax = A.get_nz();

    Ai = A.get_Ai(); Ap = A.get_Ap(); Ax = A.get_Ax();
    for(int jj=0;jj<A.get_nz();jj++)
    {
        chol_Ai[jj] = Ai[jj];
        chol_Ax[jj] = Ax[jj];
    }
    for(int jj=0;jj<A.get_number_of_columns()+1;jj++)
        chol_Ap[jj] = Ap[jj];

    /* perform the transpose */
    chol_C = cholmod_l_allocate_sparse(A.get_number_of_columns(), A.get_number_of_rows(), A.get_nz(), 0,1,0,CHOLMOD_REAL,cc);
    cholmod_l_transpose_unsym(chol_A, 1, NULL, NULL, 0 ,chol_C, cc); 

    /* Copy the result to our matrix class */
    chol_Ci = (int*)chol_C->i;
    chol_Cp = (int*)chol_C->p;
    chol_Cx = (double*) chol_C->x;
    
    /* copy the compressed column format as usual */
    Ci.resize(chol_C->nzmax);
    Cx.resize(chol_C->nzmax);
    Cp.resize(C.get_number_of_columns()+1);
    for(unsigned int jj=0;jj<chol_C->nzmax;jj++)
    {
        Ci[jj] = chol_Ci[jj];
        Cx[jj] = chol_Cx[jj];
    }

    for(int jj=0;jj<C.get_number_of_columns()+1;jj++)
        Cp[jj] = chol_Cp[jj];

    C.set_Ai(Ci);
    C.set_Ap(Cp);
    C.set_Ax(Cx);
    C.set_nz(chol_C->nzmax);

    /* free the cholmod structures */
    cholmod_l_free_sparse(&chol_A,cc);
    cholmod_l_free_sparse(&chol_C,cc);
    cholmod_l_finish(cc);
}

/* Horizontally augment two matrices */
/* C = [A B] */
//! This method horizontally augment two matrices 
//! @param A [in]: First matrix
//! @param B [in]: Second matrix
//! @param C [out]: Resulting sparse matrix 
void matrix::sparse_matrix_horizontal_augment(sparse_matrix &A, sparse_matrix &B, sparse_matrix &C)
{
    /* Check if the matrices have the same number of rows */
    if( A.get_number_of_rows() != B.get_number_of_rows() )
	    printf("Matrices should have the same number of rows.\n");
    else
    {
	    /* set the appropriate dimensions of C, and update the non-zeros*/
	    C.set_number_of_rows   (A.get_number_of_rows());
	    C.set_number_of_columns(A.get_number_of_columns() + B.get_number_of_columns());
	    C.set_nz(A.get_nz() + B.get_nz() );

	    /* start the augmentation process */
	    /* this is done in triplet format, so both matrices should be in triplet
	    * format */
        if( A.get_triplet_format() && B.get_triplet_format() )
        {
            /* make sure the triplet arrays have the same size */
            assert( A.get_Ti().size() == A.get_Tj().size() );
            assert( A.get_Ti().size() == A.get_Tx().size() );

            /* make sure the triplet arrays have the same size */
            assert (B.get_Ti().size() == B.get_Tj().size() );
            assert (B.get_Ti().size() == B.get_Tx().size() );

            C.get_Ti().resize(A.get_Ti().size() + B.get_Ti().size());
            C.get_Tj().resize(A.get_Tj().size() + B.get_Tj().size());
            C.get_Tx().resize(A.get_Tx().size() + B.get_Tx().size());

            for (unsigned int i = 0; i<A.get_Ti().size(); i++)
            {
                /* entries of A are copied in the same positions */
                C.get_Ti()[i] = A.get_Ti()[i];
                C.get_Tj()[i] = A.get_Tj()[i];
                C.get_Tx()[i] = A.get_Tx()[i];
            }

            for (unsigned int i = 0; i<B.get_Ti().size(); i++)
            {
                /* no change to the row index of the right matrix */
                C.get_Ti()[A.get_Ti().size() + i] = B.get_Ti()[i];

                /* the column index in shifted by the number of columns of the left
                * matrix */
                C.get_Tj()[A.get_Tj().size() + i] = A.get_number_of_columns() + B.get_Tj()[i];

                /* push back the value */
                C.get_Tx()[A.get_Tx().size() + i] = B.get_Tx()[i];
            }
            
            /* get the column format of the matrix */
            C.get_column_format();
        }
        else
            printf("Horizantal Augmentation Failed.\n");
    }
}

/* Vertically augment two matrices */
/* C = [A ; B] */
//! This method vertically augment two matrices 
//! @param A [in]: First matrix
//! @param B [in]: Second matrix
//! @param C [out]: Resulting sparse matrix
void matrix::sparse_matrix_vertical_augment(sparse_matrix &A, sparse_matrix &B, sparse_matrix &C)
{
    /* Check of the matrices have the same number of columns */
    if( A.get_number_of_columns() != B.get_number_of_columns() )
	    printf("Matrices should have the same number of columns.\n");
    else
    {
	    /* set the appropriate dimensions of C, and update the non-zeros*/
	    C.set_number_of_rows   (A.get_number_of_rows() + B.get_number_of_rows());
	    C.set_number_of_columns(A.get_number_of_columns());
	    C.set_nz(A.get_nz() + B.get_nz() );

	    /* start the augmentation process */
	    /* this is done in triplet format, so both matrices should be in triplet
	    * format */
        if( A.get_triplet_format() && B.get_triplet_format() )
        {
            /* make sure the triplet arrays have the same size */
            assert( A.get_Ti().size() == A.get_Tj().size() );
            assert( A.get_Ti().size() == A.get_Tx().size() );

            /* make sure the triplet arrays have the same size */
            assert (B.get_Ti().size() == B.get_Tj().size() );
            assert (B.get_Ti().size() == B.get_Tx().size() );

            C.get_Ti().resize(A.get_Ti().size() + B.get_Ti().size());
            C.get_Tj().resize(A.get_Tj().size() + B.get_Tj().size());
            C.get_Tx().resize(A.get_Tx().size() + B.get_Tx().size());

            for (unsigned int i = 0; i<A.get_Ti().size(); i++)
            {
                /* entries of A are copied in the same positions */
                C.get_Ti()[i] = A.get_Ti()[i];
                C.get_Tj()[i] = A.get_Tj()[i];
                C.get_Tx()[i] = A.get_Tx()[i];
            }

            for (unsigned int i = 0; i<B.get_Ti().size(); i++)
            {
                /* no change to the column index of the right matrix */
                C.get_Tj()[A.get_Tj().size() + i] = B.get_Tj()[i];

                /* the row index in shifted by the number of rows of the left
                * matrix */
                C.get_Ti()[A.get_Ti().size() + i] = A.get_number_of_rows() + B.get_Ti()[i];

                /* push back the value */
                C.get_Tx()[A.get_Tx().size() + i] = B.get_Tx()[i];
            }

            /* get the column format of the matrix */
            C.get_column_format();
        }
        else
            printf("Horizantal Augmentation Failed.\n");
    }
}

// vector-vector operations
vector<double> &matrix::operator+(vector<double> &v1, vector<double> &v2)
{
    // Make sure the sizes are the same
    assert(v1.size() == v2.size());

    vector<double> *rv = new vector<double>(v1.size());
    for (unsigned int i = 0; i<v1.size(); i++)
        rv->at(i) = v1[i] + v2[i];

    return *rv;
}

vector<double> &matrix::operator-(vector<double> &v1, vector<double> &v2)
{
    // Make sure the sizes are the same
    assert(v1.size() == v2.size());

    vector<double> *rv = new vector<double>(v1.size());
    for (unsigned int i = 0; i<v1.size(); i++)
        rv->at(i) = v1[i] - v2[i];

    return *rv;
}

double matrix::operator*(vector<double> &v1, vector<double> &v2)
{
    // Make sure the sizes are the same
    assert(v1.size() == v2.size());

    double rv=0;
    for (unsigned int i = 0; i<v1.size(); i++)
        rv += v1[i] * v2[i];

    return rv;
}

vector<double> &matrix::operator*(double alpha, vector<double> &v1)
{
    vector<double> *rv = new vector<double>(v1.size());
    for (unsigned int i = 0; i<v1.size(); i++)
        rv->at(i) = alpha * v1[i];

    return *rv;
}

vector<double> &matrix::operator*(vector<double> &v1, double alpha )
{
    vector<double> *rv = new vector<double>(v1.size());
    for (unsigned int i = 0; i<v1.size(); i++)
        rv->at(i) = alpha * v1[i];

    return *rv;
}

/* Overloading matrix-vector operators */
vector<double> &matrix::operator|(sparse_matrix &A, vector<double> &b)
{
    // Check if no factorization exists 
    assert(A.get_lu_numeric() != NULL || A.get_chol_L() != NULL);
       
    vector<double> *rv = new vector<double>();

    // If Chol factorization exists
    if (A.get_chol_L() != NULL)
    {
        *rv = chol_solve(A, b);
        return *rv;
    }

    // If LU factorization exists
    if (A.get_lu_numeric() != NULL)
    {
        *rv = backward_forward_subtitution(A, b);
        return *rv;
    }

    // Default
    return *rv;
}


/* Overloading matrix-matrix operator */
vector<double> &matrix::operator*(sparse_matrix &A, vector<double> &b)
{
    vector<double> *rv = new vector<double>();
    *rv = sparse_matrix_vector_mult(A, b);
    return *rv;
}

sparse_matrix &matrix::operator+(sparse_matrix &A, sparse_matrix &B)
{
    sparse_matrix *C = new sparse_matrix();;
    sparse_matrix_add(1, A, 1, B, *C);
    return *C;
}

sparse_matrix &matrix::operator-(sparse_matrix &A, sparse_matrix &B)
{
    sparse_matrix *C = new sparse_matrix();
    sparse_matrix_add(1, A, -1, B, *C);
    return *C;
}

sparse_matrix &matrix::operator*(sparse_matrix &A, sparse_matrix &B)
{
    sparse_matrix *C = new sparse_matrix();
    sparse_matrix_mult(A, B, *C);
    return *C;
}

sparse_matrix &matrix::operator*(double alpha, sparse_matrix &A)
{
    sparse_matrix *C = new sparse_matrix(A);
    for (unsigned int i = 0; i<C->get_Ax().size(); i++)
        C->get_Ax()[i] *= alpha;
    return *C;
}

sparse_matrix &matrix::operator*(sparse_matrix &A, double alpha)
{
    sparse_matrix *C = new sparse_matrix(A);
    for (unsigned int i = 0; i<C->get_Ax().size(); i++)
        C->get_Ax()[i] *= alpha;
    return *C;
}

sparse_matrix &matrix::operator||(sparse_matrix &A, sparse_matrix &B)
{
    sparse_matrix *C = new sparse_matrix();
    sparse_matrix_horizontal_augment(A, B, *C);
    return *C;
}

sparse_matrix &matrix::operator/(sparse_matrix &A, sparse_matrix &B)
{
    sparse_matrix *C = new sparse_matrix();
    sparse_matrix_vertical_augment(A, B, *C);
    return *C;
}

sparse_matrix &matrix::operator~(sparse_matrix &A)
{
    sparse_matrix *C = new sparse_matrix();
    transpose(A, *C);
    return *C;
}


/* Factorize */
void matrix::factorize(sparse_matrix &A, int method)
{
    if (method == 0)
        lu_factorize(A);
    else if (method == 1)
        chol_factorize(A);
    else
    {
        cout<<"Error. Factorization method not supported."<<endl;
        exit(1);
    }
}

sparse_matrix &matrix::inv(sparse_matrix &A, vector<int> & cols, double delta)
{
    assert(A.get_number_of_rows() > 0);

    sparse_matrix *Ainv = new sparse_matrix();
    if (cols.size() == 0)
        return *Ainv;
    int n = A.get_number_of_rows();
    int m = cols.size();

    Ainv->set_number_of_rows(n);
    Ainv->set_number_of_columns(m);
    vector<double> ei(n,0);
    vector<double> col_of_inv;
    long nz = 0;

    Ainv->get_Ap().push_back(0);

    for (int i = 0; i<m; i++)
    {
        ei[i] = 1;
        
        col_of_inv = A|ei;
        for (int j = 0; j<n; j++)
        {
            if (abs(col_of_inv[j]) >= delta)
            {
                Ainv->get_Ai().push_back(j);
                Ainv->get_Ax().push_back(col_of_inv[j]);
                nz++;
            }
        }
        Ainv->get_Ap().push_back(nz);

        ei[i] = 0;

        //Track progress
        if (m > 100)
        {
            if (i % (int)(m/100) == 0)
            {
                printf("\r%2d%% complete...", (int)(100.0*i/m));
                fflush(stdout);
                printf("\r");
            }
        }
    }
    Ainv->set_nz(nz);

    return *Ainv;
}

double matrix::average(vector<double> &v)
{
    double rv = 0;
    for (unsigned int i = 0; i<v.size(); i++)
        rv += v[i];

    return rv/v.size();
}

double matrix::norm(vector<double> &v, int p) //support p=1,2,INF
{
    if (v.size() == 0)
        return 0;

    double rv = 0;
    
    if (p == 1)
    {
        for (unsigned int i = 0; i<v.size(); i++)
		    rv += abs(v[i]);
	    return rv;
    }
    else if (p == 2)
    {
        for (unsigned int i = 0; i<v.size(); i++)
		    rv += v[i]*v[i];
	    return sqrt(rv);
    }
    else if (p == INF)
    {
        for (unsigned int i = 0; i<v.size(); i++)
            if (abs(v[i]) > rv)
                rv = abs(v[i]);
        return rv;
    }
    else
    {
        cout<<"Error. p-norm requested is not supported."<<endl;
        exit(1);
    }
}

vector<double> matrix::randv(int size)
{
    vector<double> rv(size);
    for (int i = 0; i<size; i++)
        rv[i] = 1.0*rand()/RAND_MAX;
    return rv;
}

void matrix::print_matrix(sparse_matrix &A)
{

    cout<<"Printing matrix..."<<endl;
    for (int i = 0; i<A.get_number_of_rows(); i++)
    {
        for (int j = 0; j<A.get_number_of_columns(); j++)
            printf("%6e ", A.get_element(i,j));
        printf("\n");
    }
}

void matrix::print_vector(vector<double> &v)
{
    for (unsigned int i = 0; i<v.size(); i++)
        cout<<v[i]<<" ";
    cout<<endl;
}

/* Class sparse_matrix functions */
//! Default Constructor
sparse_matrix::sparse_matrix()
{
    _global_sparse_matrix_iterator++;
    _local_sparse_matrix_iterator = _global_sparse_matrix_iterator;
}

//! Regular Constructors
//! @param number_of_rows    [in]: Number of rows of the sparse matrix
//! @param number_of_columns [in]: Number of columns of the sparse matrix
sparse_matrix::sparse_matrix(int number_of_rows, int number_of_columns)
{
    _number_of_rows     = number_of_rows;
    _number_of_columns  = number_of_columns;
    _global_sparse_matrix_iterator++;
    _local_sparse_matrix_iterator = _global_sparse_matrix_iterator;
}

//! Regular Constructors
//! @param number_of_rows    [in]: Number of rows of the sparse matrix
//! @param number_of_columns [in]: Number of columns of the sparse matrix 
//! @param nz                [in]: Number of non-zeros in the sparse matrix
sparse_matrix::sparse_matrix(int number_of_rows, int number_of_columns, int nz)
{
    _number_of_rows     = number_of_rows;
    _number_of_columns  = number_of_columns;
    _nz                 = nz;
    _global_sparse_matrix_iterator++;
    _local_sparse_matrix_iterator = _global_sparse_matrix_iterator;
}

//! Copy Constructor
//! @param sparse_matrix_to_copy [in]: Sparse matrix to copy
sparse_matrix::sparse_matrix(sparse_matrix & sparse_matrix_to_copy)
{
    _number_of_rows = sparse_matrix_to_copy.get_number_of_rows();
    _number_of_columns = sparse_matrix_to_copy.get_number_of_columns();
    _nz = sparse_matrix_to_copy.get_nz();

    _Ai = sparse_matrix_to_copy.get_Ai();
    _Ap = sparse_matrix_to_copy.get_Ap();
    _Ax = sparse_matrix_to_copy.get_Ax();

    _Ti = sparse_matrix_to_copy.get_Ti();
    _Tj = sparse_matrix_to_copy.get_Tj();
    _Tx = sparse_matrix_to_copy.get_Tx();

    _lu_factorized = sparse_matrix_to_copy.is_lu_factorized();

    _global_sparse_matrix_iterator++;
    _local_sparse_matrix_iterator = _global_sparse_matrix_iterator;

    _chol_L = sparse_matrix_to_copy.get_chol_L();
    _lu_numeric = sparse_matrix_to_copy.get_lu_numeric();
}


// Overload Operators
//! This method is an overload of the assignment operator
//! @param M [in]: Sparse matrix to assign to the new sparse matrix
//! \return sparse_matrix
sparse_matrix & sparse_matrix::operator= (sparse_matrix & sparse_matrix_to_copy)
{
    if (this == &sparse_matrix_to_copy)
        return (*this);
 
    _number_of_rows = sparse_matrix_to_copy.get_number_of_rows();
    _number_of_columns = sparse_matrix_to_copy.get_number_of_columns();
    _nz = sparse_matrix_to_copy.get_nz();

    _Ai = sparse_matrix_to_copy.get_Ai();
    _Ap = sparse_matrix_to_copy.get_Ap();
    _Ax = sparse_matrix_to_copy.get_Ax();

    _Ti = sparse_matrix_to_copy.get_Ti();
    _Tj = sparse_matrix_to_copy.get_Tj();
    _Tx = sparse_matrix_to_copy.get_Tx();

    _lu_factorized = sparse_matrix_to_copy.is_lu_factorized();

    _global_sparse_matrix_iterator++;
    _local_sparse_matrix_iterator = _global_sparse_matrix_iterator;
    
    return *this;
}

//! Destroyer
sparse_matrix::~sparse_matrix()
{
    ; 
}

// Getter Functions
//! This method returns the Ap vector of a sparse matrix in compressed column format
//! \retunr vector<int>
vector<int> & sparse_matrix::get_Ap()
{
    return _Ap;
}

//! This method returns the Ai vector of a sparse matrix in compressed column format
//! \retunr vector<int>
vector<int> & sparse_matrix::get_Ai()
{
    return _Ai;
}

//! This method returns the Ax vector of a sparse matrix in compressed column format
//! \retunr vector<double>
vector<double> & sparse_matrix::get_Ax()
{
    return _Ax;
}

//! This method returns the Ti vector of a sparse matrix in triplet format
//! \retunr vector<int>
vector<int> & sparse_matrix::get_Ti()
{
    return _Ti;
}

//! This method returns the Tj vector of a sparse matrix in triplet format
//! \retunr vector<int>
vector<int> & sparse_matrix::get_Tj()
{
    return _Tj;
}

//! This method returns the Tx vector of a sparse matrix in triplet format
//! \retunr vector<double>
vector<double> & sparse_matrix::get_Tx()
{
    return _Tx;
}

//! This method returns the number of rows of the sparse matrix
//! \return int
int sparse_matrix::get_number_of_rows()
{
    return _number_of_rows;
}

//! This method returns the number of columns of the sparse matrix
//! \return int
int sparse_matrix::get_number_of_columns ()
{
    return _number_of_columns;
}

//! This method returns the number of nonzeros of the sparse matrix
//! \return int
int sparse_matrix::get_nz()
{
    return _nz;
}

//! This method returns global index of the sparse matrix
//! \return int
int sparse_matrix::get_global_sparse_matrix_iterator()
{
    return _global_sparse_matrix_iterator;
}

//! This method returns the local index of the sparse matrix
//! \return int
int sparse_matrix::get_local_sparse_matrix_iterator()
{
    return _local_sparse_matrix_iterator;
}

//! This method builds the triplet format of the sparse matrix
//! \return bool
bool sparse_matrix::get_triplet_format()
{
    if (_nz == 0)
        return true;
    
    // Allocate memory for Tj
    int *Tj  = (int*) malloc ((_nz) * sizeof(int));

    // Get the array pointer from the vector _Ap
    int *cAp = &_Ap[0];
  
    // Get cAp from Tj
    int status = umfpack_di_col_to_triplet(get_number_of_columns(),cAp,Tj);

    // Check for errors
    if (status < 0)
    {
	    cerr<<"Failed column-to-triplet conversion"<<std::endl;
	    return false;
    }
    else
    {
        // Copy result into triplet format
	    _Ti = _Ai;
	    _Tx = _Ax;
    
        _Tj.resize(_nz);
	    for(int i =0; i< _nz;i++)
	        _Tj[i] = Tj[i];
    }

    free(Tj);
    return true;
}

//! This method builds the compressed column format of the sparse matrix
//! \return bool
bool sparse_matrix::get_column_format()
{
    if (_number_of_columns == 0 && _number_of_rows == 0)
        return true;
  
    int *    Ti = _Ti.empty()? (int*)    malloc(0) : &_Ti[0];
    int *    Tj = _Tj.empty()? (int*)    malloc(0) : &_Tj[0];
    double * Tx = _Tx.empty()? (double*) malloc(0) : &_Tx[0];

    int*    Ap = (int*)    malloc ((_number_of_columns+1) * sizeof(int));
    int*    Ai = (int*)    malloc (_nz * sizeof(int));
    double* Ax = (double*) malloc (_nz * sizeof(double));
       
    double status = umfpack_di_triplet_to_col (_number_of_rows, _number_of_columns, _Ti.size(), Ti, Tj, Tx, Ap, Ai, Ax, (int *) NULL);

    if (status < 0)
    {
        cerr << "Failed triplet-to-column conversion. Error code = "<<status<<endl;
        exit(1);
    } 
    else 
    {
        _Ai.resize(_nz);
        _Ax.resize(_nz);
        _Ap.resize(_number_of_columns + 1);
        for (int i = 0; i < _nz; ++i) 
        {
            _Ai[i] = Ai[i];
            _Ax[i] = Ax[i];
        }
        for (int i = 0; i < _number_of_columns+1; i++)
            _Ap[i] = Ap[i];
    }

    free (Ax);
    free (Ai);
    free (Ap);

    return true;
}

//! This method checks if the sparse matrix is already lu-factorized or not 
//! \return bool
bool sparse_matrix::is_lu_factorized()
{
    return _lu_factorized;
}

//! This method returns a void pointer to the LU numeric factorization of the matrix
void * sparse_matrix::get_lu_numeric()
{
    return _lu_numeric;
}

//! This method returns a cholmod_factr pointer to the Cholesky factorization of the matrix
cholmod_factor * sparse_matrix::get_chol_L()
{
    return _chol_L;
}

//! This method returns the sum of the elements of a particular column
//! @param column [in]: index of the column to sum 
//! \return double
double sparse_matrix::get_column_sum (int column)
{
    double result = 0.0;

    int start_index = _Ap[column];
    int end_index   = _Ap[column+1];

    for ( int index = start_index; index < end_index; ++index)
        result += _Ax[index];

    return result;
}

//! This method returns a particular column in a sparse form 
//! @param column [in]: index of the column to retrieve 
//! \return pair<vector<int>, vector<double> >
pair<vector<int>, vector<double> > sparse_matrix::get_sparse_column (int column) 
{
    pair<vector<int>,vector<double> > col_index_val;

    int start_index = _Ap[column];
    int end_index   = _Ap[column+1];

    for (int index = start_index; index < end_index; ++index)
	{
	    col_index_val.first .push_back(_Ai[index]);
	    col_index_val.second.push_back(_Ax[index]);
    }
    return col_index_val;
}

//! This method destroy the void pointer to the LU Numeric factorization 
void sparse_matrix::destroy_lu_numeric()
{
    free(_lu_numeric);
}

//! This method destroy the pointer to the chol L factorization 
void sparse_matrix::destroy_chol_L()
{
    free(_chol_L);
}

//! This method returns the element in the matrix at the given index
//! @param row [in]: Row index of the element
//! @param col [in]: Column index of the element
//! \return double
double sparse_matrix::get_element(int row, int col)
{
    if (row < 0 || row >= _number_of_rows || col < 0 || col >= _number_of_columns)
    {
        cout<<"Location out of bounds. Exiting."<<endl;
        exit(1);
    }
    int start_index = _Ap[col];
    int end_index = _Ap[col+1];
  
    for ( int index = start_index; index < end_index && _Ai[index] <= row; ++index)
        if (_Ai[index] == row)
            return  _Ax[index];
    return 0;
}

// Setter Functions

//! This method sets the Ap vector of the compressed column format of the sprse matrix 
//! @param Ap [in]: Ap vector to set
void sparse_matrix::set_Ap(vector<int> & Ap)
{
    _Ap = Ap;
}

//! This method sets the Ai vector of the compressed column format of the sprse matrix 
//! @param Ai [in]: Ai vector to set
void sparse_matrix::set_Ai(vector<int> & Ai)
{
    _Ai = Ai;
}

//! This method sets the Ax vector of the compressed column format of the sprse matrix 
//! @param Ax [in]: Ax vector to set
void sparse_matrix::set_Ax(vector<double> & Ax)
{
    _Ax = Ax;
}

//! This method sets the Ti vector of the triplet format of the sprse matrix 
//! @param Ti [in]: Ti vector to set
void sparse_matrix::set_Ti(vector<int> & Ti)
{
    _Ti = Ti;
}

//! This method sets the Tj vector of the triplet format of the sprse matrix 
//! @param Tj [in]: Tj vector to set
void sparse_matrix::set_Tj(vector<int> & Tj)
{
    _Tj = Tj;
}

//! This method sets the Tx vector of the triplet format of the sprse matrix 
//! @param Tx [in]: Tx vector to set
void sparse_matrix::set_Tx(vector<double> & Tx)
{
    _Tx = Tx;
}

//! This method adds to the value of the element at the given index
//! @param row [in]: Row index of the element to insert
//! @param col [in]: Column index of the element to insert
//! @param val [in]: value of the element to insert
void sparse_matrix::insert_element(int row, int col, double val)
{
    //set_valid_collumns(false); (fix later)
    _Ti.push_back(row);
    _Tj.push_back(col);
    _Tx.push_back(val);
}

//! This method sets the element at the given index with a new value
//! @param row [in]: Row index of the element to replace
//! @param col [in]: Column index of the element to replace
//! @param val [in]: value of the element to replace
void sparse_matrix::replace_nonzero_entry(int row, int col, double val)
{
    if (row < 0 || row >= _number_of_rows || col < 0 || col >= _number_of_columns)
    {
        cout<<"Location out of bounds. Exiting."<<endl;
        exit(1);
    }
    int start_index = _Ap[col];
    int end_index   = _Ap[col+1];
  
    for ( int index = start_index; index < end_index && _Ai[index] <= row; ++index)
        if (_Ai[index] == row)
        {
            _Ax[index] = val;
            return;
        }

    // If entry is not zero, ignore function call
    cout<<"Matrix entry ("<<row<<","<<col<<") is zero. Use the function insert_element() to update it"<<endl;
}

//! This method sets the number of rows of the sparse matrix 
//! @param nrows [in]: number of rows of the matrix
void sparse_matrix::set_number_of_rows(int nrows)
{
    _number_of_rows = nrows;
}

//! This method sets the number of columns of the sparse matrix 
//! @param ncols [in]: number of columns of the matrix
void sparse_matrix::set_number_of_columns(int ncols)
{
    _number_of_columns = ncols;
}

//! This method sets the number of number of nonzeros of the sparse matrix 
//! @param nz [in]: number of nonzeros of the matrix
void sparse_matrix::set_nz(int nz)
{
    _nz = nz;   
}

//! This method sets the flag that indicates if the matrix is already lu-factorized or not
//! @param lu_factorize [in]: new value of the flag
void sparse_matrix::set_lu_factorized(bool lu_factorize)
{
    _lu_factorized = lu_factorize;
}

//! This method points the LU Numeric factorization void pointer to another pointer
//! @param lu_numeric [in]: new value of the void pointer
void sparse_matrix::set_lu_numeric(void * lu_numeric)
{
    _lu_numeric = lu_numeric;
}

//! This method points the Cholesky factorization pointer to another pointer
//! @param chol_L [in]: new value of the void pointer
void sparse_matrix::set_chol_L(cholmod_factor * chol_L)
{
    _chol_L = chol_L;
}

/* function to remove the fill_ins from a sparse matrix using cholmod_l_drop
 * based on seme tolerance */
//! This function removes the entries in the matrix whose absolute value is smaller than tol
//! @param tol [in]: a tolerance that determines which elements are to be dropped
void sparse_matrix::sparse_remove_fill_ins(double tol)
{
    /* Cholmod variables */
    cholmod_common Common, *cc;

    /* Matrix A */
    cholmod_sparse *chol_A;
    int *chol_Ai, *chol_Ap;
    double *chol_Ax;

    /* start cholmod */
    cc = &Common;
    cc->print = 5;
    cholmod_l_start(cc);

    /* create the cholmod matrix A */
    chol_A = cholmod_l_allocate_sparse(get_number_of_rows(),get_number_of_columns(), get_nz(), 0,1,0,CHOLMOD_REAL,cc);

    chol_Ai = (int*)chol_A->i;
    chol_Ap = (int*)chol_A->p;
    chol_Ax = (double*)chol_A->x;
    chol_A->nzmax = get_nz();

    for(int jj=0;jj<get_nz();jj++)
    {
        chol_Ai[jj] = get_Ai()[jj];
        chol_Ax[jj] = get_Ax()[jj];
    }

    for(int jj=0;jj<get_number_of_columns()+1;jj++)
        chol_Ap[jj] = get_Ap()[jj];

    /* drop the small entries */
    cholmod_l_drop(tol, chol_A, cc);

    /* Copy the result to our matrix class */
    chol_Ai = (int*)chol_A->i;
    chol_Ap = (int*)chol_A->p;
    chol_Ax = (double*)chol_A->x;

    /* clear the vectors before copying the new data to them */
    get_Ai().clear();
    get_Ap().clear();
    get_Ax().clear();

    /* this is when the result is 0 */
    if( (abs(chol_Ax[0]) <= tol ) && (chol_A->nzmax == 1))
    {
        set_nz(0);
        for(int jj=0;jj<get_number_of_columns()+1;jj++)
            get_Ap().push_back(chol_Ap[jj]);
    }
    else
    {
        set_nz(chol_A->nzmax);

        /* copy the compressed column format */
        for(unsigned int jj=0;jj<chol_A->nzmax;jj++)
        {
            get_Ai().push_back(chol_Ai[jj]);
            get_Ax().push_back(chol_Ax[jj]);
        }

        for(int jj=0;jj<get_number_of_columns()+1;jj++)
            get_Ap().push_back(chol_Ap[jj]);
    }

    /* free the cholmod structures */
    cholmod_l_free_sparse(&chol_A,cc);
    cholmod_l_finish(cc);
}

//! This method reads a metrix from a file
//! @param filename [in]: Name of the file to read the matrix from
void sparse_matrix::read_matrix_from_file(string filename)
{       
    // Open matrix file
    ifstream inv;
    inv.open(filename.c_str());

    // Helpers
    int    entry_int = 0;
    double entry_double = 0;
    _Ai.clear(); _Ap.clear(); _Ax.clear();
    
    // Read the first row 
    inv>>_number_of_rows>>_number_of_columns>>_nz;
    _lu_factorized = false;

    // Read Ai
    while(true)
    {
        inv>>entry_int;
        if (entry_int == -999)
            break;		
        _Ai.push_back(entry_int);
    }

    // Read Ap
    while(true)
    {
        inv>>entry_int;
        if (entry_int == -999)
            break;
        _Ap.push_back(entry_int);
    }

    // Read Ax
    while (true)
    {
        inv>>entry_double;
        if (entry_double == -999)
            break;
        _Ax.push_back(entry_double);
    }

    //Close the input file
    inv.close();
}

//! This methods removes the matrix from memory
void sparse_matrix::clear_matrix_memory()
{
    set_number_of_rows(0);
    set_number_of_columns(0);
    set_nz(0);
    set_lu_factorized(false);
    vector<int>()   .swap(_Ti);
    vector<int>()   .swap(_Tj);
    vector<double>().swap(_Tx);
    vector<int>()   .swap(_Ap);
    vector<int>()   .swap(_Ai);
    vector<double>().swap(_Ax);
}

// Overloading few operators
double sparse_matrix::operator()(int row, int col)
{
    return get_element(row, col);
}
