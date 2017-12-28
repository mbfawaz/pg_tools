#ifndef _SPARSE_MATRIX_
#define _SPARSE_MATRIX_

#include <vector>
#include <string>
#include "cholmod.h"

using namespace std;

namespace matrix
{
    // Classes
    class sparse_matrix;

    /* Linear system solvers */
    bool            lu_factorize                    (sparse_matrix &A);
    void            chol_factorize                  (sparse_matrix &A);
    vector<double>  backward_forward_substitution    (sparse_matrix &A, vector<double> &b);
    vector<double>  chol_solve                      (sparse_matrix &A, vector<double> &b);
    
    /* Find C = alpha*A + beta*B */
    void            sparse_matrix_add               (double alpha, sparse_matrix &A,
                                                     double beta,  sparse_matrix &B,
                                                                   sparse_matrix &C);

    /* Find C = A*B */
    void            sparse_matrix_mult              (sparse_matrix &A,
                                                    sparse_matrix &B,
                                                    sparse_matrix &C);

    /* post-multiply a sparse matrix by vector */
    vector<double>  sparse_matrix_vector_mult       (sparse_matrix &A , vector<double> &b);

    /* pre- multiply a sparse matrix by vector */
    vector<double>  sparse_vector_matrix_mult       (vector<double> &b, sparse_matrix &A );

    /* drop either the positive or the negative elements */
    /* 1 -> drops the positive entries, -1 -> drops the negative entries */
    void            sparse_drop_entries             (sparse_matrix &A, int sign, sparse_matrix &C);

    /* transpose matrix A and put the result in C */
    void            transpose                       (sparse_matrix &A , sparse_matrix &C);

    void            sparse_matrix_horizontal_augment(sparse_matrix &A, sparse_matrix &B, sparse_matrix &C);
    void            sparse_matrix_vertical_augment  (sparse_matrix &A, sparse_matrix &B, sparse_matrix &C);
};


class matrix::sparse_matrix
{
public:
    // Default Constructor
    sparse_matrix();

    // Regular Constructors
    sparse_matrix(int number_of_rows, int number_of_columns);
    sparse_matrix(int number_of_rows, int number_of_columns, int nz);

    // Copy Constructor
    sparse_matrix(sparse_matrix & sparse_matrix_to_copy);
    
    // Overload Operators
    virtual sparse_matrix & operator= (sparse_matrix & sparse_matrix_to_copy);

    // Destroyer
    ~sparse_matrix();

    // Getter Functions
    virtual vector<int>     &       get_Ap                            ();
    virtual vector<int>     &       get_Ai                            ();
    virtual vector<double>  &       get_Ax                            ();
    virtual vector<int>     &       get_Ti                            ();
    virtual vector<int>     &       get_Tj                            ();
    virtual vector<double>  &       get_Tx                            ();
    virtual int                     get_number_of_rows                ();
    virtual int                     get_number_of_columns             ();
    virtual int                     get_nz                            ();
    virtual double                  get_element                       (int row, int column);
    static  int                     get_global_sparse_matrix_iterator ();
    virtual int                     get_local_sparse_matrix_iterator  ();
    virtual bool                    get_triplet_format                ();
    virtual bool                    get_column_format                 ();
    virtual bool                    is_lu_factorized                  ();
    virtual bool                    is_chol_factorized                ();
    virtual void *                  get_lu_numeric                    ();
    virtual cholmod_factor *        get_chol_L                        ();
    virtual double                  get_column_sum                    (int column);
    virtual pair<vector<int>, 
                 vector<double> >   get_sparse_column                 (int column);
   
    // Setter Functions
    virtual void                    set_Ap                (vector<int> & Ap);
    virtual void                    set_Ai                (vector<int> & Ai);
    virtual void                    set_Ax                (vector<double> & Ax);
    virtual void                    set_Ti                (vector<int> & Ti);
    virtual void                    set_Tj                (vector<int> & Tj);
    virtual void                    set_Tx                (vector<double> & Tx);
    virtual void                    insert_element        (int row, int col, double val);
    virtual void                    replace_nonzero_entry (int row, int col, double val);
    virtual void                    set_number_of_rows    (int nrows);
    virtual void                    set_number_of_columns (int ncols);
    virtual void                    set_nz                (int nz);
    virtual void                    set_lu_factorized     (bool lu_factorized);
    virtual void                    set_chol_factorized   (bool chol_factorized);
    virtual void                    set_lu_numeric        (void * lu_numeric);
    virtual void                    set_chol_L            (cholmod_factor * chol_L);
    virtual void                    sparse_remove_fill_ins(double tol);
    virtual void                    destroy_lu_numeric    ();
    virtual void                    destroy_chol_L        ();
    virtual void                    read_matrix_from_file (string filename);
    virtual void                    clear_matrix_memory   ();

protected:
    // Dimensions
    int _number_of_rows;
    int _number_of_columns;
    int _nz;

    // Matrix iterator
    static int _global_sparse_matrix_iterator;
           int _local_sparse_matrix_iterator;

    // LU factorization numeric pointer
    void           * _lu_numeric;
    cholmod_factor * _chol_L;

    // Flags
    bool   _lu_factorized;
    bool   _chol_factorized;

    // CCS format
    vector<int>    _Ap;
    vector<int>    _Ai;
    vector<double> _Ax;

    // Triplet format
    vector<int>    _Ti;
    vector<int>    _Tj;
    vector<double> _Tx;
};

#endif
