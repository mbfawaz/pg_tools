#include <sstream>
#include <fstream>
#include "rc_vectorless.h"
#include "cholmod.h"
#include "string.h"
#include "omp.h"
#include "math.h"
#include "mosek_interface.h"
#include "rc_grid.h"
#include <cmath>

using namespace std;

vector<pair<double, string> > rc_vectorless_stamps;
void rc_vectorless_stamp(std::string op);

rc_vectorless::rc_vectorless()
{
    cholmod_start(&c);
}

rc_vectorless::rc_vectorless(rc_grid &g, string options_file)
{
    cholmod_start(&c);

    // Parameters of an RC grid
    // Read straight from the variable g
    G                   = g.G;
    C                   = g.C;
    B                   = cholmod_copy_sparse(g.C, &c);
    csi                 = g.csi;
    NODES               = g.number_of_nodes;
    SOURCES             = g.number_of_sources;

    // Loading options from the rc vectorless options file
    ifstream infile;
    infile.open(options_file.c_str());
    if (!infile)
    {
        cout<<"ERROR: could not open rc_vectorless options file."<<endl;
        exit(1);
    }
    while (!infile.eof())
    {
        char str[255];
        infile.getline(str, 255);
        if (strcmp(str, ".end") == 0)
            break;

        char *  param_name         = strtok(str,  " =");
        char *  param_value        = strtok(NULL, " =");

        if (param_value != NULL)
        {   
            double  param_value_double = atof(param_value);

                 if (strcmp(param_name, "CONSTRAINTS")     == 0) CONSTRAINTS     = param_value_double;
            else if (strcmp(param_name, "RANDOM_SEED")     == 0) RANDOM_SEED     = param_value_double;
//          else if (strcmp(param_name, "ERROR_THRESHOLD") == 0) ERROR_THRESHOLD = param_value_double;
            else if (strcmp(param_name, "MULTITHREADING")  == 0) MULTITHREADING  = param_value_double;
            else if (strcmp(param_name, "THREADS")         == 0) THREADS         = param_value_double;
        }
    }
    cout<<endl<<"Problem parameters as read from the problem options file"<<endl;
    cout      <<"--------------------------------------------------------"<<endl;
    cout<<"CONSTRAINTS                  = "<<CONSTRAINTS<<endl;
    cout<<"RANDOM_SEED                  = "<<RANDOM_SEED<<endl;
//  cout<<"ERROR_THRESHOLD              = "<<ERROR_THRESHOLD<<endl;
    cout<<"MULTITHREADING               = "<<MULTITHREADING<<endl;
    cout<<"THREADS                      = "<<THREADS<<endl;
    cout<<endl;

    cout<<"Computing the bound..."<<endl;
    compute_upper_bound();

    // Free Memory
    cholmod_free_sparse(&G, &c);
    cholmod_free_sparse(&C, &c);
    cholmod_free_sparse(&B, &c);

    cholmod_finish(&c);
    cholmod_finish(&g.c);
}

void rc_vectorless::compute_upper_bound()
{
    rc_vectorless_stamp("Start");

    cout<<"Factorizing G..."<<endl;
    G_L = cholmod_analyze(G, &c);  // Symbolic factorization
    cholmod_factorize(G, G_L, &c); // Numeric  factorization
    rc_vectorless_stamp("Factorizing G");

    cout<<"Finding RC time step..."<<endl;
    find_rc_time_step();
    rc_vectorless_stamp("Finding RC time step");  
    cout<<"RC time step found = "<<dtb<<" s"<<endl;

    cout<<"Generating F..."<<endl;
    srand(RANDOM_SEED);
    F.generate(SOURCES, CONSTRAINTS, 0.002/(3*35), 5);
    rc_vectorless_stamp("Generating F");

    if (!check_validity(F))
        cout<<"The container F is not valid"<<endl;
    else 
        cout<<"The container F are valid"<<endl;

    cout<<"Computing A..."<<endl;
    compute_A();
    rc_vectorless_stamp("Computing A");

//  cout<<"Computing the sparisification threshold xi..."<<endl;
//  find_rc_sparsification_threshold(ERROR_THRESHOLD);
//  rc_vectorless_stamp("Computing xi");
//  cout<<"Sparsification threshold computed xi = "<<xi<<endl;

    cout<<"Factorizing A..."<<endl;
    A_L = cholmod_analyze(A, &c);  // Symbolic factorization
    cholmod_factorize(A, A_L, &c); // Numeric factorization
    rc_vectorless_stamp("Factorizing A");

    cout<<"Computing the emax vector..."<<endl;
    cholmod_dense*  emax;
    emax = compute_Ainv_and_emax(&F);
    rc_vectorless_stamp("Computing the emax vector");

    cout<<"Computing the upper bound vector..."<<endl;

    cholmod_dense * Aemax = smatrix_dmatrix_multiply(1, A, emax);

    cholmod_dense * v_bar = cholmod_solve(CHOLMOD_A, G_L, Aemax, &c);
    double * v_bar_x = (double*) v_bar->x;

    // Find max upper bound
    double max_ub = 0;
    for (int i = 0; i<NODES; i++)
    {   
        if (max_ub < v_bar_x[i])
            max_ub = v_bar_x[i];
    }
    cout<<"Max upper bound = "<<max_ub<<endl;

    cholmod_free_dense(&emax, &c);
    cholmod_free_sparse(&F.U, &F.c);
    cholmod_free_sparse(&A, &c);
    cholmod_finish(&F.c);
    cholmod_free_dense(&Aemax, &c);
    cholmod_free_factor(&A_L, &c);
    cholmod_free_factor(&G_L, &c);
    rc_vectorless_stamp("Computing the upper bound vector");

    cholmod_free_dense(&v_bar, &c);

    printf("\nPrinting the profiler details [timer records]\
            \n---------------------------------------------\n");
    printf("%60s %11s\n", "Task Description", "Wall (s)");
    printf("%60s %11s\n", "----------------", "--------");
    for (unsigned int i = 1; i<rc_vectorless_stamps.size(); i++)
    {
        printf("%60s %11.2f\n", rc_vectorless_stamps[i].second.c_str(), rc_vectorless_stamps[i].first-rc_vectorless_stamps[i-1].first);
    }
}

cholmod_dense* rc_vectorless::compute_Ainv_and_emax(container * F)
{
    cholmod_dense * emax;
    double * emaxx;
    emax  = cholmod_zeros(NODES, 1, CHOLMOD_REAL, &c);
    emaxx = (double*) emax->x;

    if (MULTITHREADING)
    { 
        // We need to find an optimal distribution:
        // Assign floor(NODES/THREADS) + 1 to the first NODES%THREADS threads
        // Assign floor(NODES/THREADS)     to the rest
        vector<int> task_loads(THREADS, NODES/THREADS);
        for (int i = 0; i<NODES%THREADS; i++)
            task_loads[i]++;

        vector<int> start_points(THREADS); start_points[0] = 0;
        for (unsigned int i = 1; i<start_points.size(); i++)
        {
            start_points[i] = start_points[i-1] + task_loads[i-1];
        }

        int sum = 0;
        #pragma omp parallel num_threads(THREADS)
        {
            int tid = omp_get_thread_num();

            // Setup RHS
            cholmod_dense * RHS = cholmod_zeros(NODES, 1, CHOLMOD_REAL, &c); 
            cholmod_dense * col_ainv;
            double * RHSx = (double*) RHS->x;

            // Setup LP dimensions
            lp_opt_problem lp;

            lp.set_numvar(F->DIMENSION);
            lp.set_numcon(F->CONSTRAINTS);

            // Create environment and task
            lp.create_env();
            lp.create_opt_task();

            // Choose the optimization method
            lp.set_opt_method(SIMPLEX);

            // Set up constraints
            int * Ap = (int*) F->U->p;
            int * Ai = (int*) F->U->i;
            double * Ax = (double*) F->U->x;
            lp.set_A_2(Ap, Ai, Ax);

            vector<double> lower_local(F->DIMENSION, 0);
            vector<double> lower_global(F->CONSTRAINTS, 0);
            lp.set_bounds_on_constraints(lower_global, F->gb);
            lp.set_bounds_on_variables  (lower_local , F->lb );

            // Load data into task
            lp.load_constraints_into_task();

            // Set optimization sense (max or min)
            lp.set_opt_sense(MAXIMIZE);
        
            vector<double> obj;

            for (int i = start_points[tid]; i<start_points[tid] + task_loads[tid]; i++)
            {
                RHSx[i] = 1;
                col_ainv = cholmod_solve(CHOLMOD_A, A_L, RHS, &c);
                double * col_ainvx = (double*)col_ainv->x;

                RHSx[i] = 0;

                obj.clear();
                obj.resize(F->DIMENSION, 0);
                for (int k = 0; k<SOURCES; k++)
//                     if (col_ainvx[csi[k]] - xi >= 0)
//                     {
                        obj[k] = col_ainvx[csi[k]];
//                         nz[tid]++;
//                     }
                lp.set_obj(obj);

                lp.load_objective_into_task();

                // Optimize
                lp.optimize();
                
                // Collect the results
                emaxx[i] = lp.get_optimized_objective();

                #pragma omp critical(PRINT)
                {
                    sum++;
                    if (sum % (int)(NODES/100) == 0)
                    {
                        printf("\r%2d%% complete...", (int)(100.0*sum/NODES));
                        fflush(stdout);
                        printf("\r");
                    } 
                }
                cholmod_free_dense(&col_ainv, &c);
            }
            
            lp.delete_task_and_env();
                
            cholmod_free_dense(&col_ainv, &c);
            cholmod_free_dense(&RHS, &c);
        }
        // Nonzeros:
//         long total_nz = 0;
//         for (unsigned int i = 0; i<nz.size(); i++)
//             total_nz += nz[i];
//         cout<<"Percentage of nonzeros left = "<<100.0*total_nz/(NODES*SOURCES)<<"%"<<endl;
    }
    else
    {
        // Setup LP dimensions
        lp_opt_problem lp;

        lp.set_numvar(F->DIMENSION);
        lp.set_numcon(F->CONSTRAINTS);

        // Create environment and task
        lp.create_env();
        lp.create_opt_task();

        // Choose the optimization method
        lp.set_opt_method(SIMPLEX);

        // Set up constraints
        int * Ap = (int*) F->U->p;
        int * Ai = (int*) F->U->i;
        double * Ax = (double*) F->U->x;
        lp.set_A_2(Ap, Ai, Ax);

        vector<double> lower_local(F->DIMENSION, 0);
        vector<double> lower_global(F->CONSTRAINTS, 0);
        lp.set_bounds_on_constraints(lower_global, F->gb);
        lp.set_bounds_on_variables  (lower_local , F->lb );

        // Load data into task
        lp.load_constraints_into_task();

        // Set optimization sense (max or min)
        lp.set_opt_sense(MAXIMIZE);
    
        vector<double> obj;

        // Setup RHS
        cholmod_dense * RHS = cholmod_zeros(NODES, 1, CHOLMOD_REAL, &c); 
        cholmod_dense * col_ainv;
        double * RHSx = (double*) RHS->x;
        for (int i = 0; i<NODES; i++)
        {
            RHSx[i] = 1;
            col_ainv = cholmod_solve(CHOLMOD_A, A_L, RHS, &c);
            double * col_ainvx = (double*)col_ainv->x;
            
            RHSx[i] = 0;

            obj.clear();
            obj.resize(F->DIMENSION, 0);
            for (int k = 0; k<SOURCES; k++)
//              if (col_ainvx[csi[k]] - xi >= 0)
                    obj[k] = col_ainvx[csi[k]];
            
            lp.set_obj(obj);

            lp.load_objective_into_task();

            // Optimize
            lp.optimize();
            
            // Collect the results
            emaxx[i] = lp.get_optimized_objective();

            printf("\r%2d%% complete...", (int)(100.0*i/NODES));
            fflush(stdout);
            printf("\r");

            cholmod_free_dense(&col_ainv, &c);
        }
        
        lp.delete_task_and_env();
            
        cholmod_free_dense(&col_ainv, &c);
        cholmod_free_dense(&RHS, &c);
    }
    return emax;
}

bool rc_vectorless::check_validity(container &F)
{
    cholmod_dense * corner = cholmod_zeros(F.DIMENSION, 1, CHOLMOD_REAL, &c);
    for (int i = 0; i<F.DIMENSION; i++)
    {
        double * cx = (double*) corner->x;
        cx[i] = F.lb[i];
        cholmod_dense* temp = smatrix_dmatrix_multiply(1, F.U, corner);  
        double * tx = (double*) temp->x;
        for (int j = 0; j<F.CONSTRAINTS; j++)
            if (tx[j] > F.gb[j])
                return false;
        cx[i] = 0;
        cholmod_free_dense(&temp, &c);
    }   
    cholmod_free_dense(&corner, &c);
    return true;
}

void rc_vectorless::compute_A()
{
    // Find B
    double * Bx = (double*)B->x; 
    for (int i = 0; i<NODES; i++)
        Bx[i]/=dtb;

    // Find A
    A = smatrix_smatrix_add(1, G, 1, B);
}

void rc_vectorless::find_rc_time_step()
{
    int nv = NODES;
    double l = 0, l_prev = 0;
    double * xk, * C_xk, * xkp1;
    cholmod_dense *x;
    cholmod_dense *b = cholmod_ones(nv, 1, CHOLMOD_REAL, &c);
    double * bx = (double* ) b->x; 
    
    // Allocation
    xk   = (double*) calloc    (nv, sizeof(double));
    C_xk = (double*) calloc    (nv, sizeof(double));
    xkp1 = (double*) malloc    (nv* sizeof(double));
 
    for (int i = 0; i<nv; i++)   xk[i] = 1.0*rand()/RAND_MAX;
    
    int    * Cp = (int*)   C->p;
    int    * Ci = (int*)   C->i;
    double * Cx = (double*)C->x;

    // Main loop
    while (true)
    {
        double x_norm_2 = 0; 
        for (int i = 0; i<nv; i++)
            x_norm_2 += xk[i]*xk[i];

        for (int i = 0; i<nv; i++)
            xk[i] *= 1.0/x_norm_2;            // xk = xk / |xk|_2
        
        for(int j = 0; j < nv; j++)
            for(int q = Cp[j]; q < Cp[j+1]; q++)
                C_xk[Ci[q]] += Cx[q]*xk[j];

        for (int i = 0; i<nv; i++)
            bx[i] = C_xk[i];
        x = cholmod_solve(CHOLMOD_A, G_L, b, &c); 
        xkp1 = (double*) x->x;

        double num = 0, den = 0;
        for (int i =0; i<nv; i++)
        {
            num += xkp1[i]*xk[i];
            den += xk  [i]*xk[i];
        }
        l = num/den;                          // l = (xkp1*xk)/(xk*xk);
        if (abs(l - l_prev)/l < 1e-3)          // Check convergence - free and return
        {   
            cholmod_free_dense(&b, &c);
            cholmod_free_dense(&x, &c);
            free(xk); 
            free(C_xk); 
            dtb = l;
            break;
        }
        l_prev = l;

        for (int i = 0; i<nv; i++)             
        {
            xk[i]   = xkp1[i];                // New guess for the evector
            C_xk[i] = 0;                      // Reinitialize C_xk
        }
        cholmod_free_dense(&x, &c);
    }
}

// void rc_vectorless::find_rc_sparsification_threshold(double zeta) //zeta = voltage error threshold (in Volts)
// {
//     // Prep
//     int nv = NODES;
// 
//     // Find delta
//     cholmod_dense* ones = cholmod_ones(nv, 1, CHOLMOD_REAL, &c);
//     cholmod_dense* A_ones = smatrix_dmatrix_multiply(1, A, ones);
//     cholmod_dense* y = cholmod_solve(CHOLMOD_A, G_L, A_ones , &c); cholmod_free_dense(&A_ones, &c);
//     double inf_norm = 0;
//     double * y_x = (double*)y->x;
//     for (int i = 0; i<nv; i++)
//         if (inf_norm < abs(y_x[i]))
//             inf_norm = abs(y_x[i]);
//     cholmod_free_dense(&y, &c); //we're done with y
// 
//     double delta = zeta/inf_norm;
// 
//     // Find xi, which is the sparsification threshold essentially
//     double sum_il = 0;
//     for (unsigned int i = 0; i<F.lb.size(); i++)
//         sum_il += F.lb[i];
// 
//     xi = delta/sum_il;
//     
//     cholmod_free_dense(&ones, &c);
//     cholmod_free_dense(&A_ones, &c);
//     cholmod_free_dense(&y, &c);
// }


cholmod_sparse * rc_vectorless::smatrix_smatrix_add(double a1, cholmod_sparse *A1,
                                                    double a2, cholmod_sparse *A2)
{
    double alpha[1] = {a1}; 
    double beta [1] = {a2}; 
    cholmod_sparse * A3 = cholmod_add(A1, A2, alpha, beta, true, false, &c);
    return A3;
}

cholmod_dense  * rc_vectorless::smatrix_dmatrix_multiply(double a1, cholmod_sparse *A1,
                                                         cholmod_dense  *A2)       
{
    cholmod_dense *A3 = cholmod_zeros(A1->nrow, A2->ncol, CHOLMOD_REAL, &c);
    double alpha[1] = {a1};
    double beta [1] = {0};
    cholmod_sdmult(A1, 0, alpha, beta, A2, A3, &c);
    return A3;
}

void rc_vectorless_stamp(string op)
{
    pair<double, string> st;
    st.first  = omp_get_wtime();
    st.second = op;
    rc_vectorless_stamps.push_back(st);
}




/*    // Now, find the exact upper bound without sparsification
    xi = 0;
    cout<<"Computing the (exact) emax vectors..."<<endl;
    cholmod_dense*  emax_exact;
    emax_exact = compute_Ainv_and_emax(&F);
    rc_vectorless_stamp("Computing the (exact) emax_vectors");

    cout<<"Computing the (exact) upper bound vector..."<<endl;
    cholmod_dense * v_bar_exact;

    cholmod_dense *Aemax_exact = cholmod_zeros(NODES, 1, CHOLMOD_REAL, &c);
    cholmod_sdmult(A, 0, alpha, beta, emax_exact, Aemax_exact, &c);

    v_bar_exact = cholmod_solve(CHOLMOD_A, G_L, Aemax_exact, &c);
    double * v_bar_exact_x = (double*) v_bar_exact->x;

    // Finding the error vector
    vector<double> error(NODES);
    for (int i = 0; i<NODES; i++)
        error[i] = abs(v_bar_exact_x[i] - v_bar_x[i]);

    // Find the max of the error
    double max_error = 0 ;
    for (int i = 0; i<NODES; i++)
        if (max_error < error[i])
            max_error = error[i];

    cout<<"Max absolute error = "<<max_error<<endl;*/


