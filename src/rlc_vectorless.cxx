#include <sstream>
#include <fstream>
#include "rlc_vectorless.h"
#include "cholmod.h"
#include "string.h"
#include "omp.h"
#include "math.h"
#include "mosek_interface.h"
#include "rlc_grid.h"
#include <cmath>

using namespace std;

vector<pair<double, string> > rlc_vectorless_stamps;
void rlc_vectorless_stamp(std::string op);

rlc_vectorless::rlc_vectorless()
{
    cholmod_start(&c);
}

rlc_vectorless::rlc_vectorless(rlc_grid &g, string options_file)
{
    cholmod_start(&c);

    // Parameters of an RLC grid
    // Read straight from the variable g
    G                   = g.G;
    C                   = g.C;
    L                   = g.L;
    M                   = g.M;
    csi                 = g.csi;
    NODES               = g.number_of_nodes;
    SOURCES             = g.number_of_sources;
    INDUCTORS           = g.number_of_inductors;

    // Loading options from the rlc vectorless options file
    ifstream infile;
    infile.open(options_file.c_str());
    if (!infile)
    {
        cout<<"ERROR: could not open rlc_vectorless options file."<<endl;
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

    eye_nv = cholmod_speye(NODES,     NODES,     CHOLMOD_REAL, &c);
    eye_nl = cholmod_speye(INDUCTORS, INDUCTORS, CHOLMOD_REAL, &c);

    cout<<"Computing the bound..."<<endl;
    compute_upper_and_lower_bounds();

    // Free Memory
    cholmod_free_sparse(&eye_nv, &c);
    cholmod_free_sparse(&eye_nl, &c);
    cholmod_free_sparse(&G, &c);
    cholmod_free_sparse(&C, &c);
    cholmod_free_sparse(&L, &c);
    cholmod_free_sparse(&M, &c);
    cholmod_free_sparse(&Mt, &c);
    cholmod_free_sparse(&B, &c);
    cholmod_free_sparse(&Einv, &c);

    printf("\nPrinting the profiler details [timer records]\
            \n---------------------------------------------\n");
    printf("%60s %11s\n", "Task Description", "Wall (s)");
    printf("%60s %11s\n", "----------------", "--------");
    for (unsigned int i = 1; i<rlc_vectorless_stamps.size(); i++)
    {
        printf("%60s %11.2f\n", rlc_vectorless_stamps[i].second.c_str(), rlc_vectorless_stamps[i].first-rlc_vectorless_stamps[i-1].first);
    }

    cholmod_finish(&c);
    cholmod_finish(&g.c);
}

void rlc_vectorless::compute_upper_and_lower_bounds()
{
    rlc_vectorless_stamp("Start");

    cout<<"Finding RLC time step..."<<endl;
    find_rlc_time_step(-10, -8);
    rlc_vectorless_stamp("Finding RLC time step");  
    cout<<"RLC time step found = "<<dtb<<" s"<<endl;

    cout<<"Generating the feasible space F..."<<endl;
    srand(RANDOM_SEED);
    F.generate(SOURCES, CONSTRAINTS, 0.002, 5);
    rlc_vectorless_stamp("Generating the feasible space F");

//  if (!check_validity(F))
//      cout<<"The container F is not valid"<<endl;
//  else 
//      cout<<"The container F are valid"<<endl;
    
    // Find F22p and F22n
    cout<<"Finding F22p and F22n..."<<endl;
    F22p = smatrix_smatrix_add(0.5, F22,  0.5, abs_F22); //cholmod_drop(1e-7, F22p, &c);
    F22m = smatrix_smatrix_add(0.5, F22, -0.5, abs_F22); //cholmod_drop(1e-7, F22m, &c);
    rlc_vectorless_stamp("Finding F22p and F22n");

//  cout<<"Computing the sparisification threshold xi..."<<endl;
//  find_rlc_sparsification_threshold(ERROR_THRESHOLD);
//  rlc_vectorless_stamp("Computing xi");
//  cout<<"Sparsification threshold computed xi = "<<xi<<endl;

    cout<<"Computing the eopt vector..."<<endl;
    cholmod_dense*  eopt;
    eopt = compute_Dinv_and_eopt(&F);
    rlc_vectorless_stamp("Computing the eopt vector");

    double * eopt_x = (double*) eopt->x;
    int nv = NODES;
    int nl = INDUCTORS;
    int n = nv + nl;
       
    cholmod_free_sparse(&F22, &c);
    cholmod_free_sparse(&abs_F22, &c);

    cholmod_free_sparse(&F22, &c);
    cholmod_free_sparse(&abs_F22, &c);

    // Iterative Solver
    cout<<"Computing the final bound..."<<endl;
    cholmod_dense * et1 = cholmod_zeros(nv, 1, CHOLMOD_REAL, &c);
    cholmod_dense * et2 = cholmod_zeros(nl, 1, CHOLMOD_REAL, &c);
    cholmod_dense * eb1 = cholmod_zeros(nv, 1, CHOLMOD_REAL, &c);
    cholmod_dense * eb2 = cholmod_zeros(nl, 1, CHOLMOD_REAL, &c);
    
    double * et1x = (double*) et1->x;
    double * et2x = (double*) et2->x;
    double * eb1x = (double*) eb1->x;
    double * eb2x = (double*) eb2->x;

    for (int i = 0; i<nv; i++)
    {       
        et1x[i] = eopt_x[i];
        eb1x[i] = eopt_x[n+i];
    }   
    for (int i = 0; i<nl; i++)
    {       
        et2x[i] = eopt_x[nv+i];
        eb2x[i] = eopt_x[n+nv+i];
        cout<<eb2x[i]<<" ";
    }   

    cholmod_dense * v_bar_u, * i_bar_u, * v_bar_l, * i_bar_l;
    solve_eye_minus_F_tilde_eq_f( et1,  et2,  eb1,  eb2, 
                                  &v_bar_u,  &i_bar_u,  &v_bar_l,  &i_bar_l,
                                  et1,  et2,  eb1,  eb2);

    rlc_vectorless_stamp("Computing the upper bound vector");

//  double * v_bar_ux = (double*)v_bar_u->x;
//  double * v_bar_lx = (double*)v_bar_l->x;

    cholmod_free_dense(&eopt, &c);
    cholmod_free_sparse(&F.U, &F.c);
    cholmod_free_sparse(&D, &c);
    cholmod_free_factor(&D_L, &c);
    cholmod_free_sparse(&A, &c);
    cholmod_free_sparse(&Y, &c);
    cholmod_free_sparse(&F22, &c);
    cholmod_free_sparse(&abs_F22, &c);
    cholmod_free_sparse(&F22m, &c);
    cholmod_free_sparse(&F22p, &c);
    cholmod_finish(&F.c);

    cholmod_free_dense(&v_bar_u, &c);
    cholmod_free_dense(&i_bar_u, &c);
    cholmod_free_dense(&v_bar_l, &c);
    cholmod_free_dense(&i_bar_l, &c);
    cholmod_free_dense(&et1, &c);
    cholmod_free_dense(&et2, &c);
    cholmod_free_dense(&eb1, &c);
    cholmod_free_dense(&eb2, &c);
    

}

cholmod_dense* rlc_vectorless::compute_Dinv_and_eopt(container * F)
{
    cholmod_dense * eopt;
    double * eoptx;
    eopt  = cholmod_zeros(2*NODES+2*INDUCTORS, 1, CHOLMOD_REAL, &c);
    eoptx = (double*) eopt->x;

    if (MULTITHREADING)
    { 
//      vector<long> nz(THREADS, 0);
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
            cholmod_dense * col_dinv;
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
                col_dinv = cholmod_solve(CHOLMOD_A, D_L, RHS, &c);
                double * col_dinvx = (double*)col_dinv->x;

                RHSx[i] = 0;

                obj.clear();
                obj.resize(F->DIMENSION, 0);
                for (int k = 0; k<SOURCES; k++)
//                 if (col_dinvx[csi[k]] - 1e-5 >= 0)
//                 {
                       obj[k] = col_dinvx[csi[k]];
//                     nz[tid]++;
//                 }
                lp.set_obj(obj);

                lp.load_objective_into_task();

                // Optimize
                lp.optimize();
                
                // Collect the results
                eoptx[i] = lp.get_optimized_objective();
                eoptx[NODES+INDUCTORS+i] = 0;

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
                cholmod_free_dense(&col_dinv, &c);
            }
            
            lp.delete_task_and_env();
                
            cholmod_free_dense(&col_dinv, &c);
            cholmod_free_dense(&RHS, &c);
        }
        
        sum = 0;
        task_loads.clear(); start_points.clear();
        task_loads.resize(THREADS, INDUCTORS/THREADS);
        for (int i = 0; i<INDUCTORS%THREADS; i++)
            task_loads[i]++;

        start_points.resize(THREADS); start_points[0] = 0;
        for (unsigned int i = 1; i<start_points.size(); i++)
        {
            start_points[i] = start_points[i-1] + task_loads[i-1];
        }
        
        THREADS = min(THREADS, INDUCTORS);
        #pragma omp parallel num_threads(THREADS)
        {
            int tid = omp_get_thread_num();

            //Setup RH2
            cholmod_dense * RHS2 = cholmod_zeros(INDUCTORS, 1, CHOLMOD_REAL, &c);
            cholmod_dense * col_y;
            double * RHS2x = (double*) RHS2->x;
            double * Einvx = (double*) Einv->x; // this is the same as the diagonal
 
            // Second part of eopt (opt of -E^{-1} M^T D^{-1} H i)
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

            lp.set_opt_sense(MINIMIZE);  //will only do minimizations here

            vector<double> obj;
           
            for (int i = start_points[tid]; i<start_points[tid] + task_loads[tid]; i++)
            {
                eoptx[NODES + i] = 0; // we know that the maximization will yield a zero
                RHS2x[i] = 1;
                col_y = smatrix_dmatrix_multiply(1, Y, RHS2);

                double * col_yx = (double*)col_y->x;

                RHS2x[i] = 0;
                
                obj.clear();
                obj.resize(F->DIMENSION, 0);
                for (int k = 0; k<SOURCES; k++)
                    obj[k] = - Einvx[i] * col_yx[csi[k]];
                
                lp.set_obj(obj);

                lp.load_objective_into_task();

                // Optimize
                lp.optimize();
                
                // Collect the results
                eoptx[2*NODES + INDUCTORS + i] = lp.get_optimized_objective();

                #pragma omp critical(PRINT)
                {
                    sum++;
                    printf("\r%2d%% complete...", (int)(100.0*sum/INDUCTORS));
                    fflush(stdout);
                    printf("\r");
                }

                cholmod_free_dense(&col_y, &c);
            }
            lp.delete_task_and_env();
            cholmod_free_dense(&col_y, &c);
            cholmod_free_dense(&RHS2, &c);
        }

        // Nonzeros:
//      long total_nz = 0;
//      for (unsigned int i = 0; i<nz.size(); i++)
//          total_nz += nz[i];
//      cout<<"Percentage of nonzeros left = "<<100.0*total_nz/(NODES*SOURCES)<<"%"<<endl;
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

        // First part of each half of eopt (opt of D^{-1} H i)
        // Setup RHS1
        cholmod_dense * RHS1 = cholmod_zeros(NODES, 1, CHOLMOD_REAL, &c); 
        cholmod_dense * col_dinv;
        double * RHS1x = (double*) RHS1->x;
        for (int i = 0; i<NODES; i++)
        {
            RHS1x[i] = 1;
            col_dinv = cholmod_solve(CHOLMOD_A, D_L, RHS1, &c);
            double * col_dinvx = (double*)col_dinv->x;
            
            RHS1x[i] = 0;

            obj.clear();
            obj.resize(F->DIMENSION, 0);
            for (int k = 0; k<SOURCES; k++)
//              if (col_dinvx[csi[k]] - xi >= 0)
                    obj[k] = col_dinvx[csi[k]];
            
            lp.set_obj(obj);

            lp.load_objective_into_task();

            // Optimize
            lp.optimize();
            
            // Collect the results
            eoptx[i] = lp.get_optimized_objective();
            eoptx[NODES + INDUCTORS + i] = 0; // we know that the minimization will yield a zero

            printf("\r%2d%% complete...", (int)(100.0*i/NODES));
            fflush(stdout);
            printf("\r");

            cholmod_free_dense(&col_dinv, &c);
        }
        cholmod_free_dense(&RHS1, &c);
        cholmod_free_dense(&col_dinv, &c);

        // Second part of eopt (opt of -E^{-1} M^T D^{-1} H i)
        lp.set_opt_sense(MINIMIZE);  //will only do minimizations here
        cholmod_dense * RHS2 = cholmod_zeros(INDUCTORS, 1, CHOLMOD_REAL, &c);
        cholmod_dense * col_y;
        double * RHS2x = (double*) RHS2->x;
        double * Einvx = (double*) Einv->x; // this is the same as the diagonal
        for (int i = 0; i<INDUCTORS; i++)
        {
            eoptx[NODES + i] = 0; // we know that the maximization will yield a zero
            RHS2x[i] = 1;
            col_y = smatrix_dmatrix_multiply(1, Y, RHS2);
            double * col_yx = (double*)col_y->x;

            RHS2x[i] = 0;
            
            obj.clear();
            obj.resize(F->DIMENSION, 0);
            for (int k = 0; k<SOURCES; k++)
                obj[k] = - Einvx[i] * col_yx[csi[k]];
            
            lp.set_obj(obj);

            lp.load_objective_into_task();

            // Optimize
            lp.optimize();
            
            // Collect the results
            eoptx[2*NODES + INDUCTORS + i] = lp.get_optimized_objective();

            printf("\r%2d%% complete...", (int)(100.0*i/INDUCTORS));
            fflush(stdout);
            printf("\r");

            cholmod_free_dense(&col_y, &c);
        }

        lp.delete_task_and_env();
            
        cholmod_free_dense(&col_y, &c);
        cholmod_free_dense(&RHS2, &c);
    }
    return eopt;
}


bool rlc_vectorless::check_validity(container &F)
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

void rlc_vectorless::compute_D()
{
    // Find B 
    B = cholmod_copy_sparse(C, &c);
    double * Bx = (double*)B->x; 
    for (int i = 0; i<NODES; i++)
        Bx[i]/=dtb;

    // Find E
    Einv = cholmod_copy_sparse(L, &c);
    double * Einvx = (double*)Einv->x; 
    for (int i = 0; i<INDUCTORS; i++)
        Einvx[i] = dtb/Einvx[i];

    // Find Mt
    Mt = transpose(M);

    // Find A
    A = smatrix_smatrix_add(1, G, 1, B);

    // Find M_E_inv
    cholmod_sparse * M_Einv    = smatrix_smatrix_multiply(M, Einv);
    cholmod_sparse * M_Einv_Mt = smatrix_smatrix_multiply(M_Einv, Mt);

    // Find D
    D = smatrix_smatrix_add(1, A, 1, M_Einv_Mt);
    D->stype = 1;

    cholmod_free_sparse(&M_Einv, &c);
    cholmod_free_sparse(&M_Einv_Mt, &c);
}

/* UTILITIES */
void rlc_vectorless::solve_eye_minus_F_tilde_eq_f
            (cholmod_dense *  ft1, cholmod_dense  * ft2, cholmod_dense *  fb1, cholmod_dense *  fb2,
             cholmod_dense ** xt1, cholmod_dense ** xt2, cholmod_dense ** xb1, cholmod_dense ** xb2,
             cholmod_dense * x0t1, cholmod_dense * x0t2, cholmod_dense * x0b1, cholmod_dense * x0b2)
{
    
   
    double one[1] = {1}; double minus_one[1] = {-1};

    *xt1 = cholmod_copy_dense(x0t1, &c); // cholmod_ones(nv, 1, CHOLMOD_REAL, &c);
    *xt2 = cholmod_copy_dense(x0t2, &c); // cholmod_ones(l , 1, CHOLMOD_REAL, &c);
    *xb1 = cholmod_copy_dense(x0b1, &c); // cholmod_ones(nv, 1, CHOLMOD_REAL, &c);
    *xb2 = cholmod_copy_dense(x0b2, &c); // cholmod_ones(l , 1, CHOLMOD_REAL, &c);

    cholmod_dense * xt1_prev = cholmod_copy_dense(x0t1, &c);
    cholmod_dense * xt2_prev = cholmod_copy_dense(x0t2, &c);
    cholmod_dense * xb1_prev = cholmod_copy_dense(x0b1, &c);
    cholmod_dense * xb2_prev = cholmod_copy_dense(x0b2, &c);

    cholmod_dense * zt1, *zt2, *zb1, *zb2;

    double sum_norm_of_f =  cholmod_norm_dense(ft1, 1, &c)
                          + cholmod_norm_dense(ft2, 1, &c)
                          + cholmod_norm_dense(fb1, 1, &c)
                          + cholmod_norm_dense(fb2, 1, &c);


    for (int k = 0; k<20; k++)
    {
        // x = F_tilde * x
        F_tilde_times_x(*xt1,*xt2, *xb1, *xb2, &zt1, &zt2, &zb1, &zb2); 

        cholmod_free_dense(xt1, &c);
        cholmod_free_dense(xt2, &c);
        cholmod_free_dense(xb1, &c);
        cholmod_free_dense(xb2, &c);

        *xt1 = cholmod_copy_dense(zt1, &c);
        *xt2 = cholmod_copy_dense(zt2, &c);
        *xb1 = cholmod_copy_dense(zb1, &c);
        *xb2 = cholmod_copy_dense(zb2, &c);

        // x = f + x
        cholmod_sdmult(eye_nv, 0, one, one, ft1, *xt1, &c);
        cholmod_sdmult(eye_nl , 0, one, one, ft2, *xt2, &c);
        cholmod_sdmult(eye_nv, 0, one, one, fb1, *xb1, &c);
        cholmod_sdmult(eye_nl , 0, one, one, fb2, *xb2, &c);
        
        // Find residue
//      double *xt1x = (double*)((*xt1)->x);
//      double *xt1x_prev = (double*)(xt1_prev->x);
        cholmod_sdmult(eye_nv, 0, one, minus_one, *xt1, xt1_prev, &c);
        cholmod_sdmult(eye_nl , 0, one, minus_one, *xt2, xt2_prev, &c);
        cholmod_sdmult(eye_nv, 0, one, minus_one, *xb1, xb1_prev, &c);
        cholmod_sdmult(eye_nl , 0, one, minus_one, *xb2, xb2_prev, &c);



        double residue =  cholmod_norm_dense(xt1_prev, 1, &c)
                        + cholmod_norm_dense(xt2_prev, 1, &c)
                        + cholmod_norm_dense(xb1_prev, 1, &c)
                        + cholmod_norm_dense(xb2_prev, 1, &c);
        

//      cout<<1.0*residue/sum_norm_of_f<<endl;
    
        if (1.0*residue/sum_norm_of_f < 1e-3)
        {
            cholmod_free_dense(&xt1_prev, &c);       
            cholmod_free_dense(&xt2_prev, &c);       
            cholmod_free_dense(&xb1_prev, &c);       
            cholmod_free_dense(&xb2_prev, &c);       
            
            cholmod_free_dense(&zt1, &c);       
            cholmod_free_dense(&zt2, &c);       
            cholmod_free_dense(&zb1, &c);       
            cholmod_free_dense(&zb2, &c);       

            break;
        }
        cholmod_free_dense(&xt1_prev, &c);       
        cholmod_free_dense(&xt2_prev, &c);       
        cholmod_free_dense(&xb1_prev, &c);       
        cholmod_free_dense(&xb2_prev, &c);       

        cholmod_free_dense(&zt1, &c);       
        cholmod_free_dense(&zt2, &c);       
        cholmod_free_dense(&zb1, &c);       
        cholmod_free_dense(&zb2, &c);       

        // x_prev = x 
        xt1_prev = cholmod_copy_dense(*xt1, &c);
        xt2_prev = cholmod_copy_dense(*xt2, &c);
        xb1_prev = cholmod_copy_dense(*xb1, &c);
        xb2_prev = cholmod_copy_dense(*xb2, &c);
    }
//     cout<<endl;
}


void  rlc_vectorless::F_tilde_times_x(cholmod_dense *  xt1, cholmod_dense *  xt2, cholmod_dense *  xb1, cholmod_dense *  xb2,
                                      cholmod_dense ** yt1, cholmod_dense ** yt2, cholmod_dense ** yb1, cholmod_dense ** yb2)
{
    double one[1] = {1}; double zero[1] = {0};
    int nl = INDUCTORS;

    // Finding yt1 //
    cholmod_dense * D_yt1 = smatrix_dmatrix_multiply(1, B, xt1); // B*xt1
    cholmod_sdmult(M, 0, one, one, xt2, D_yt1, &c); // D_yt1 = B*xt1 + M*xt2
    *yt1 = cholmod_solve(CHOLMOD_A, D_L, D_yt1, &c); // yt1 = inv(D)*D_yt1
    cholmod_free_dense(&D_yt1, &c);

    // Finding yb1 //
    cholmod_dense * D_yb1 = smatrix_dmatrix_multiply(1, B, xb1); // B*xb1
    cholmod_sdmult(M, 0, one, one, xb2, D_yb1, &c); // D_yb1 = B*xb1 + M*xb2
    *yb1 = cholmod_solve(CHOLMOD_A, D_L, D_yb1, &c); // yb1 = inv(D)*D_yb1 
    cholmod_free_dense(&D_yb1, &c);

    cholmod_dense * B_xt1 = smatrix_dmatrix_multiply(1, B, xt1);
    cholmod_dense * B_xb1 = smatrix_dmatrix_multiply(1, B, xb1);

    // Finding yt2 //
    // Set yt2 = -Ei*Y'*B*xb1
    cholmod_dense * Y_B_xb1 = cholmod_zeros(nl, xb1->ncol, CHOLMOD_REAL, &c);
    cholmod_sdmult(Y, 1, one, zero, B_xb1, Y_B_xb1, &c);
    *yt2 = smatrix_dmatrix_multiply(-1, Einv, Y_B_xb1); 
    cholmod_free_dense(&Y_B_xb1, &c);
    // Set yt2 = F22m*xb2 + yt2
    cholmod_sdmult(F22m, 0, one, one, xb2, *yt2, &c); // yt2 = F22m*xb2 + yt2
    // Set yt2 = F22p*xt2 + yt2
    cholmod_sdmult(F22p, 0, one, one, xt2, *yt2, &c); // yt2 = F22p*xt2 + yt2

    // Finding yb2 //
    // Set yb2 = -Einv*Y'*B*xt1
    cholmod_dense * Y_B_xt1 = cholmod_zeros(nl, xt1->ncol, CHOLMOD_REAL, &c);
    cholmod_sdmult(Y, 1, one, zero, B_xt1, Y_B_xt1, &c);
    *yb2 = smatrix_dmatrix_multiply(-1, Einv, Y_B_xt1); 
    cholmod_free_dense(&Y_B_xt1, &c);
    // Set yb1 = F22m*xt2 + yb1
    cholmod_sdmult(F22m, 0, one, one, xt2, *yb2, &c); // yb2 = F22m*xt2 + yb2
    // Set yb1 = F22p*xt2 + yb1
    cholmod_sdmult(F22p, 0, one, one, xb2, *yb2, &c); // yb2 = F22p*xb2 + yb2

    // Free //
    cholmod_free_dense(&B_xb1, &c);
    cholmod_free_dense(&B_xt1, &c);
}


void rlc_vectorless::find_rlc_time_step(double log_dtb_min, double log_dtb_max)
{
    // bisection method
    srand(NODES);
    double a1 = log_dtb_min;
    double b1 = log_dtb_max;
    double p1 = a1 + (b1 - a1)*rand()/RAND_MAX;
    
    dtb = pow(10, 1);
    while (true)
    {
        if (is_spec_radius_less_than(p1) == 0)
        {   
            a1 = p1;
        }
        else
            return;
        p1 = (a1 + b1)/2;
        dtb = pow(10, p1);
    }
}

int rlc_vectorless::is_spec_radius_less_than(double target)
{
    int nl = INDUCTORS;
    int nv = NODES;
    compute_D();

    D_L = cholmod_analyze(D, &c);  // Symbolic factorization
    cholmod_factorize(D, D_L, &c); // Numeric factorization

    cholmod_dense * M_dense = cholmod_sparse_to_dense(M, &c);
    cholmod_dense * Y_dense = cholmod_solve(CHOLMOD_A, D_L, M_dense, &c);
    cholmod_free_dense(&M_dense, &c);
    Y  = cholmod_dense_to_sparse(Y_dense, 1, &c);
    cholmod_free_dense(&Y_dense, &c);

    // Find F22
    cholmod_sparse *Einv_Mt  = smatrix_smatrix_multiply(Einv, Mt);
    cholmod_sparse * Mt_Y    = smatrix_smatrix_multiply(Mt, Y);
    cholmod_sparse * Einv_Mt_Y = smatrix_smatrix_multiply(Einv, Mt_Y); cholmod_free_sparse(&Mt_Y, &c);
    F22     = smatrix_smatrix_add     (1, eye_nl, -1, Einv_Mt_Y); 
    cholmod_free_sparse(&Einv_Mt_Y, &c);
    cholmod_free_sparse(&Einv_Mt, &c);

    // Find abs_F22
    abs_F22 = cholmod_copy_sparse(F22, &c);
    double * abs_F22x = (double*)abs_F22->x;
    for (unsigned int i = 0; i<abs_F22->nzmax; i++)
        abs_F22x[i] = abs(abs_F22x[i]);
   
    cholmod_dense * q1 = cholmod_ones(nv, 1, CHOLMOD_REAL, &c);
    cholmod_dense * q2 = cholmod_ones(nl, 1, CHOLMOD_REAL, &c);
    cholmod_dense * q1_next, * q2_next;
    double lambda_min, lambda_max, ratio;    
    while (true)
    {
        /*First step - Multiply absolue_F by [q1; q2] */
        abs_F_times_x(q1, q2, &q1_next, &q2_next);

        /* Now, we find min (q_next/q) and max(q_next,q) */
        double * q1x = (double*) q1->x;
        double * q2x = (double*) q2->x;
        double * q1_nextx = (double*) q1_next->x;
        double * q2_nextx = (double*) q2_next->x;

        lambda_min = q1_nextx[0]/q1x[0];
        lambda_max = q1_nextx[0]/q1x[0];

        for (int i = 1; i<nv; i++)
        {
            ratio = q1_nextx[i]/q1x[i];
            if (lambda_min >= ratio)
                lambda_min = ratio;
            if (lambda_max <= ratio)
                lambda_max = ratio;
        }
        for (int i = 0; i<nl; i++)
        {
            ratio = q2_nextx[i]/q2x[i];
            if (lambda_min >= ratio)
                lambda_min = ratio;
            if (lambda_max <= ratio)
                lambda_max = ratio;
        }
//      cout<<lambda_min<<" to "<<lambda_max<<endl;
 
        if (lambda_min >= target)
        {
            cholmod_free_dense(&q1, &c);
            cholmod_free_dense(&q2, &c);
            cholmod_free_dense(&q1_next, &c);
            cholmod_free_dense(&q2_next, &c);
            return 0;
        }
        if (target >= lambda_max)
        {
            cholmod_free_dense(&q1, &c);
            cholmod_free_dense(&q2, &c);
            cholmod_free_dense(&q1_next, &c);
            cholmod_free_dense(&q2_next, &c);
            return 1;
        }

        q1 = cholmod_copy_dense(q1_next, &c); // cholmod_ones(nv, 1, CHOLMOD_REAL, &c);
        q2 = cholmod_copy_dense(q2_next, &c); // cholmod_ones(nv, 1, CHOLMOD_REAL, &c);
    }
}

void  rlc_vectorless::abs_F_times_x(cholmod_dense *  x1, cholmod_dense *  x2,
                                    cholmod_dense ** y1, cholmod_dense ** y2)
{
    double one[1] = {1}; double zero[1] = {0};
    int nl = INDUCTORS;

    // Finding y1 //
    cholmod_dense * D_y1 = smatrix_dmatrix_multiply(1, B, x1); // B*x1
    cholmod_sdmult(M, 0, one, one, x2, D_y1, &c); // D_y1 = B*x1 + M*x2
    *y1 = cholmod_solve(CHOLMOD_A, D_L, D_y1, &c); // y1 = inv(D)*D_y1

    // Finding y2 //
    // Set y2 = -Ei*Y'*B*x1
    cholmod_dense * B_x1 = smatrix_dmatrix_multiply(1, B, x1);
    cholmod_dense * Y_B_x1 = cholmod_zeros(nl, x1->ncol, CHOLMOD_REAL, &c);
    cholmod_sdmult(Y, 1, one, zero, B_x1, Y_B_x1, &c);
    *y2 = smatrix_dmatrix_multiply(1, Einv, Y_B_x1); 

    // Set y2 = abs_F22*x2 + y2
    cholmod_sdmult(abs_F22, 0, one, one, x2, *y2, &c); // y2 = F22m*x2 + y2

    // Free //
    cholmod_free_dense(&D_y1, &c);
    cholmod_free_dense(&B_x1, &c);
    cholmod_free_dense(&Y_B_x1, &c);
}


// void rlc_vectorless::find_rlc_sparsification_threshold(double zeta) //zeta = voltage error threshold (in Volts)
// {
//     // Prep
//     int nv = NODES;
//     int nl = INDUCTORS;
// 
//     // Find delta
// 
//     cholmod_dense * ut1 = cholmod_ones(nv, 1, CHOLMOD_REAL, &c);
//     cholmod_dense * ut2 = cholmod_ones(nl, 1, CHOLMOD_REAL, &c);
//     cholmod_dense * ub1 = cholmod_ones(nv, 1, CHOLMOD_REAL, &c);
//     cholmod_dense * ub2 = cholmod_ones(nl, 1, CHOLMOD_REAL, &c);
//     double * ub1x = (double*) ub1->x;
//     double * ub2x = (double*) ub2->x;
//     for (int i = 0; i<nv; i++)
//         ub1x[i] *= -1;
//     for (int i = 0; i<nl; i++)
//         ub2x[i] *= -1;
// 
//     cholmod_dense * zt1, * zt2, * zb1, * zb2;
//     F_tilde_times_x(ut1, ut2, ub1, ub2, &zt1, &zt2, &zb1, &zb2); 
// 
//     double inf_norm = 0;
//     double * zt1x = (double*) zt1->x;
//     double * zt2x = (double*) zt2->x;
//     double * zb1x = (double*) zb1->x;
//     double * zb2x = (double*) zb2->x;
// 
//     for (int i = 0; i<nv; i++)
//         if (inf_norm < abs(zt1x[i]))
//             inf_norm = abs(zt1x[i]);
// 
//     for (int i = 0; i<nl; i++)
//         if (inf_norm < abs(zt2x[i]))
//             inf_norm = abs(zt2x[i]);
// 
//     for (int i = 0; i<nv; i++)
//         if (inf_norm < abs(zb1x[i]))
//             inf_norm = abs(zb1x[i]);
// 
//     for (int i = 0; i<nl; i++)
//         if (inf_norm < abs(zb2x[i]))
//             inf_norm = abs(zb2x[i]);
// 
//     cholmod_free_dense(&zt1, &c); 
//     cholmod_free_dense(&zt2, &c); 
//     cholmod_free_dense(&zb1, &c); 
//     cholmod_free_dense(&zb2, &c); 
//     cholmod_free_dense(&ut1, &c); 
//     cholmod_free_dense(&ut2, &c); 
//     cholmod_free_dense(&ub1, &c); 
//     cholmod_free_dense(&ub2, &c); 
// 
//     double delta = zeta*(1-inf_norm);
// 
//     // Find xi, which is the sparsification threshold essentially
//     double sum_il = 0;
//     for (unsigned int i = 0; i<F.lb.size(); i++)
//         sum_il += F.lb[i];
// 
//     xi = delta/sum_il;
// }


cholmod_sparse * rlc_vectorless::smatrix_smatrix_add(double a1, cholmod_sparse *A1,
                                                     double a2, cholmod_sparse *A2)
{
    double alpha[1] = {a1}; 
    double beta [1] = {a2}; 
    cholmod_sparse * A3 = cholmod_add(A1, A2, alpha, beta, true, false, &c);
    return A3;
}

cholmod_dense  * rlc_vectorless::smatrix_dmatrix_multiply(double a1, cholmod_sparse *A1,
                                                          cholmod_dense  *A2)       
{
    cholmod_dense *A3 = cholmod_zeros(A1->nrow, A2->ncol, CHOLMOD_REAL, &c);
    double alpha[1] = {a1};
    double beta [1] = {0};
    cholmod_sdmult(A1, 0, alpha, beta, A2, A3, &c);
    return A3;
}

cholmod_sparse * rlc_vectorless::smatrix_smatrix_multiply(cholmod_sparse *A1, 
                                                          cholmod_sparse *A2)
{
    cholmod_sparse * A3 = cholmod_ssmult(A1, A2, 0, 1, 1, &c);
    return A3;
}

cholmod_sparse * rlc_vectorless::transpose(cholmod_sparse * A1)
{
    cholmod_sparse * A1t = cholmod_allocate_sparse(A1->ncol, A1->nrow, A1->nzmax, 0, 1, 0, CHOLMOD_REAL, &c);
    cholmod_transpose_unsym(A1, 1, NULL, NULL, 0, A1t, &c);
    return A1t;
}

void rlc_vectorless_stamp(string op)
{
    pair<double, string> st;
    st.first  = omp_get_wtime();
    st.second = op;
    rlc_vectorless_stamps.push_back(st);
}



/*    cout<<"Computing D..."<<endl;
    compute_D();
    rlc_vectorless_stamp("Computing D");
        
    cout<<"Factorizing D..."<<endl;
    D_L = cholmod_analyze(D, &c);  // Symbolic factorization
    cholmod_factorize(D, D_L, &c); // Numeric factorization
    rlc_vectorless_stamp("Factorizing D");

    cout<<"Computing Y = D^{-1}M..."<<endl; // might parallelize later
    cholmod_dense * M_dense = cholmod_sparse_to_dense(M, &c);
    cholmod_dense * Y_dense = cholmod_solve(CHOLMOD_A, D_L, M_dense, &c);
    cholmod_free_dense(&M_dense, &c);
    Y  = cholmod_dense_to_sparse(Y_dense, 1, &c);
    cholmod_free_dense(&Y_dense, &c);
    rlc_vectorless_stamp("Computing Y = D^{-1}M");

    // Find F22
    cholmod_sparse *Einv_Mt  = smatrix_smatrix_multiply(Einv, Mt);
    cout<<"Finding F22..."<<endl;
    cholmod_sparse * Mt_Y    = smatrix_smatrix_multiply(Mt, Y);
    cholmod_sparse * Einv_Mt_Y = smatrix_smatrix_multiply(Einv, Mt_Y); cholmod_free_sparse(&Mt_Y, &c);
    cholmod_sparse * F22     = smatrix_smatrix_add     (1, eye_nl, -1, Einv_Mt_Y); 
    cholmod_free_sparse(&Einv_Mt_Y, &c);
    cholmod_free_sparse(&Einv_Mt, &c);
    rlc_vectorless_stamp("Finding F22");

    // Find abs_F22
    cout<<"Finding abs_F22..."<<endl;
    abs_F22 = cholmod_copy_sparse(F22, &c);
    double * abs_F22x = (double*)abs_F22->x;
    for (unsigned int i = 0; i<abs_F22->nzmax; i++)
        abs_F22x[i] = abs(abs_F22x[i]);
    rlc_vectorless_stamp("Finding abs_F22");
*/
