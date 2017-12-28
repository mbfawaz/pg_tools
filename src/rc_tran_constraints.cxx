#include <sstream>
#include <fstream>
#include "rc_tran_constraints.h"
#include "cholmod.h"
#include "string.h"
#include "omp.h"
#include "math.h"
#include "mosek_interface.h"
#include "rc_grid.h"
#include <cmath>

using namespace std;

vector<pair<double, string> > rc_tran_constraints_stamps;
void rc_tran_constraints_stamp(std::string op);

rc_tran_constraints::rc_tran_constraints()
{
    cholmod_start(&c);
}

rc_tran_constraints::rc_tran_constraints(rc_grid &g, string options_file)
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
        cout<<"ERROR: could not open rc_tran_constraints options file."<<endl;
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
            else if (strcmp(param_name, "CONTAINERS")      == 0) CONTAINERS      = param_value_double;
            else if (strcmp(param_name, "RANDOM_SEED")     == 0) RANDOM_SEED     = param_value_double;
            else if (strcmp(param_name, "MIN_T")           == 0) MIN_T           = param_value_double;
            else if (strcmp(param_name, "MAX_T")           == 0) MAX_T           = param_value_double;
            else if (strcmp(param_name, "ERROR_THRESHOLD") == 0) ERROR_THRESHOLD = param_value_double;
            else if (strcmp(param_name, "THREADS")         == 0) THREADS         = param_value_double;
        }
    }
    cout<<endl<<"Problem parameters as read from the problem options file"<<endl;
    cout      <<"--------------------------------------------------------"<<endl;
    cout<<"CONSTRAINTS                  = "<<CONSTRAINTS<<endl;
    cout<<"RANDOM_SEED                  = "<<RANDOM_SEED<<endl;
    cout<<"ERROR_THRESHOLD              = "<<ERROR_THRESHOLD<<endl;
    cout<<"MULTITHREADING               = "<<MULTITHREADING<<endl;
    cout<<"THREADS                      = "<<THREADS<<endl;
    cout<<endl;

    cout<<"Computing the transient bound..."<<endl;
    compute_transient_bound();

    // Free Memory
    cholmod_free_sparse(&G, &c);
    cholmod_free_sparse(&C, &c);
    cholmod_free_sparse(&B, &c);

    cholmod_finish(&c);
    cholmod_finish(&g.c);
}

void rc_tran_constraints::compute_transient_bound()
{
    rc_tran_constraints_stamp("Start");

    cout<<"Factorizing G..."<<endl;
    G_L = cholmod_analyze(G, &c);  // Symbolic factorization
    cholmod_factorize(G, G_L, &c); // Numeric  factorization
    rc_tran_constraints_stamp("Factorizing G");

    cout<<"Finding RC time step..."<<endl;
    find_rc_time_step();
    dtb = 2e-10;
    rc_tran_constraints_stamp("Finding RC time step");  
    cout<<"RC time step found = "<<dtb<<" s"<<endl;

    cout<<"Generating F0..."<<endl;
    F0.generate(SOURCES, CONSTRAINTS, 0.002, 5);
    rc_tran_constraints_stamp("Generating F0");

    cout<<"Generating other containers..."<<endl;
    generate_list_of_containers(Ft, F0, CONTAINERS);
    rc_tran_constraints_stamp("Generating other containers");

    time_steps.resize(CONTAINERS);
    for (int i = 0; i<CONTAINERS; i++)
    {
        time_steps[i] = floor((MIN_T + (MAX_T - MIN_T)*1.0*rand()/RAND_MAX/dtb)); 
    }
        
    cout<<"Computing A..."<<endl;
    compute_A();
    rc_tran_constraints_stamp("Computing A");

    cout<<"Computing the sparisification threshold xi..."<<endl;
    find_rc_sparsification_threshold(ERROR_THRESHOLD);
    rc_tran_constraints_stamp("Computing xi");
    cout<<"Sparsification threshold computed xi = "<<xi<<endl;

    cout<<"Factorizing A..."<<endl;
    A_L = cholmod_analyze(A, &c);  // Symbolic factorization
    cholmod_factorize(A, A_L, &c); // Numeric factorization
    rc_tran_constraints_stamp("Factorizing A");

    cout<<"Computing the emax vectors..."<<endl;
    vector<cholmod_dense* >  emax_vec(Ft.size());
    emax_vec = compute_Ainv_and_emax(Ft);
    c.print = 5;
    for (int i = 0; i<emax_vec.size(); i++)
    {
        cholmod_print_dense(emax_vec[i], "e", &c);
    }
    rc_tran_constraints_stamp("Computing the emax_vectors");

    cholmod_free_sparse(&F0.U, &F0.c);
    cholmod_finish(&F0.c);

    for (unsigned int i = 0; i<Ft.size(); i++)
           delete Ft[i];

    /* Compute v_bar_0 */
    cout<<"Computing the first upper bound vector..."<<endl;
    cholmod_dense * v_bar_0;
    cholmod_dense * Aemax_0 = smatrix_dmatrix_multiply(1, A, emax_vec[0]);

    v_bar_0 = cholmod_solve(CHOLMOD_A, G_L, Aemax_0, &c);
    cholmod_free_dense(&Aemax_0, &c);
    vector<cholmod_dense * > v_bar;
    v_bar.push_back(v_bar_0);
    rc_tran_constraints_stamp("Computing the first upper bound vector");

    cout<<"Computing the other upper bound vectors..."<<endl;
    compute_v_bar(v_bar, emax_vec);
    rc_tran_constraints_stamp("Computing the other upper bound vectors");
    
    for (unsigned int i = 0; i<emax_vec.size(); i++)
        cholmod_free_dense(&emax_vec[i], &c);
        
    FILE * fout;
    stringstream ss1;
    ss1<<"grid_"<<NODES<<".m";
    fout = fopen(ss1.str().c_str(), "w");

    fprintf(fout, "v_bar = [");
    for (int j = 0; j<NODES; j++)
    {
        for (unsigned int i = 0; i<v_bar.size(); i++)
        {
            double * x = (double*)v_bar[i]->x;
            fprintf(fout, "%.5f ", x[j]);
        }
        fprintf(fout, "\n");
    }
    fprintf(fout, "];\n\n\n");

    fprintf(fout, "t = [");
    for (unsigned int i = 0; i<time_steps.size(); i++)
    {
        fprintf(fout, "%i ", time_steps[i]);
    }
    fprintf(fout, "];\n\n\n");

    fclose(fout);

    printf("\nPrinting the profiler details [timer records]\
            \n---------------------------------------------\n");
    printf("%60s %11s\n", "Task Description", "Wall (s)");
    printf("%60s %11s\n", "----------------", "--------");
    for (unsigned int i = 1; i<rc_tran_constraints_stamps.size(); i++)
    {
        printf("%60s %11.2f\n", rc_tran_constraints_stamps[i].second.c_str(), rc_tran_constraints_stamps[i].first-rc_tran_constraints_stamps[i-1].first);
    }
    for (unsigned int i = 0; i<v_bar.size(); i++)
        cholmod_free_dense(&v_bar[i], &c);
    cholmod_free_sparse(&A, &c);
    cholmod_free_sparse(&G, &c);
    cholmod_free_sparse(&C, &c);
    cholmod_free_sparse(&B, &c);
    cholmod_free_factor(&A_L, &c);
    cholmod_free_factor(&G_L, &c);
    
}
void rc_tran_constraints::compute_v_bar(vector<cholmod_dense*> & v_bar, vector<cholmod_dense*> & emax_vec)
{
    cholmod_dense * delta_p;
    cholmod_dense * rhs;
    double * delta_px;
    double * emax_vecx;

    for (unsigned int i = 0; i<time_steps.size(); i++)
    {
        for (int j = 0; j< time_steps[i]; j++)
        {
            rhs = smatrix_dmatrix_multiply(1, B, v_bar[v_bar.size()-1]);

            delta_p = cholmod_solve(CHOLMOD_A, A_L, rhs, &c);

            delta_px  = (double*)delta_p->x;
            emax_vecx = (double*)emax_vec[i+1]->x;
            
            cholmod_dense * v_bar_p = cholmod_zeros(NODES, 1, CHOLMOD_REAL, &c);
            double * v_bar_px;
            v_bar_px  = (double*)v_bar_p->x;
            for (int k = 0; k<NODES; k++)
                v_bar_px[k] = emax_vecx[k] + delta_px[k];

            v_bar.push_back(v_bar_p);
            cholmod_free_dense(&rhs, &c);
            cholmod_free_dense(&delta_p, &c);
        }
    }
}


vector<cholmod_dense*> rc_tran_constraints::compute_Ainv_and_emax(vector<container *> &F) //Always parallel
{
    vector<cholmod_dense *> emax_vec(Ft.size());
    vector<double *> emax_vecx(Ft.size());
    for (unsigned int i = 0;i<Ft.size(); i++)
    {
        emax_vec[i]  = cholmod_zeros(NODES, 1, CHOLMOD_REAL, &c);
        emax_vecx[i] = (double*) emax_vec[i]->x;
    }

    int MAX_TASK_LOAD = ceil(1.0*NODES/THREADS);
    long sum = 0;

    #pragma omp parallel num_threads(THREADS)
    {
        int tid = omp_get_thread_num();
        int LOAD = MAX_TASK_LOAD;
        if (tid == THREADS-1)
            LOAD = NODES - (THREADS-1)*MAX_TASK_LOAD;
       
        // Setup RHS
        cholmod_dense * RHS = cholmod_zeros(NODES, 1, CHOLMOD_REAL, &c); 
        cholmod_dense * col_ainv;
        double * RHSx = (double*) RHS->x;

        // Setup LP dimensions
        vector<lp_opt_problem> lp(Ft.size());

        for (unsigned int j = 0; j<Ft.size(); j++)
        {
            lp[j].set_numvar(Ft[j]->DIMENSION);
            lp[j].set_numcon(Ft[j]->CONSTRAINTS);

            // Create environment and task
            lp[j].create_env();
            lp[j].create_opt_task();

            // Choose the optimization method
            lp[j].set_opt_method(SIMPLEX);

            // Set up constraints
            int * Ap = (int*) Ft[j]->U->p;
            int * Ai = (int*) Ft[j]->U->i;
            double * Ax = (double*) Ft[j]->U->x;
            lp[j].set_A_2(Ap, Ai, Ax);

            vector<double> lower_local(Ft[j]->DIMENSION, 0);
            vector<double> lower_global(Ft[j]->CONSTRAINTS, 0);
            lp[j].set_bounds_on_constraints(lower_global, Ft[j]->gb);
            lp[j].set_bounds_on_variables  (lower_local , Ft[j]->lb );

            // Load data into task
            lp[j].load_constraints_into_task();

            // Set optimization sense (max or min)
            lp[j].set_opt_sense(MAXIMIZE);
        }
        
        vector<double> obj;

        for (int i = tid*LOAD; i<(tid+1)*LOAD; i++)
        {
            RHSx[i] = 1;
            col_ainv = cholmod_solve(CHOLMOD_A, A_L, RHS, &c);
            double * col_ainvx = (double*)col_ainv->x;
            RHSx[i] = 0;

            for (unsigned int j = 0; j<Ft.size(); j++)
            {
                obj.resize(Ft[j]->DIMENSION, 0);
                for (int k = 0; k<SOURCES; k++)
                {
//                     if (abs(col_ainvx[csi[k]]) > xi)
                        obj[k] = col_ainvx[csi[k]];
                        cout<<obj[k]<<endl;
                }
                lp[j].set_obj(obj);

                lp[j].load_objective_into_task();

                // Optimize
                lp[j].optimize();
                
                // Collect the results
                emax_vecx[j][i] = lp[j].get_optimized_objective();

//                 #pragma omp critical(PRINT)
//                 {
//                     sum++;
//                     if (sum % (int)(NODES*Ft.size()/100) == 0)
//                     {
//                         printf("\r%2d%% complete...", (int)(100.0*sum/(NODES*Ft.size())));
//                         fflush(stdout);
//                         printf("\r");
//                     }
//                 }
            }
            cholmod_free_dense(&col_ainv, &c);
        }
        
        for (unsigned int j = 0; j<lp.size(); j++)
            lp[j].delete_task_and_env();
            
        cholmod_free_dense(&col_ainv, &c);
        cholmod_free_dense(&RHS, &c);
    }

    return emax_vec;
}

void rc_tran_constraints::generate_list_of_containers(vector<container*> &Ft, container &F0, int N)
{
    //push back F0
    container * Fnew = new container(F0);
    Ft.push_back(Fnew);

    //create the other spaces
    for (int k = 0; k<N; k++)
    {
        // Constructor
        container * F = new container(F0);
        
        // Copy local bounds 
        for (int i = 0; i<F->DIMENSION; i++)
            F->lb[i] = (0.75 + 0.5*rand()/RAND_MAX) * F->lb[i];

        // Copy global bounds 
        for (int i = 0; i<F->CONSTRAINTS; i++)
            F->gb[i] = (0.75 + 0.5*rand()/RAND_MAX) * F->gb[i];

        Ft.push_back(F);
    }
}


bool rc_tran_constraints::check_validity(container &F)
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

void rc_tran_constraints::compute_A()
{
    // Find B
    double * Bx = (double*)B->x; 
    for (int i = 0; i<NODES; i++)
        Bx[i]/=dtb;
        
    // Find A
    A = smatrix_smatrix_add(1, G, 1, B);
    A->stype = 1;
}

void rc_tran_constraints::find_rc_time_step()
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

void rc_tran_constraints::find_rc_sparsification_threshold(double zeta) //zeta = voltage error threshold (in Volts)
{
    // Prep
    int nv = NODES;

    // Find delta
    cholmod_dense* ones = cholmod_ones(nv, 1, CHOLMOD_REAL, &c);
    cholmod_dense* A_ones = smatrix_dmatrix_multiply(1, A, ones);
    cholmod_dense* y = cholmod_solve(CHOLMOD_A, G_L, A_ones , &c); cholmod_free_dense(&A_ones, &c);
    double inf_norm = 0;
    double * y_x = (double*)y->x;
    for (int i = 0; i<nv; i++)
        if (inf_norm < abs(y_x[i]))
            inf_norm = abs(y_x[i]);
    cholmod_free_dense(&y, &c); //we're done with y

    double delta = zeta/inf_norm;

    // Find xi, which is the sparsification threshold essentially
    double sum_il = 0;
    for (unsigned int i = 0; i<F0.lb.size(); i++)
        sum_il += F0.lb[i];

    xi = delta/sum_il;
    
    cholmod_free_dense(&ones, &c);
    cholmod_free_dense(&A_ones, &c);
    cholmod_free_dense(&y, &c);
}


cholmod_sparse * rc_tran_constraints::smatrix_smatrix_add(double a1, cholmod_sparse *A1,
                                                    double a2, cholmod_sparse *A2)
{
    double alpha[1] = {a1}; 
    double beta [1] = {a2}; 
    cholmod_sparse * A3 = cholmod_add(A1, A2, alpha, beta, true, false, &c);
    return A3;
}

cholmod_dense  * rc_tran_constraints::smatrix_dmatrix_multiply(double a1, cholmod_sparse *A1,
                                                         cholmod_dense  *A2)       
{
    cholmod_dense *A3 = cholmod_zeros(A1->nrow, A2->ncol, CHOLMOD_REAL, &c);
    double alpha[1] = {a1};
    double beta [1] = {0};
    cholmod_sdmult(A1, 0, alpha, beta, A2, A3, &c);
    return A3;
}

void rc_tran_constraints_stamp(string op)
{
    pair<double, string> st;
    st.first  = omp_get_wtime();
    st.second = op;
    rc_tran_constraints_stamps.push_back(st);
}




/*    // Now, find the exact upper bound without sparsification
    xi = 0;
    cout<<"Computing the (exact) emax vectors..."<<endl;
    cholmod_dense*  emax_exact;
    emax_exact = compute_Ainv_and_emax(&F);
    rc_tran_constraints_stamp("Computing the (exact) emax_vectors");

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


