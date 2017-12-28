#include <sstream>
#include <fstream>
#include "rc_simulation.h"
#include "cholmod.h"
#include "string.h"
#include "omp.h"
#include "math.h"
#include "mosek_interface.h"
#include "rc_grid.h"
#include "waveforms.h"
#include "rc_tr.h"
#include <cmath>

using namespace std;

vector<pair<double, string> > rc_simulation_stamps;
void rc_simulation_stamp(std::string op);

rc_simulation::rc_simulation()
{
    cholmod_start(&c);
}

rc_simulation::rc_simulation(rc_grid &g, string options_file)
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

    // Loading options from the rc simulation options file
    ifstream infile;
    infile.open(options_file.c_str());
    if (!infile)
    {
        cout<<"ERROR: could not open rc_simulation options file."<<endl;
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

                 if (strcmp(param_name, "RANDOM_SEED")     == 0) RANDOM_SEED     = param_value_double;
            else if (strcmp(param_name, "NUM_TIME_POINTS") == 0) NUM_TIME_POINTS = param_value_double;
            else if (strcmp(param_name, "MULTITHREADING")  == 0) MULTITHREADING  = param_value_double;
            else if (strcmp(param_name, "THREADS")         == 0) THREADS         = param_value_double;
            else if (strcmp(param_name, "ENV_TYPE")        == 0) ENV_TYPE        = param_value_double;
            else if (strcmp(param_name, "TR")              == 0) TR              = param_value_double;
        }
    }
    cout<<endl<<"Problem parameters as read from the problem options file"<<endl;
    cout      <<"--------------------------------------------------------"<<endl;
    cout<<"RANDOM_SEED                  = "<<RANDOM_SEED<<endl;
    cout<<"NUM_TIME_POINTS              = "<<NUM_TIME_POINTS<<endl;
    cout<<"MULTITHREADING               = "<<MULTITHREADING<<endl;
    cout<<"THREADS                      = "<<THREADS<<endl;
    cout<<"TR                           = "<<TR<<endl;

    // Waveforms
    cout<<endl<<"Generating waveforms"<<endl;
    cout<<      "--------------------"<<endl;
    W.generate(SOURCES, NUM_TIME_POINTS, 20e-11, 20.9e-10, 0.005, 0.5, RANDOM_SEED);
    cout<<"MIN_I            = "<<W.MIN_I<<endl;
    cout<<"MAX_I            = "<<W.MAX_I<<endl;
    cout<<"MIN_TIME_STEP    = "<<W.MIN_TIME_STEP<<endl;
    cout<<"MAX_TIME_STEP    = "<<W.MAX_TIME_STEP<<endl;
    cout<<"NUM_TIME_POINTS  = "<<W.NUM_TIME_POINTS<<endl;
    W.cholmod_currents(currents, THREADS, NODES, csi);
    rc_simulation_stamp("Generating waveforms");
    
    // TR
    rc_tr * tr;
    if (TR == 1)
    {
        cout<<endl<<"Computing the exact voltage drop waveforms"<<endl;
        cout<<      "------------------------------------------"<<endl;
        tr = new rc_tr(g, W);
        tr->compute_voltage_drops_sm(W.time_points, W.current_vectors);
    }

    cholmod_dense * dc_envelopes;
    vector<cholmod_dense *> tran_envelopes;
    // Find Envelopes
    if (ENV_TYPE == 0)
    {
        cout<<endl<<"Finding the DC envelopes"<<endl;
        cout<<      "------------------------"<<endl;
        dc_envelopes = find_dc_envelope();
    
        if (TR == 1)
        {   
            // Compare
            vector<double> exact = tr->get_max_voltage_drops();
            double * bound = (double*)dc_envelopes->x;
            double max_error = 0, max_exact = 0;
            for (int i = 0; i<NODES; i++)
               if (max_exact < exact[i])
                    max_exact = exact[i];
            double scale = 55/max_exact; 
            for (int i = 0; i<NODES; i++)
            {
                exact[i] *= scale;
                bound[i] *= scale;
                if (max_error < abs(exact[i] - bound[i]))
                    max_error = abs(exact[i] - bound[i]);
            }

            cout<<"Max exact = "<<max_exact*scale<<endl;
            cout<<"Max error = "<<max_error<<endl;
        }
    }
    else if (ENV_TYPE == 1)
    {
        cout<<endl<<"Finding the TRAN envelopes"<<endl;
        cout<<      "--------------------------"<<endl;
        tran_envelopes = find_tran_envelope();
    }

    // Free Memory
    if (TR == 1)
    {   
        cholmod_finish(&tr->c);
        delete tr;
    }
    cholmod_free_dense(&dc_envelopes, &c);
    for (unsigned int i = 0; i<tran_envelopes.size(); i++)
        cholmod_free_dense(&tran_envelopes[i], &c);
    cholmod_free_sparse(&G, &c);
    cholmod_free_sparse(&C, &c);
    cholmod_free_sparse(&B, &c);
    
    cholmod_finish(&c);
    cholmod_finish(&g.c);
}

vector<cholmod_dense *> rc_simulation::find_tran_envelope()
{
    vector<cholmod_dense *> transient_envelopes;

    rc_simulation_stamp("Start");

    cout<<"Factorizing G..."<<endl;
    G_L = cholmod_analyze(G, &c);  // Symbolic factorization
    cholmod_factorize(G, G_L, &c); // Numeric factorization
    rc_simulation_stamp("Factorizing G");

    cout<<"Finding RC time step..."<<endl;
    find_rc_time_step();
    cout<<"RC time step found = "<<dtb<<" s"<<endl;
    rc_simulation_stamp("Finding RC time step");

    cout<<"Computing A..."<<endl;
    compute_A();
    rc_simulation_stamp("Computing A");

    cout<<"Factorizing A..."<<endl;
    A_L = cholmod_analyze(A, &c);  // Symbolic factorization
    cholmod_factorize(A, A_L, &c); // Numeric factorization
    rc_simulation_stamp("Factorizing A");

    int nv = NODES;
    int N  = NUM_TIME_POINTS;
    int TASK_LOAD = N/THREADS;
    vector<cholmod_dense* > w1(THREADS);    

    cout<<"Computing w1 - First set of linear systems..."<<endl;
    if (MULTITHREADING) 
    {
        #pragma omp parallel num_threads(THREADS)
        {
            int tid = omp_get_thread_num();
            w1[tid] = cholmod_solve(CHOLMOD_A, A_L, currents[tid], &c); 
        }
    }
    else
        for (int pid = 0; pid<THREADS; pid++)
            w1[pid] = cholmod_solve(CHOLMOD_A, A_L, currents[pid], &c); 
    
    rc_simulation_stamp("Computing w1 - First set of linear systems");
    
    for (int k = 0; k<THREADS; k++)
        cholmod_free_dense(&currents[k], &c);
    
    cout<<"Finding the RC history..."<<endl;
    find_rc_history(0.001);
//     psi = psi/4;
    cout<<"RC history found = "<<psi<<" s"<<endl;
    rc_simulation_stamp("Finding the RC history");

    vector<cholmod_dense *> et1_vec(THREADS);
    double * et1x;

    for (int i = 0; i<THREADS; i++)
        et1_vec[i] = cholmod_zeros(nv, TASK_LOAD, CHOLMOD_REAL, &c);
    
    transient_envelopes.resize(THREADS);
    
    if (MULTITHREADING) 
    {
        cout<<"Computing the et1_vec..."<<endl;
        #pragma omp parallel num_threads(THREADS)
        {
            double * w1x;
            int km1_u = 0, jtid;

            int tid = omp_get_thread_num();
            for (int k = tid*TASK_LOAD; k<(tid+1)*TASK_LOAD; k++)
            {
                et1x = (double*) et1_vec[tid]->x;

                if (k == 0)
                    km1_u = 0;
                else
                    km1_u = get_ku(k-1, psi);
//                 cout<<k - km1_u+1<<endl;
                int col = k%TASK_LOAD; 
                for (int j = k; j>0 && j>=km1_u; j--)
                {
                    jtid =  floor(1.0*j/TASK_LOAD);
                    w1x  = (double*) w1[jtid]->x;

                    // emax
                    for (int q = 0; q<nv; q++)
                    {
                        if (et1x[col*nv + q ] < w1x[(j - jtid*TASK_LOAD)*nv + q])
                            et1x[col*nv + q ] = w1x[(j - jtid*TASK_LOAD)*nv + q];
                    }
                }
            }
        }
        rc_simulation_stamp("Computing the et1_vec");
    
        for (int k = 0; k<THREADS; k++)
            cholmod_free_dense(&w1[k], &c);

        cholmod_free_factor(&A_L, &c);

        cout<<"Computing the final transient envelopes..."<<endl;
        #pragma omp parallel num_threads(THREADS)
        {
            int tid = omp_get_thread_num();
            cholmod_dense * A_et1 = smatrix_dmatrix_multiply(1, A, et1_vec[tid]);
            transient_envelopes[tid] = cholmod_solve(CHOLMOD_A, G_L, A_et1, &c);
            cholmod_free_dense(&A_et1, &c);
        }
        rc_simulation_stamp("Computing the final transient envelopes");
    }
    else
    {
        cout<<"Computing the et1_vec..."<<endl;
        for (int pid = 0; pid<THREADS; pid++)
        {
            double * w1x;
            int km1_u = 0, jpid;

            for (int k = pid*TASK_LOAD; k<(pid+1)*TASK_LOAD; k++)
            {
                et1x = (double*) et1_vec[pid]->x;

                if (k == 0)
                    km1_u = 0;
                else
                    km1_u = get_ku(k-1, psi);

//                 cout<<k - km1_u+1<<endl;
                int col = k%TASK_LOAD;
                for (int j = k; j>0 && j>=km1_u; j--)
                {
                    jpid =  floor(1.0*j/TASK_LOAD);
                    w1x  = (double*) w1[jpid]->x;

                    // emax
                    for (int q = 0; q<nv; q++)
                    {
                        if (et1x[col*nv + q ] < w1x[(j - jpid*TASK_LOAD)*nv + q])
                            et1x[col*nv + q ] = w1x[(j - jpid*TASK_LOAD)*nv + q];
                    }
                }
            }
        }
        rc_simulation_stamp("Computing the et1_vec");

        for (int k = 0; k<THREADS; k++)
            cholmod_free_dense(&w1[k], &c);
        
        cholmod_free_factor(&A_L, &c);

        cout<<"Computing the final transient envelopes..."<<endl;
        for (int pid = 0; pid<THREADS; pid++)
        {
            cholmod_dense * A_et1 = smatrix_dmatrix_multiply(1, A, et1_vec[pid]);
            transient_envelopes[pid] = cholmod_solve(CHOLMOD_A, G_L, A_et1, &c);
//             double * xx = (double*)transient_envelopes[pid]->x;
//             for (int i = 0; i<TASK_LOAD; i++)
//                 cout<<rc_simulation_time_points[i + pid*TASK_LOAD]<<" "<<xx[i*nv + 20000]<<endl;
            cholmod_free_dense(&A_et1, &c);
        }
        rc_simulation_stamp("Computing the final transient envelopes");
    }

    for (int k = 0; k<THREADS; k++)
        cholmod_free_dense(&et1_vec[k], &c);

    cholmod_free_sparse(&A, &c);
    cholmod_free_factor(&A_L, &c);
    cholmod_free_factor(&G_L, &c);
    
    printf("\nPrinting the rc_simulation profiler details [timer records]\
            \n-------------------------------------------------\n");
    printf("%60s %11s\n", "Task Description", "Wall (s)");
    printf("%60s %11s\n", "----------------", "--------");
    for (unsigned int i = 1; i<rc_simulation_stamps.size(); i++)
    {
        printf("%60s %11.2f\n", rc_simulation_stamps[i].second.c_str(), rc_simulation_stamps[i].first-rc_simulation_stamps[i-1].first);
    }
    return transient_envelopes;
}


cholmod_dense * rc_simulation::find_dc_envelope()
{
    rc_simulation_stamps.clear();
    rc_simulation_stamp("Start");
    cout<<"Factorizing G..."<<endl;
    G_L = cholmod_analyze(G, &c);  // Symbolic factorization
    cholmod_factorize(G, G_L, &c); // Numeric factorization
    rc_simulation_stamp("Factorizing G");  

    cout<<"Finding RC time step..."<<endl;
    find_rc_time_step();
    rc_simulation_stamp("Finding RC time step");  
    cout<<"RC time step found = "<<dtb<<" s"<<endl;

    cout<<"Computing A..."<<endl;
    compute_A();
    rc_simulation_stamp("Computing A");

    cout<<"Factorizing A..."<<endl;
    A_L = cholmod_analyze(A, &c);  // Symbolic factorization
    cholmod_factorize(A, A_L, &c); // Numeric factorization
    rc_simulation_stamp("Factorizing A");

    int nv = NODES;
    int N  = NUM_TIME_POINTS;
    int TASK_LOAD = N/THREADS;
    vector<cholmod_dense* > w(THREADS);    
    vector<double* > w_bars(THREADS);  // emax for each task
    for (int k = 0; k<THREADS; k++)
        w_bars[k]= (double*) calloc(nv, sizeof(double));

    cout<<"Computing w and w_bars - First set of linear systems..."<<endl;
    if (MULTITHREADING) 
    {
        #pragma omp parallel num_threads(THREADS)
        {
            int tid = omp_get_thread_num();
            w[tid] = cholmod_solve(CHOLMOD_A, A_L, currents[tid], &c); 

            double *wx = (double*)w[tid]->x;
 
            // emax
            for (int i = 0; i<nv; i++)
                for (int k = 0; k<TASK_LOAD; k++)
                    if (w_bars[tid][i]   < wx[nv*k+i])
                        w_bars[tid][i]   = wx[nv*k+i];
        }
    }
    else
    {
        for (int pid = 0; pid<THREADS; pid++)
        {   
            w[pid] = cholmod_solve(CHOLMOD_A, A_L, currents[pid], &c); 
        
            double *wx = (double*)w[pid]->x;
 
            // emax
            for (int i = 0; i<nv; i++)
                for (int k = 0; k<TASK_LOAD; k++)
                    if (w_bars[pid][i]   < wx[nv*k+i])
                        w_bars[pid][i]   = wx[nv*k+i];
        }
    }
    rc_simulation_stamp("Computing w and w_bars");

    for (int k = 0; k<THREADS; k++)
        cholmod_free_dense(&currents[k], &c);
 
    // Compute the final emax vector
    cout<<"Computing the final emax vector..."<<endl;
    double * W_bar = (double*) calloc(nv, sizeof(double));
    for (int i = 0; i<nv; i++)
        for (int k = 0; k<THREADS; k++)
            if (W_bar[i]   < w_bars[k][i])
                W_bar[i]   = w_bars[k][i];
    for (int k = 0; k<THREADS; k++)
    {
        free(w_bars[k]);
        cholmod_free_dense(&w[k], &c);
    }
    rc_simulation_stamp("Computing the final emax vector");

    //Compute Final bound
    cout<<"Computing the final bound..."<<endl;
    cholmod_dense * et1 = cholmod_zeros(nv, 1, CHOLMOD_REAL, &c);
    double * et1x = (double*) et1->x;
    for (int i = 0; i<nv; i++)
        et1x[i] = W_bar[i];
    
    free(W_bar);

    cholmod_dense * A_et1 = smatrix_dmatrix_multiply(1, A, et1);

    cholmod_dense * dc_envelopes = cholmod_solve(CHOLMOD_A, G_L, A_et1, &c);

    rc_simulation_stamp("Computing the final bound");
    cholmod_free_dense(&et1, &c);
    cholmod_free_dense(&A_et1, &c);
    cholmod_free_sparse(&A, &c);
    cholmod_free_factor(&A_L, &c);
    cholmod_free_factor(&G_L, &c);

//     printf("\nPrinting the profiler details [timer records]
//             \n---------------------------------------------\n");
    printf("%60s %11s\n", "Task Description", "Wall (s)");
    printf("%60s %11s\n", "----------------", "--------");
    for (unsigned int i = 1; i<rc_simulation_stamps.size(); i++)
    {
        printf("%60s %11.2f\n", rc_simulation_stamps[i].second.c_str(), rc_simulation_stamps[i].first-rc_simulation_stamps[i-1].first);
    }

    return dc_envelopes;
}

void rc_simulation::compute_A()
{
    // Find B
    double * Bx = (double*)B->x; 
    for (int i = 0; i<NODES; i++)
        Bx[i]/=dtb;

    // Find A
    A = smatrix_smatrix_add(1, G, 1, B);
}

void rc_simulation::find_rc_time_step()
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

void rc_simulation::find_rc_history(double eta)
{
    double * Cx = (double*)C->x;
    double cmax = Cx[0];
    double cmin = Cx[0];
    for (int i = 1; i<NODES; i++)
    {
        if (cmax < Cx[i])
            cmax = Cx[i];
        if (cmin > Cx[i])
            cmin = Cx[i];
    }
    
    double two_norm_ub = 0;
    for (int i = 0; i<NODES; i++)
        two_norm_ub += 0.1*0.1;
    two_norm_ub = sqrt(two_norm_ub);

    // return the history
    psi = dtb * log(sqrt(cmax/cmin)/(eta/two_norm_ub));
}

int rc_simulation::get_ku(int k, double psi)
{
    if (W.time_points[k] - psi <= 0)
        return 0;

    int ku = 0;
    for (int i = 0; i<=k; i++)
    {
        if (W.time_points[ku] > W.time_points[k] - psi)
            break;
        ku++;
    }
    return ku - 1;
}


cholmod_sparse * rc_simulation::smatrix_smatrix_add(double a1, cholmod_sparse *A1,
                                                    double a2, cholmod_sparse *A2)
{
    double alpha[1] = {a1}; 
    double beta [1] = {a2}; 
    cholmod_sparse * A3 = cholmod_add(A1, A2, alpha, beta, true, false, &c);
    return A3;
}

cholmod_dense  * rc_simulation::smatrix_dmatrix_multiply(double a1, cholmod_sparse *A1,
                                                         cholmod_dense  *A2)       
{
    cholmod_dense *A3 = cholmod_zeros(A1->nrow, A2->ncol, CHOLMOD_REAL, &c);
    double alpha[1] = {a1};
    double beta [1] = {0};
    cholmod_sdmult(A1, 0, alpha, beta, A2, A3, &c);
    return A3;
}

void rc_simulation_stamp(string op)
{
    pair<double, string> st;
    st.first  = omp_get_wtime();
    st.second = op;
    rc_simulation_stamps.push_back(st);
}
