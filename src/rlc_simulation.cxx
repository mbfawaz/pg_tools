#include <sstream>
#include <fstream>
#include "rlc_simulation.h"
#include "cholmod.h"
#include "string.h"
#include "omp.h"
#include "math.h"
#include "mosek_interface.h"
#include "rlc_grid.h"
#include "waveforms.h"
#include "rlc_tr.h"
#include <cmath>

using namespace std;

vector<pair<double, string> > rlc_simulation_stamps;
void rlc_simulation_stamp(std::string op);

rlc_simulation::rlc_simulation()
{
    cholmod_start(&c);
}

rlc_simulation::rlc_simulation(rlc_grid &g, string options_file)
{
    cholmod_start(&c);

    // Parameters of an RC grid
    // Read straight from the variable g
    G                   = g.G;
    C                   = g.C;
    L                   = g.L;
    M                   = g.M;
    csi                 = g.csi;
    NODES               = g.number_of_nodes;
    SOURCES             = g.number_of_sources;
    INDUCTORS           = g.number_of_inductors;

    // Loading options from the rc simulation options file
    ifstream infile;
    infile.open(options_file.c_str());
    if (!infile)
    {
        cout<<"ERROR: could not open rlc_simulation options file."<<endl;
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

    eye_nv = cholmod_speye(NODES,     NODES,     CHOLMOD_REAL, &c);
    eye_nl = cholmod_speye(INDUCTORS, INDUCTORS, CHOLMOD_REAL, &c);

    // Waveforms
    cout<<endl<<"Generating waveforms"<<endl;
    cout<<      "--------------------"<<endl;
    W.generate(SOURCES, NUM_TIME_POINTS, 20e-11, 20.9e-10, 0.005, 0.4, RANDOM_SEED);
    cout<<"MIN_I            = "<<W.MIN_I<<endl;
    cout<<"MAX_I            = "<<W.MAX_I<<endl;
    cout<<"MIN_TIME_STEP    = "<<W.MIN_TIME_STEP<<endl;
    cout<<"MAX_TIME_STEP    = "<<W.MAX_TIME_STEP<<endl;
    cout<<"NUM_TIME_POINTS  = "<<W.NUM_TIME_POINTS<<endl;
    W.cholmod_currents(currents, THREADS, NODES, csi);
    rlc_simulation_stamp("Generating waveforms");
    
    // TR
    rlc_tr * tr = new rlc_tr();
    if (TR == 1)
    {
        cout<<endl<<"Computing the exact voltage drop waveforms"<<endl;
        cout<<      "------------------------------------------"<<endl;
        tr = new rlc_tr(g, W);
        tr->compute_voltage_drops_sm(W.time_points, W.current_vectors);
    }

    // Find Envelopes
    if (ENV_TYPE == 0)
    {
        cout<<endl<<"Finding the DC envelopes"<<endl;
        cout<<      "------------------------"<<endl;
        find_dc_envelope();
        
        double * boundu = (double*)dc_envelopes_vu->x;
        double * boundl = (double*)dc_envelopes_vl->x;
        cout<<"Bound at 0 = "<<boundu[NODES-1]<<endl;
        cout<<"Bound at 0 = "<<boundl[NODES-1]<<endl;
    
        if (TR == 1)
        {   
            // Compare
            vector<double> exact = tr->get_max_voltage_drops();
            double * bound = (double*)dc_envelopes_vu->x;
            double max_error = 0, max_exact = 0;
            for (int i = 0; i<NODES; i++)
               if (max_exact < exact[i])
                    max_exact = exact[i];
            double scale = 1;
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

        if (TR == 1)
        {   
            // Compare
            vector<double> exact = tr->get_min_voltage_drops();
            double * bound = (double*)dc_envelopes_vl->x;
            double max_error = 0, min_exact = 0;
            for (int i = 0; i<NODES; i++)
               if (min_exact > exact[i])
                    min_exact = exact[i];
            double scale = 1; 
            for (int i = 0; i<NODES; i++)
            {
                exact[i] *= scale;
                bound[i] *= scale;
                if (max_error < abs(exact[i] - bound[i]))
                    max_error = abs(exact[i] - bound[i]);
            }

            cout<<"Min exact = "<<min_exact*scale<<endl;
            cout<<"Max error = "<<max_error<<endl;
        }
    }
    else if (ENV_TYPE == 1)
    {
        cout<<endl<<"Finding the TRAN envelopes"<<endl;
        cout<<      "--------------------------"<<endl;
        find_tran_envelope();
    }

    // Free Memory
//     if (TR == 1)
//     {   
//         cholmod_finish(&tr->c);
//         delete tr;
//     }
    cholmod_free_dense(&dc_envelopes_vu, &c);
    cholmod_free_dense(&dc_envelopes_iu, &c);
    cholmod_free_dense(&dc_envelopes_vl, &c);
    cholmod_free_dense(&dc_envelopes_il, &c);

    for (unsigned int i = 0; i<transient_envelopes_vu.size(); i++)
    {
        cholmod_free_dense(&transient_envelopes_vu[i], &c);
        cholmod_free_dense(&transient_envelopes_iu[i], &c);
        cholmod_free_dense(&transient_envelopes_vl[i], &c);
        cholmod_free_dense(&transient_envelopes_il[i], &c);
    }

    cholmod_free_sparse(&G, &c);
    cholmod_free_sparse(&C, &c);
    cholmod_free_sparse(&L, &c);
    cholmod_free_sparse(&M, &c);

    printf("\nPrinting the profiler details [timer records]\
            \n---------------------------------------------\n");
    printf("%60s %11s\n", "Task Description", "Wall (s)");
    printf("%60s %11s\n", "----------------", "--------");
    for (unsigned int i = 1; i<rlc_simulation_stamps.size(); i++)
    {
        printf("%60s %11.2f\n", rlc_simulation_stamps[i].second.c_str(), rlc_simulation_stamps[i].first-rlc_simulation_stamps[i-1].first);
    }
   
    cholmod_finish(&c);
    cholmod_finish(&g.c);
}

void rlc_simulation::find_tran_envelope()
{
    rlc_simulation_stamp("Start");

    cout<<"Finding RLC time step..."<<endl;
    find_rlc_time_step(-8.05, -8.05);
    cout<<"RLC time step found = "<<dtb<<" s"<<endl;
    rlc_simulation_stamp("Finding RLC time step");

    int nv = NODES;
    int nl = INDUCTORS;
    int N  = NUM_TIME_POINTS;
    int TASK_LOAD = N/THREADS;

    vector<cholmod_dense* > w1(THREADS), w2(THREADS);

    cholmod_sparse *Einv_Mt    = smatrix_smatrix_multiply(Einv, Mt);

    cout<<"Computing w1 and w2 - First set of linear systems..."<<endl;
    if (MULTITHREADING) 
    {
        #pragma omp parallel num_threads(THREADS)
        {
            int tid = omp_get_thread_num();
            w1[tid] = cholmod_solve(CHOLMOD_A, D_L, currents[tid], &c); 
            w2[tid] = smatrix_dmatrix_multiply(-1, Einv_Mt, w1[tid]); 
        }
    }
    else
    {
        for (int pid = 0; pid<THREADS; pid++)
        {
            w1[pid] = cholmod_solve(CHOLMOD_A, D_L, currents[pid], &c); 
            w2[pid] = smatrix_dmatrix_multiply(-1, Einv_Mt, w1[pid]); 
        }
    }
    rlc_simulation_stamp("Computing w1 - First set of linear systems");
    
    for (int k = 0; k<THREADS; k++)
        cholmod_free_dense(&currents[k], &c);

    // Find F22p and F22n
    cout<<"Finding F22p and F22n..."<<endl;
    F22p = smatrix_smatrix_add(0.5, F22,  0.5, abs_F22); cholmod_drop(1e-7, F22p, &c);
    F22m = smatrix_smatrix_add(0.5, F22, -0.5, abs_F22); cholmod_drop(1e-7, F22m, &c);
    rlc_simulation_stamp("Finding F22p and F22n");
    cholmod_free_sparse(&F22, &c);
    cholmod_free_sparse(&abs_F22, &c);

    cout<<"Finding the RLC history..."<<endl;
    find_rlc_history(0.001);
    cout<<"RLC history found = "<<psi<<" s"<<endl;
    rlc_simulation_stamp("Finding the RLC history");
    
    int NUM_BLOCKS = 12;
    int BLOCK_LOAD = N/NUM_BLOCKS;
    vector<cholmod_dense *> et1_vec(NUM_BLOCKS), et2_vec(NUM_BLOCKS), eb1_vec(NUM_BLOCKS), eb2_vec(NUM_BLOCKS);
    double * et1x, * et2x, * eb1x, * eb2x;

    for (int i = 0; i<NUM_BLOCKS; i++)
    {
        et1_vec[i] = cholmod_zeros(nv, BLOCK_LOAD, CHOLMOD_REAL, &c);
        et2_vec[i] = cholmod_zeros(nl, BLOCK_LOAD, CHOLMOD_REAL, &c);
        eb1_vec[i] = cholmod_zeros(nv, BLOCK_LOAD, CHOLMOD_REAL, &c);
        eb2_vec[i] = cholmod_zeros(nl, BLOCK_LOAD, CHOLMOD_REAL, &c);
    }
    transient_envelopes_vu.resize(NUM_BLOCKS);
    transient_envelopes_iu.resize(NUM_BLOCKS); 
    transient_envelopes_vl.resize(NUM_BLOCKS);
    transient_envelopes_il.resize(NUM_BLOCKS); 
    
    if (MULTITHREADING) 
    {
        cout<<"Computing et1_vec, et2_vec, eb1_vec, and eb2_vec..."<<endl;
        #pragma omp parallel num_threads(THREADS)
        {
            double * w1x, * w2x;
            int km1_u = 0, jtid;

            int tid = omp_get_thread_num();
            for (int k = tid*TASK_LOAD; k<(tid+1)*TASK_LOAD; k++)
            {
                et1x = (double*) et1_vec[k%NUM_BLOCKS]->x;
                et2x = (double*) et2_vec[k%NUM_BLOCKS]->x;
                eb1x = (double*) eb1_vec[k%NUM_BLOCKS]->x;
                eb2x = (double*) eb2_vec[k%NUM_BLOCKS]->x;

                if (k == 0)
                    km1_u = 0;
                else
                    km1_u = get_ku(k-1, psi);
                int col = k%BLOCK_LOAD;
                for (int j = k; j>0 && j>=km1_u; j--)
                {
                    jtid =  floor(1.0*j/TASK_LOAD);
                    w1x  = (double*) w1[jtid]->x;
                    w2x  = (double*) w2[jtid]->x;
                    // top part
                    // emax
                    for (int q = 0; q<nv; q++)
                    {
                        if (et1x[col*nv + q ] < w1x[(j - jtid*TASK_LOAD)*nv + q])
                            et1x[col*nv + q ] = w1x[(j - jtid*TASK_LOAD)*nv + q];
                    }
                    //emin
                    for (int q = 0; q<nv; q++)
                    {
                        if (eb1x[col*nv + q ] > w1x[(j - jtid*TASK_LOAD)*nv + q])
                            eb1x[col*nv + q ] = w1x[(j - jtid*TASK_LOAD)*nv + q];
                    }
                    // bottom part
                    // emax
                    for (int q = 0; q<nl; q++)
                    {
                        if (et2x[col*nl + q ] < w2x[(j - jtid*TASK_LOAD)*nl + q])
                            et2x[col*nl + q ] = w2x[(j - jtid*TASK_LOAD)*nl + q];
                    }
                    //emin
                    for (int q = 0; q<nl; q++)
                    {
                        if (eb2x[col*nl + q ] > w2x[(j - jtid*TASK_LOAD)*nl + q])
                            eb2x[col*nl + q ] = w2x[(j - jtid*TASK_LOAD)*nl + q];
                    }
                }
            }
        }
        rlc_simulation_stamp("Computing et1_vec, et2_vec, eb1_vec, and eb2_vec");

        cout<<"Computing the final transient envelopes..."<<endl;
        #pragma omp parallel num_threads(THREADS)
        {
            int tid = omp_get_thread_num();
            int blocks_per_thread = NUM_BLOCKS/THREADS; // assuming integer value
            for (int bid = tid*blocks_per_thread; bid<(tid+1)*blocks_per_thread; bid++)
            {
                if (bid == tid*blocks_per_thread)
                    solve_eye_minus_F_tilde_eq_f
                                (  et1_vec[bid]  ,  et2_vec[bid]  ,  eb1_vec[bid]  ,  eb2_vec[bid],
                                  &transient_envelopes_vu[bid]  , &transient_envelopes_iu[bid]  , &transient_envelopes_vl[bid]  , &transient_envelopes_il[bid],
                                   et1_vec[bid]  ,  et2_vec[bid]  ,  eb1_vec[bid]  ,  eb2_vec[bid]);
                else
                    solve_eye_minus_F_tilde_eq_f
                                (  et1_vec[bid]    ,  et2_vec[bid]    ,  eb1_vec[bid]    ,  eb2_vec[bid],
                                  &transient_envelopes_vu[bid]    , &transient_envelopes_iu[bid]    , &transient_envelopes_vl[bid]    , &transient_envelopes_il[bid],
                                   transient_envelopes_vu[bid-1]  ,  transient_envelopes_iu[bid-1]  ,  transient_envelopes_vl[bid-1]  ,  transient_envelopes_il[bid-1]);
            }
        }
        rlc_simulation_stamp("Computing the final transient envelopes");
    }
    else
    {
        double * w1x, * w2x;
        int km1_u = 0, jpid;
        cout<<"Computing et1_vec, et2_vec, eb1_vec, and eb2_vec..."<<endl;
        for (int pid = 0; pid<THREADS; pid++)
        {
            for (int k = pid*TASK_LOAD; k<(pid+1)*TASK_LOAD; k++)
            {
//                 cout<<k<<endl;
                et1x = (double*) et1_vec[floor(k/BLOCK_LOAD)]->x;
                et2x = (double*) et2_vec[floor(k/BLOCK_LOAD)]->x;
                eb1x = (double*) eb1_vec[floor(k/BLOCK_LOAD)]->x;
                eb2x = (double*) eb2_vec[floor(k/BLOCK_LOAD)]->x;

                if (k == 0)
                    km1_u = 0;
                else
                    km1_u = get_ku(k-1, psi);
                int col = k%BLOCK_LOAD;
                for (int j = k; j>0 && j>=km1_u; j--)
                {
                    jpid =  floor(1.0*j/TASK_LOAD);
                    w1x  = (double*) w1[jpid]->x;
                    w2x  = (double*) w2[jpid]->x;
                    // top part
                    // emax
                    for (int q = 0; q<nv; q++)
                    {
                        if (et1x[col*nv + q ] < w1x[(j - jpid*TASK_LOAD)*nv + q])
                            et1x[col*nv + q ] = w1x[(j - jpid*TASK_LOAD)*nv + q];
//                         if (k%NUM_BLOCKS == 0)
//                             cout<<"col = "<<col<<" and: "<<"At "<<col*nv+q<<" "<<et1x[col*nv+q]<<endl;
                    }
                    //emin
                    for (int q = 0; q<nv; q++)
                    {
                        if (eb1x[col*nv + q ] > w1x[(j - jpid*TASK_LOAD)*nv + q])
                            eb1x[col*nv + q ] = w1x[(j - jpid*TASK_LOAD)*nv + q];
                    }
                    // bottom part
                    // emax
                    for (int q = 0; q<nl; q++)
                    {
                        if (et2x[col*nl + q ] < w2x[(j - jpid*TASK_LOAD)*nl + q])
                            et2x[col*nl + q ] = w2x[(j - jpid*TASK_LOAD)*nl + q];
                    }
                    //emin
                    for (int q = 0; q<nl; q++)
                    {
                        if (eb2x[col*nl + q ] > w2x[(j - jpid*TASK_LOAD)*nl + q])
                            eb2x[col*nl + q ] = w2x[(j - jpid*TASK_LOAD)*nl + q];
                    }
                }
            }
        }
        rlc_simulation_stamp("Computing et1_vec, et2_vec, eb1_vec, and eb2_vec");
        cout<<"Computing the final transient envelopes..."<<endl;
        for (int bid = 0; bid<NUM_BLOCKS; bid++)
        {
            if (bid <= 1)
                solve_eye_minus_F_tilde_eq_f
                                (  et1_vec[bid]  ,  et2_vec[bid]  ,  eb1_vec[bid]  ,  eb2_vec[bid],
                                  &transient_envelopes_vu[bid]  , &transient_envelopes_iu[bid]  , &transient_envelopes_vl[bid]  , &transient_envelopes_il[bid],
                                   et1_vec[bid]  ,  et2_vec[bid]  ,  eb1_vec[bid]  ,  eb2_vec[bid]);
            else
                solve_eye_minus_F_tilde_eq_f
                                (  et1_vec[bid]    ,  et2_vec[bid]    ,  eb1_vec[bid]    ,  eb2_vec[bid],
                                  &transient_envelopes_vu[bid]    , &transient_envelopes_iu[bid]    , &transient_envelopes_vl[bid]    , &transient_envelopes_il[bid],
                                   transient_envelopes_vu[bid-1]  ,  transient_envelopes_iu[bid-1]  ,  transient_envelopes_vl[bid-1]  ,  transient_envelopes_il[bid-1]);
        }
        rlc_simulation_stamp("Computing the final transient envelopes");
    }
    int TEST = NODES-1;
    vector<double> rv;
//     ss.str(string());
//     ss<<"mat_data_"<<NODES<<".m";

    FILE * mat_file = fopen("mat_data.m", "wb");

//     fprintf(mat_file, "vlb = [");
//     for (unsigned int i = 0; i<vlb.size(); i++)
//         fprintf(mat_file, "%.10f\n", vlb[i]);
//     fprintf(mat_file, "];\n\n");

 
    fprintf(mat_file, "time_points = [");
    for (unsigned int i = 0; i<NUM_TIME_POINTS; i++)
        fprintf(mat_file, "%.10f\n", W.time_points[i]);
    fprintf(mat_file, "];\n\n");



    for (int i = 0; i<transient_envelopes_vu.size(); i++)
    {
        double * x = (double*) transient_envelopes_vu[i]->x;

        for (int j = 0; j<NUM_TIME_POINTS/transient_envelopes_vu.size(); j++)
            rv.push_back(x[j*NODES + TEST]);
    }
    
    fprintf(mat_file, "vub = [");
    for (unsigned int i = 0; i<NUM_TIME_POINTS; i++)
        fprintf(mat_file, "%.10f\n", rv[i]);
    fprintf(mat_file, "];\n\n");

    rv.clear();
    for (int i = 0; i<transient_envelopes_vl.size(); i++)
    {
        double * x = (double*) transient_envelopes_vl[i]->x;

        for (int j = 0; j<NUM_TIME_POINTS/transient_envelopes_vl.size(); j++)
            rv.push_back(x[j*NODES + TEST]);
    }
    fprintf(mat_file, "vlb = [");
    for (unsigned int i = 0; i<NUM_TIME_POINTS; i++)
        fprintf(mat_file, "%.10f\n", rv[i]);
    fprintf(mat_file, "];\n\n");

    fclose(mat_file);

    cout<<"Done printing..."<<endl;
    

    for (int k = 0; k<NUM_BLOCKS; k++)
    {
        cholmod_free_dense(&et1_vec[k], &c);
        cholmod_free_dense(&et2_vec[k], &c);
        cholmod_free_dense(&eb1_vec[k], &c);
        cholmod_free_dense(&eb2_vec[k], &c);
    }
    for (int k = 0; k<THREADS; k++)
    {
        cholmod_free_dense(&w1[k], &c);
        cholmod_free_dense(&w2[k], &c);
    }

    for (int k = 0; k<THREADS; k++)
        cholmod_free_dense(&et1_vec[k], &c);

    cholmod_free_factor(&D_L, &c);


/*    rlc_simulation_stamp("Start");

    cout<<"Factorizing G..."<<endl;
    G_L = cholmod_analyze(G, &c);  // Symbolic factorization
    cholmod_factorize(G, G_L, &c); // Numeric factorization
    rlc_simulation_stamp("Factorizing G");

    cout<<"Finding RC time step..."<<endl;
    find_rlc_time_step();
    cout<<"RC time step found = "<<dtb<<" s"<<endl;
    rlc_simulation_stamp("Finding RC time step");

    cout<<"Computing A..."<<endl;
    compute_A();
    rlc_simulation_stamp("Computing A");

    cout<<"Factorizing A..."<<endl;
    A_L = cholmod_analyze(A, &c);  // Symbolic factorization
    cholmod_factorize(A, A_L, &c); // Numeric factorization
    rlc_simulation_stamp("Factorizing A");

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
    
    rlc_simulation_stamp("Computing w1 - First set of linear systems");
    
    for (int k = 0; k<THREADS; k++)
        cholmod_free_dense(&currents[k], &c);
    
    cout<<"Finding the RC history..."<<endl;
    find_rlc_history(0.001);
//     psi = psi/4;
    cout<<"RC history found = "<<psi<<" s"<<endl;
    rlc_simulation_stamp("Finding the RC history");

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
        rlc_simulation_stamp("Computing the et1_vec");
    
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
        rlc_simulation_stamp("Computing the final transient envelopes");
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
        rlc_simulation_stamp("Computing the et1_vec");

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
//                 cout<<rlc_simulation_time_points[i + pid*TASK_LOAD]<<" "<<xx[i*nv + 20000]<<endl;
            cholmod_free_dense(&A_et1, &c);
        }
        rlc_simulation_stamp("Computing the final transient envelopes");
    }

    for (int k = 0; k<THREADS; k++)
        cholmod_free_dense(&et1_vec[k], &c);

    cholmod_free_sparse(&A, &c);
    cholmod_free_factor(&A_L, &c);
    cholmod_free_factor(&G_L, &c);
    
    printf("\nPrinting the rlc_simulation profiler details [timer records]\
            \n-------------------------------------------------\n");
    printf("%60s %11s\n", "Task Description", "Wall (s)");
    printf("%60s %11s\n", "----------------", "--------");
    for (unsigned int i = 1; i<rlc_simulation_stamps.size(); i++)
    {
        printf("%60s %11.2f\n", rlc_simulation_stamps[i].second.c_str(), rlc_simulation_stamps[i].first-rlc_simulation_stamps[i-1].first);
    }
    return transient_envelopes;
    */
}


void rlc_simulation::find_dc_envelope()
{
    rlc_simulation_stamps.clear();
    rlc_simulation_stamp("Start");
    
    cout<<"Finding RLC time step..."<<endl;
    find_rlc_time_step(-8.2, -8.2);
    rlc_simulation_stamp("Finding RLC time step");  
    cout<<"RLC time step found = "<<dtb<<" s"<<endl;

    int nv = NODES;
    int nl = INDUCTORS;
    int n = nv + nl;
    int N  = NUM_TIME_POINTS;
    int TASK_LOAD = N/THREADS;



    vector<cholmod_dense* > w1(THREADS), w2(THREADS);    
    vector<double* > w_bars(THREADS);  // eopt for each task
    for (int k = 0; k<THREADS; k++)
        w_bars[k]= (double*) calloc(2*n, sizeof(double));

    cholmod_sparse *Einv_Mt    = smatrix_smatrix_multiply(Einv, Mt);

    cout<<"Computing w1, w2, and w_bars - First set of linear systems..."<<endl;
    if (MULTITHREADING) 
    {
        #pragma omp parallel num_threads(THREADS)
        {
            int tid = omp_get_thread_num();
            w1[tid] = cholmod_solve(CHOLMOD_A, D_L, currents[tid], &c); 
            w2[tid] = smatrix_dmatrix_multiply(-1, Einv_Mt, w1[tid]); 

            double *w1x = (double*)w1[tid]->x;
            double *w2x = (double*)w2[tid]->x; 
 
            // top part
            for (int i = 0; i<nv; i++)
            {
                // emax
                for (int k = 0; k<TASK_LOAD; k++)
                    if (w_bars[tid][i]   < w1x[nv*k+i])
                        w_bars[tid][i]   = w1x[nv*k+i];

                // emin
                for (int k = 0; k<TASK_LOAD; k++)
                    if (w_bars[tid][n+i] > w1x[nv*k+i])
                        w_bars[tid][n+i] = w1x[nv*k+i];
            }

            // bottom part
            for (int j = 0; j<nl; j++)
            {
                // emax
                for (int k = 0; k<TASK_LOAD; k++)
                    if (w_bars[tid][nv+j]   < w2x[nl*k+j])
                        w_bars[tid][nv+j]   = w2x[nl*k+j];

                // emin
                for (int k = 0; k<TASK_LOAD; k++)
                    if (w_bars[tid][n+nv+j] > w2x[nl*k+j])
                        w_bars[tid][n+nv+j] = w2x[nl*k+j];
            }
        }
    }
    else
    {
        for (int pid = 0; pid<THREADS; pid++)
        {   
            w1[pid] = cholmod_solve(CHOLMOD_A, D_L, currents[pid], &c); 
            w2[pid] = smatrix_dmatrix_multiply(-1, Einv_Mt, w1[pid]); 

            double *w1x = (double*)w1[pid]->x;
            double *w2x = (double*)w2[pid]->x; 
 
            // top part
            for (int i = 0; i<nv; i++)
            {
                // emax
                for (int k = 0; k<TASK_LOAD; k++)
                    if (w_bars[pid][i]   < w1x[nv*k+i])
                        w_bars[pid][i]   = w1x[nv*k+i];

                // emin
                for (int k = 0; k<TASK_LOAD; k++)
                    if (w_bars[pid][n+i] > w1x[nv*k+i])
                        w_bars[pid][n+i] = w1x[nv*k+i];
            }

            // bottom part
            for (int j = 0; j<nl; j++)
            {
                // emax
                for (int k = 0; k<TASK_LOAD; k++)
                    if (w_bars[pid][nv+j]   < w2x[nl*k+j])
                        w_bars[pid][nv+j]   = w2x[nl*k+j];

                // emin
                for (int k = 0; k<TASK_LOAD; k++)
                    if (w_bars[pid][n+nv+j] > w2x[nl*k+j])
                        w_bars[pid][n+nv+j] = w2x[nl*k+j];
            }
        }
    }
    rlc_simulation_stamp("Computing w1, w2, and w_bars - First set of linear systems");
    for (int k = 0; k<THREADS; k++)
        cholmod_free_dense(&currents[k], &c);

    // Compute the final eopt vector
    cout<<"Computing the final eopt vector..."<<endl;
    double * W_bar = (double*) calloc(2*n, sizeof(double));
    for (int i = 0; i<n; i++)
        for (int k = 0; k<THREADS; k++)
            if (W_bar[i]   < w_bars[k][i])
                W_bar[i]   = w_bars[k][i];
    for (int i = 0; i<n; i++)
        for (int k = 0; k<THREADS; k++)
            if (W_bar[n+i] > w_bars[k][n+i])
                W_bar[n+i] = w_bars[k][n+i];
    for (int k = 0; k<THREADS; k++)
        free(w_bars[k]);
    rlc_simulation_stamp("Computing the final eopt vector");

    // Find F22p and F22n
    cout<<"Finding F22p and F22n..."<<endl;
    F22p = smatrix_smatrix_add(0.5, F22,  0.5, abs_F22); cholmod_drop(1e-7, F22p, &c);
    F22m = smatrix_smatrix_add(0.5, F22, -0.5, abs_F22); cholmod_drop(1e-7, F22m, &c);
    rlc_simulation_stamp("Finding F22p and F22n");

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
        et1x[i] = W_bar[i];
        eb1x[i] = W_bar[n+i];
    }   
    for (int i = 0; i<nl; i++)
    {       
        et2x[i] = W_bar[nv+i];
        eb2x[i] = W_bar[n+nv+i];
    }   
    free(W_bar); 

    solve_eye_minus_F_tilde_eq_f( et1,  et2,  eb1,  eb2, 
                                 &dc_envelopes_vu, &dc_envelopes_iu, &dc_envelopes_vl, &dc_envelopes_il,
                                  et1,  et2,  eb1,  eb2);

    rlc_simulation_stamp("Computing the final bound");

    cholmod_free_dense(&et1, &c);
    cholmod_free_dense(&et2, &c);
    cholmod_free_dense(&eb1, &c);
    cholmod_free_dense(&eb2, &c);
}

void rlc_simulation::compute_D()
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
void rlc_simulation::solve_eye_minus_F_tilde_eq_f
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
    
        if (1.0*residue/sum_norm_of_f < 1e-6)
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

void  rlc_simulation::F_tilde_times_x(cholmod_dense *  xt1, cholmod_dense *  xt2, cholmod_dense *  xb1, cholmod_dense *  xb2,
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




void rlc_simulation::find_rlc_time_step(double log_dtb_min, double log_dtb_max)
{
    // bisection method
    srand(NODES);
    double a1 = log_dtb_min;
    double b1 = log_dtb_max;
    double p1 = a1 + (b1 - a1)*rand()/RAND_MAX;
    
    dtb = pow(10, p1);
    while (true)
    {
        if (is_spec_radius_less_than(1) == 0)
        {   
            a1 = p1;
        }
        else
            return;
        p1 = (a1 + b1)/2;
        dtb = pow(10, p1);
    }
}

int rlc_simulation::is_spec_radius_less_than(double target)
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
        cout<<lambda_min<<" to "<<lambda_max<<endl;
 
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

void  rlc_simulation::abs_F_times_x(cholmod_dense *  x1, cholmod_dense *  x2,
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


void rlc_simulation::find_rlc_history(double eta)
{
    
    // Finding the history 
    psi = 1e-8/4;
    return;

    int nv = NODES;
    G = smatrix_smatrix_add(1, G, 0.1, eye_nv);
    G_L = cholmod_analyze(G, &c);  // Symbolic factorization of G
    cholmod_factorize(G, G_L, &c); // Numeric factorization of G
    
    // Compute K
    cholmod_sparse * Linv = cholmod_copy_sparse(L, &c);
    double * Linv_x = (double*) Linv -> x;
    for (unsigned int i = 0; i<Linv->ncol; i++)
        Linv_x[i] = 1.0/Linv_x[i];
    cholmod_sparse * M_Linv = smatrix_smatrix_multiply(M, Linv);
    cholmod_sparse * K = smatrix_smatrix_multiply(M_Linv, Mt);
    K = smatrix_smatrix_add(1, K, 100, eye_nv);
    cholmod_factor * K_L;
    K_L = cholmod_analyze(K, &c);  // Symbolic factorization of K
    cholmod_factorize(K, K_L, &c); // Numeric factorization of K

    // Step 1 Finding alpha^*. For that, find lambda_1 using the power method
    double alpha_star;
    double l = 0;  int idx = 0;
    cholmod_dense * xk = cholmod_ones(nv, 1, CHOLMOD_REAL, &c);
    cholmod_dense * xkp1 = cholmod_ones(nv, 1, CHOLMOD_REAL, &c);
    while (true)
    {
        cholmod_dense * Ginv_xk = cholmod_solve(CHOLMOD_A, G_L, xk, &c); 
        cholmod_dense * CGinv_xk = smatrix_dmatrix_multiply(1, C, Ginv_xk);

         
        cholmod_dense * Kinv_xk = cholmod_solve(CHOLMOD_A, K_L, xk, &c);
        cholmod_dense * G_Kinv_xk = smatrix_dmatrix_multiply(0.25, G, Kinv_xk);

       
        double * xkp1_x = (double*) xkp1->x;
        double * CGinv_xk_x = (double*) CGinv_xk->x;
        double * G_Kinv_xk_x = (double*) G_Kinv_xk->x;
        for (int i = 0; i<nv; i++)
            xkp1_x[i] = CGinv_xk_x[i] + G_Kinv_xk_x[i];

//         c.print = 5;
//         cholmod_print_dense(xkp1, "this", &c);

        double num = 0, den = 0;
        double * xk_x = (double*) xk->x;
        for (int i =0; i<nv; i++)
        {
            num += xkp1_x[i]*xk_x[i];
            den += xk_x  [i]*xk_x[i];
        }
        l = num/den;

        cout<<l<<endl;
        
        for (int i = 0; i<nv; i++)
            xk_x[i] = xkp1_x[i];
        idx++;
        if (idx == 10)
        {
            alpha_star = 1.0/l;
            break;
        }
    }
    double alpha = 0.5*alpha_star;
    cout<<"Alpha chosen = "<<alpha<<endl;

    // Step 2: Compute E_tilde and its Cholesky
    cholmod_sparse * alpha_C = cholmod_copy_sparse(C, &c);
    double * alpha_Cx = (double*) alpha_C->x;
    for (int i = 0; i<alpha_C->nzmax; i++)
        alpha_Cx[i] *= alpha;
    
    cholmod_sparse * row1_of_E_tilde = cholmod_horzcat(K,       alpha_C, true, &c);
    cholmod_sparse * row2_of_E_tilde = cholmod_horzcat(alpha_C, C      , true, &c);
    cholmod_sparse * E_tilde         = cholmod_vertcat(row1_of_E_tilde, row2_of_E_tilde, true, &c);
    cholmod_factor * E_tilde_L;
    E_tilde_L = cholmod_analyze(E_tilde, &c);  // Symbolic factorization of E_tilde
    cholmod_factorize(E_tilde, E_tilde_L, &c); // Numeric factorization of E_tilde


    // Step 3: Find the dominant eigenvalue of E_tilde
    xk = cholmod_ones(2*nv, 1, CHOLMOD_REAL, &c);
    xkp1 = cholmod_ones(2*nv, 1, CHOLMOD_REAL, &c);
    double e_tilde_max;
    idx = 1;
    while (true)
    {
        cholmod_dense * xkp1 = smatrix_dmatrix_multiply(1, E_tilde, xk); 
       
        double * xkp1_x = (double*) xkp1->x;

        double num = 0, den = 0;
        double * xk_x = (double*) xk->x;
        for (int i =0; i<2*nv; i++)
        {
            num += xkp1_x[i]*xk_x[i];
            den += xk_x  [i]*xk_x[i];
        }
        l = num/den;

        cout<<l<<endl;
        
        for (int i = 0; i<2*nv; i++)
            xk_x[i] = xkp1_x[i];
        idx++;
        if (idx == 10)
        {
            e_tilde_max = l;
            break;
        }
    }

    // Step 3: Find the dominant eigenvalue of inv(E_tilde)
    xk = cholmod_ones(2*nv, 1 , CHOLMOD_REAL, &c);
    xkp1 = cholmod_ones(2*nv, 1, CHOLMOD_REAL, &c);
    double e_tilde_min;
    idx = 1;
    while (true)
    {
        xkp1 = cholmod_solve(CHOLMOD_A, E_tilde_L, xk, &c); 

        double * xkp1_x = (double*) xkp1->x;

        double num = 0, den = 0;
        double * xk_x = (double*) xk->x;
        for (int i =0; i<2*nv; i++)
        {
            num += xkp1_x[i]*xk_x[i];
            den += xk_x  [i]*xk_x[i];
        }
        l = num/den;

        cout<<1.0/l<<endl;
        
        for (int i = 0; i<2*nv; i++)
            xk_x[i] = xkp1_x[i];
        idx++;
        if (idx == 3)
        {
            e_tilde_min = 1.0/l;
            break;
        }
    }

    // Step 4: Build A_tilde
    cholmod_sparse * minus_alpha_K = cholmod_copy_sparse(K, &c);
    double * minus_alpha_Kx = (double*)minus_alpha_K->x;
    for (int i = 0; i<minus_alpha_K->nzmax; i++)
        minus_alpha_Kx[i] *= -1.0*alpha;
    
    cholmod_sparse * K_minus_alpha_G = smatrix_smatrix_add(1, K, -1.0*alpha, G);
   
    cholmod_sparse * minus_K = cholmod_copy_sparse(K, &c);
    double * minus_Kx = (double*)minus_K->x;
    for (int i = 0; i<minus_K->nzmax; i++)
        minus_Kx[i] *= -1.0;
    
    cholmod_sparse * minus_G_plus_alpha_C = smatrix_smatrix_add(-1.0, G, alpha, C);
 
    cholmod_sparse * row1_of_A_tilde = cholmod_horzcat(minus_alpha_K, K_minus_alpha_G     , true, &c);
    cholmod_sparse * row2_of_A_tilde = cholmod_horzcat(minus_K      , minus_G_plus_alpha_C, true, &c);
    cholmod_sparse * A_tilde         = cholmod_vertcat(row1_of_A_tilde, row2_of_A_tilde, true, &c);
   
    // Step 5: Find lambda_min
    cholmod_sparse * A_tilde_t = transpose(A_tilde);
    cholmod_sparse * Asym = smatrix_smatrix_add(0.5, A_tilde, 0.5, A_tilde_t);
    cholmod_factor * Asym_L; 
    Asym_L = cholmod_analyze(Asym, &c);  // Symbolic factorization of Asym 
    cholmod_factorize(Asym, Asym_L, &c); // Numeric factorization of Asym


    xk = cholmod_ones(2*nv, 1 , CHOLMOD_REAL, &c);
    xkp1 = cholmod_ones(2*nv, 1, CHOLMOD_REAL, &c);
    double lmin;
    idx = 1;
    while (true)
    {
        xkp1 = cholmod_solve(CHOLMOD_Lt, E_tilde_L, xk, &c); 
        xkp1 = smatrix_dmatrix_multiply(1, E_tilde, xkp1);
        xkp1 = cholmod_solve(CHOLMOD_A, Asym_L, xkp1, &c);
        xkp1 = smatrix_dmatrix_multiply(-1, E_tilde, xkp1);
        xkp1 = cholmod_solve(CHOLMOD_L, E_tilde_L, xkp1, &c);

        double * xkp1_x = (double*) xkp1->x;

        double num = 0, den = 0;
        double * xk_x = (double*) xk->x;
        for (int i =0; i<2*nv; i++)
        {
            num += xkp1_x[i]*xk_x[i];
            den += xk_x  [i]*xk_x[i];
        }

        l = num/den;

        cout<<l<<endl;
        
        for (int i = 0; i<2*nv; i++)
            xk_x[i] = xkp1_x[i];
        idx++;
        if (idx == 10)
        {
            lmin = 1/l;
            break;
        }
    }
    e_tilde_min *= 10; 

    psi = (1/lmin)*log(sqrt(e_tilde_max/e_tilde_min)/(0.001/1));
    cout<<"history = "<<psi<<endl;

    exit(0);
}

int rlc_simulation::get_ku(int k, double psi)
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

cholmod_sparse * rlc_simulation::smatrix_smatrix_multiply(cholmod_sparse *A1, cholmod_sparse *A2)
{
    cholmod_sparse * A3 = cholmod_ssmult(A1, A2, 0, 1, 1, &c);
    return A3;
}

cholmod_sparse * rlc_simulation::transpose(cholmod_sparse * A1)
{
    cholmod_sparse * A1t = cholmod_allocate_sparse(A1->ncol, A1->nrow, A1->nzmax, 0, 1, 0, CHOLMOD_REAL, &c);
    cholmod_transpose_unsym(A1, 1, NULL, NULL, 0, A1t, &c);
    return A1t;
}

cholmod_sparse * rlc_simulation::smatrix_smatrix_add(double a1, cholmod_sparse *A1,
                                                     double a2, cholmod_sparse *A2)
{
    double alpha[1] = {a1}; 
    double beta [1] = {a2}; 
    cholmod_sparse * A3 = cholmod_add(A1, A2, alpha, beta, true, false, &c);
    return A3;
}

cholmod_dense  * rlc_simulation::smatrix_dmatrix_multiply(double a1, cholmod_sparse *A1,
                                                         cholmod_dense  *A2)       
{
    cholmod_dense *A3 = cholmod_zeros(A1->nrow, A2->ncol, CHOLMOD_REAL, &c);
    double alpha[1] = {a1};
    double beta [1] = {0};
    cholmod_sdmult(A1, 0, alpha, beta, A2, A3, &c);
    return A3;
}

void rlc_simulation_stamp(string op)
{
    pair<double, string> st;
    st.first  = omp_get_wtime();
    st.second = op;
    rlc_simulation_stamps.push_back(st);
}
