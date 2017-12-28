#include <sstream>
#include <fstream>
#include "rc_tr.h"
#include "cholmod.h"
#include "string.h"
#include "omp.h"
#include "math.h"
#include "mosek_interface.h"
#include "rc_grid.h"
#include "waveforms.h"
#include <cmath>
#include "assert.h"

using namespace std;

bool SAVE_MEMORY = 1;

vector<pair<double, string> > rc_tr_stamps;
void rc_tr_stamp(std::string op);

rc_tr::rc_tr()
{
    cholmod_start(&c);
}

rc_tr::rc_tr(rc_grid &g, waveforms &W)
{
    cholmod_start(&c);

    // Parameters of an RC grid
    // Read straight from the variable g
    G                   = g.G;
    C                   = g.C;
    csi                 = g.csi;
    NODES               = g.number_of_nodes;
    SOURCES             = g.number_of_sources;
    TIME_STEP           = W.MIN_TIME_STEP;

    if (!SAVE_MEMORY)
        get_current_points(W.time_points, W.current_vectors);
}

void rc_tr::compute_voltage_drops_sm(vector<double> & time_points, vector<vector<double> > & current_vectors)
{   
    rc_tr_stamp("Start");
    // Find B_TR 
    cout<<"Finding B_TR..."<<endl;
    B_TR = smatrix_smatrix_add(1.0/TIME_STEP, C, -0.5, G); 
    rc_tr_stamp("Finding B_TR");

    // Find A_TR
    cout<<"Finding A_TR..."<<endl;
    A_TR = smatrix_smatrix_add(0.5, G, 1.0/TIME_STEP, C);
    rc_tr_stamp("Finding A_TR");

    // Factorize A_TR
    cout<<"Factorizing A_TR_L..."<<endl;
    A_TR_L = cholmod_analyze(A_TR, &c);     // symbolic factorization
    cholmod_factorize(A_TR, A_TR_L, &c);    // numeric factorization
    rc_tr_stamp("Factorizing A_TR_L");

    // Factorize G
    cout<<"Factorizing G..."<<endl;
    G_L = cholmod_analyze(G, &c);           // symbolic factorization
    cholmod_factorize(G, G_L, &c);          // numeric factorization
    rc_tr_stamp("Factorizing G");

    // DC operating point 
    cout<<"Performing simulation..."<<endl;
    sim_time = time_points[time_points.size()-1];
    cout<<"Sim time = "<<sim_time<<endl;
    int time_index = 0;
    cholmod_dense * prev_cv;
    cholmod_dense * current_cv;
    cholmod_dense * prev_vd;
    cholmod_dense * current_vd;
    prev_cv    = cholmod_zeros(NODES, 1, CHOLMOD_REAL, &c);
    current_cv = cholmod_zeros(NODES, 1, CHOLMOD_REAL, &c);

    double * prev_cv_x    = (double*) prev_cv->x;
    double * current_cv_x = (double*) current_cv->x;
    double * prev_vd_x    = (double*)prev_vd->x; 
    double * current_vd_x = (double*)current_vd->x; 
 
    for (int i = 0; i<SOURCES; i++)
        prev_cv_x[csi[i]] = current_vectors[0][i];

    max_voltage_drops.resize(NODES);
    
    prev_vd = cholmod_solve(CHOLMOD_A, G_L, prev_cv, &c);
    prev_vd_x = (double*)prev_vd->x;
//     cout<<0<<" "<<prev_vd_x[20000]<<endl;
    for (int i = 0; i<NODES; i++)
        max_voltage_drops[i] = prev_vd_x[i];

    cholmod_dense * u, * B_TR_x, * lfs;
    u   = cholmod_ones(NODES, 1, CHOLMOD_REAL, &c);
    lfs = cholmod_ones(NODES, 1, CHOLMOD_REAL, &c);
    double * u_x, * B_TR_x_x, * lfs_x;
    int idx = 1;
    NUM_TIME_POINTS = ceil(sim_time/TIME_STEP);
    cout<<"Number of time points = "<<NUM_TIME_POINTS<<endl;
    double alpha1[1] = {1};
    double alpha2[1] = {0};
   
    B_TR_x = cholmod_zeros(NODES, 1, CHOLMOD_REAL, &c);
    for (double time = TIME_STEP; time <= sim_time; time += TIME_STEP)
    {
        while (true)
        {
            if (time_points[time_index] <= time && time < time_points[time_index+1])
                break;
            time_index++;
        }

        for (int i = 0; i<SOURCES; i++)
        {
            current_cv_x[csi[i]] = current_vectors[time_index][i] 
            +( ( current_vectors[time_index+1][i] - current_vectors[time_index][i])
                    /(time_points[time_index+1]-time_points[time_index])
                * (time - time_points[time_index]) );
        }
        u_x = (double*) u->x;
        prev_cv_x    = (double*)prev_cv->x;
        current_cv_x = (double*)current_cv->x;
        for (int i = 0; i<NODES; i++)
            u_x[i] = 0.5*prev_cv_x[i] + 0.5*current_cv_x[i];
        
//         B_TR_x  = smatrix_dmatrix_multiply(1, B_TR, prev_vd);
        B_TR_x_x = (double*) B_TR_x->x;
        for (int i = 0; i<NODES; i++)
            B_TR_x_x[i] = 0;
        cholmod_sdmult(B_TR, 0, alpha1, alpha2, prev_vd, B_TR_x, &c);


        B_TR_x_x = (double*) B_TR_x->x;
        lfs_x = (double*)lfs->x;
        for (int i = 0; i<NODES; i++)
            lfs_x[i] = u_x[i] + B_TR_x_x[i];
//         lfs     = dmatrix_dmatrix_add(1, u, 1, B_TR_x);
        
        current_vd = cholmod_solve(CHOLMOD_A, A_TR_L, lfs, &c);
        current_vd_x = (double*)current_vd->x;
//         cout<<time<<" "<<current_vd_x[20000]<<endl;
        for (int i = 0; i<NODES; i++)
            if (max_voltage_drops[i] < current_vd_x[i])
                max_voltage_drops[i] = current_vd_x[i];

        // Re-assign
        cholmod_free_dense(&prev_cv, &c);
        cholmod_free_dense(&prev_vd, &c);
 
        prev_cv = cholmod_copy_dense(current_cv, &c);
        prev_vd = cholmod_copy_dense(current_vd, &c);
    
//         cholmod_free_dense(&current_cv, &c);
        cholmod_free_dense(&current_vd, &c);
 
        idx++;
        if (idx % (int)(NUM_TIME_POINTS/100) == 0)
        {
            printf("\r%2d%% complete...", (int)(100.0*idx/NUM_TIME_POINTS));
            fflush(stdout);
            printf("\r");
        }
    }
    NUM_TIME_POINTS = idx;
    rc_tr_stamp("Performing simulation");

    // Memory management
    cholmod_free_dense(&prev_cv, &c);
    cholmod_free_dense(&prev_vd, &c);
    cholmod_free_dense(&current_cv, &c);
    cholmod_free_dense(&current_vd, &c);
    cholmod_free_dense(&u, &c);
    cholmod_free_dense(&B_TR_x, &c);
    cholmod_free_dense(&lfs, &c);
    cholmod_free_sparse(&B_TR, &c);
    cholmod_free_sparse(&A_TR, &c);
    cholmod_free_factor(&G_L, &c);
    cholmod_free_factor(&A_TR_L, &c);

//     printf("\nPrinting the profiler details [timer records]
//             \n---------------------------------------------\n");
    printf("%60s %11s\n", "Task Description", "Wall (s)");
    printf("%60s %11s\n", "----------------", "--------");
    for (unsigned int i = 1; i<rc_tr_stamps.size(); i++)
    {
        printf("%60s %11.2f\n", rc_tr_stamps[i].second.c_str(), rc_tr_stamps[i].first-rc_tr_stamps[i-1].first);
    }

    cholmod_finish(&c);
}

vector<double> rc_tr::get_max_voltage_drops()
{
    if (!SAVE_MEMORY)
    {
        vector<double> rv(NODES, 0);
        for (unsigned int i = 0; i<voltage_drop_vectors.size(); i++)
        {
            double * vd = (double*)voltage_drop_vectors[i]->x;
            for (int q = 0; q<NODES; q++)
                if (rv[q] < vd[q])
                    rv[q] = vd[q];
        }
        return rv;
    }
    else
        return max_voltage_drops;
}


void rc_tr::compute_voltage_drops()
{
    // Find B_TR 
    B_TR = smatrix_smatrix_add(1.0/TIME_STEP, C, -0.5, G); 

    // Find A_TR
    A_TR = smatrix_smatrix_add(0.5, G, 1.0/TIME_STEP, C);

    // Factorize A_TR
    A_TR_L = cholmod_analyze(A_TR, &c);     // symbolic factorization
    cholmod_factorize(A_TR, A_TR_L, &c);    // numeric factorization
    
    // Factorize G
    G_L = cholmod_analyze(G, &c);           // symbolic factorization
    cholmod_factorize(G, G_L, &c);          // numeric factorization

    // Allocate memory for the result
    voltage_drop_vectors.resize(NUM_TIME_POINTS);

    // DC operating point 
    voltage_drop_vectors[0] = cholmod_solve(CHOLMOD_A, G_L, fts_currents[0], &c);
    double * xx = (double*) voltage_drop_vectors[0]->x;
    cout<<xx[0]<<" ";
    // Recursion 
    cholmod_dense * u, * B_TR_x, * lfs;
    for (int i = 1; i<NUM_TIME_POINTS; i++)
    {
        u       = dmatrix_dmatrix_add(0.5, fts_currents[i-1], 0.5, fts_currents[i]);
        B_TR_x  = smatrix_dmatrix_multiply(1, B_TR, voltage_drop_vectors[i-1]);
        lfs     = dmatrix_dmatrix_add(1, u, 1, B_TR_x);
        voltage_drop_vectors[i] = cholmod_solve(CHOLMOD_A, A_TR_L, lfs, &c);
        xx = (double*)voltage_drop_vectors[i]->x;
        cout<<xx[0]<<" ";
    }

    cout<<"Number of time points = "<<NUM_TIME_POINTS<<endl;

    // Memory management
    cholmod_free_sparse(&B_TR, &c);
    cholmod_free_sparse(&A_TR, &c);
    cholmod_free_factor(&G_L, &c);
    cholmod_free_factor(&A_TR_L, &c);
    for (unsigned int i = 0; i<fts_currents.size(); i++)
        cholmod_free_dense(&fts_currents[i], &c);
}

void rc_tr::get_current_points(vector<double> & time_points, vector<vector<double> > & current_vectors)
{
    sim_time = time_points[time_points.size()-1];
    int time_index = 0;
    for (double time = 0; time <= sim_time; time += TIME_STEP)
    {
        while (true)
        {
            if (time_points[time_index] <= time && time < time_points[time_index+1])
                break;
            time_index++;
        }
        cholmod_dense * new_current_vectors = cholmod_zeros(NODES, 1, CHOLMOD_REAL, &c);
        double * cvx = (double*)new_current_vectors->x;
        for (int i = 0; i<SOURCES; i++)
        {
            cvx[csi[i]] = current_vectors[time_index][i] 
            +( ( current_vectors[time_index+1][i] - current_vectors[time_index][i])
                    /(time_points[time_index+1]-time_points[time_index])
                * (time - time_points[time_index]) );
            cout<<cvx[csi[i]]<<endl;
        }
        fts_currents.push_back(new_current_vectors);
    }
    NUM_TIME_POINTS = fts_currents.size();
}



cholmod_sparse * rc_tr::smatrix_smatrix_add(double a1, cholmod_sparse *A1,
                                                    double a2, cholmod_sparse *A2)
{
    double alpha[1] = {a1}; 
    double beta [1] = {a2}; 
    cholmod_sparse * A3 = cholmod_add(A1, A2, alpha, beta, true, false, &c);
    return A3;
}

cholmod_dense  * rc_tr::smatrix_dmatrix_multiply(double a1, cholmod_sparse *A1,
                                                         cholmod_dense  *A2)       
{
    cholmod_dense *A3 = cholmod_zeros(A1->nrow, A2->ncol, CHOLMOD_REAL, &c);
    double alpha[1] = {a1};
    double beta [1] = {0};
    cholmod_sdmult(A1, 0, alpha, beta, A2, A3, &c);
    return A3;
}

cholmod_dense * rc_tr::dmatrix_dmatrix_add(double a1, cholmod_dense *A1,
                                            double a2, cholmod_dense *A2)
{
    // Sanity check
    assert((A1->nrow == A2->nrow) && (A1->ncol == A2->ncol));

    // Add
    cholmod_dense *A3 = cholmod_zeros(A1->nrow, A1->ncol, CHOLMOD_REAL, &c);
    double * A1x = (double*)A1->x;
    double * A2x = (double*)A2->x;
    double * A3x = (double*)A3->x;
    for (unsigned int i = 0; i<(A1->nrow * A1->ncol); i++)
        A3x[i] = a1*A1x[i] + a2*A2x[i];

    return A3;
}
void rc_tr_stamp(string op)
{
    pair<double, string> st;
    st.first  = omp_get_wtime();
    st.second = op;
    rc_tr_stamps.push_back(st);
}
