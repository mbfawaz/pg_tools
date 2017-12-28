#include <sstream>
#include <fstream>
#include "rlc_tr.h"
#include "cholmod.h"
#include "string.h"
#include "omp.h"
#include "math.h"
#include "mosek_interface.h"
#include "rlc_grid.h"
#include "waveforms.h"
#include <cmath>
#include "assert.h"
#include "umfpack.h"

using namespace std;

bool RLC_SAVE_MEMORY = 1;

vector<pair<double, string> > rlc_tr_stamps;
void rlc_tr_stamp(std::string op);

rlc_tr::rlc_tr()
{
    cholmod_start(&c);
}

rlc_tr::rlc_tr(rlc_grid &g, waveforms &W)
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
    TIME_STEP           = W.MIN_TIME_STEP;

    if (!RLC_SAVE_MEMORY)
        get_current_points(W.time_points, W.current_vectors);

    // Build new G
    cholmod_sparse * M_T = cholmod_allocate_sparse(INDUCTORS, NODES, M->nzmax, 0, 1, 0, CHOLMOD_REAL, &c);
    cholmod_transpose_unsym(M, 1, NULL, NULL, 0, M_T, &c);
    
    cholmod_sparse * minus_M = cholmod_copy_sparse(M, &c);
    double * minus_M_x = (double*) minus_M->x;
    for (unsigned int i = 0; i<minus_M->nzmax; i++)
        minus_M_x[i] *= -1;
    
    cholmod_sparse * zeros = cholmod_allocate_sparse(INDUCTORS, INDUCTORS, 0, 0, 1, 0, CHOLMOD_REAL, &c);

    cholmod_sparse * row1 = cholmod_horzcat(G, minus_M, true, &c);
    cholmod_sparse * row2 = cholmod_horzcat(M_T, zeros, true, &c);
    
    GG = cholmod_vertcat(row1, row2, true, &c); 
    cholmod_free_sparse(&M_T, &c);
    cholmod_free_sparse(&minus_M, &c);
    cholmod_free_sparse(&row1, &c);
    cholmod_free_sparse(&row2, &c);
    cholmod_free_sparse(&zeros, &c);


    // Build new C
    CC = cholmod_speye(NODES+INDUCTORS, NODES+INDUCTORS, CHOLMOD_REAL, &c);
    double * CCx = (double*)CC->x; 
    double * Cx  = (double*)C->x;
    double * Lx  = (double*)L->x;
    for (int i = 0; i<NODES; i++)
        CCx[i] = Cx[i];
    for (int i = 0; i<INDUCTORS; i++)
        CCx[NODES + i] = Lx[i];
}

void rlc_tr::compute_voltage_drops_sm(vector<double> & time_points, vector<vector<double> > & current_vectors)
{   
    rlc_tr_stamp("Start");

    // Find B_TR 
    cout<<"Finding B_TR..."<<endl;
    B_TR = smatrix_smatrix_add(1.0/TIME_STEP, CC, -0.5, GG); 
    rlc_tr_stamp("Finding B_TR");

    // Find A_TR
    cout<<"Finding A_TR..."<<endl;
    A_TR = smatrix_smatrix_add(0.5, GG, 1.0/TIME_STEP, CC);
    rlc_tr_stamp("Finding A_TR");

    // Factorize GG
    cout<<"Factorizing GG..."<<endl;
    void * GG_sym;
    (void) umfpack_di_symbolic(NODES+INDUCTORS, NODES+INDUCTORS, (int*)GG->p, (int*)GG->i, (double*)GG->x, &GG_sym, NULL, NULL);
    (void) umfpack_di_numeric ((int*)GG->p, (int*)GG->i, (double*)GG->x, GG_sym, &GG_num, NULL, NULL);
    umfpack_di_free_symbolic(&GG_sym);
    rlc_tr_stamp("Factorizing GG");

    double * prev_cv;
    double * current_cv;
    double * prev_vd; 
    double * current_vd;
    prev_cv    = (double*)calloc(NODES+INDUCTORS, sizeof(double));
    prev_vd    = (double*)calloc(NODES+INDUCTORS, sizeof(double));
    current_cv = (double*)calloc(NODES+INDUCTORS, sizeof(double));
    current_vd = (double*)calloc(NODES+INDUCTORS, sizeof(double));

    for (int i = 0; i<SOURCES; i++)
        prev_cv[csi[i]] = current_vectors[0][i];

    (void) umfpack_di_solve(UMFPACK_A, (int*)GG->p, (int*)GG->i, (double*)GG->x, prev_vd, prev_cv, GG_num, NULL, NULL);
    umfpack_di_free_numeric(&GG_num);

    // Factorize A_TR
    cout<<"Factorizing A_TR..."<<endl;
    void * A_TR_sym;
    (void) umfpack_di_symbolic(NODES+INDUCTORS, NODES+INDUCTORS, (int*)A_TR->p, (int*)A_TR->i, (double*)A_TR->x, &A_TR_sym, NULL, NULL);
    (void) umfpack_di_numeric ((int*)A_TR->p, (int*)A_TR->i, (double*)A_TR->x, A_TR_sym, &A_TR_num, NULL, NULL);
    umfpack_di_free_symbolic(&A_TR_sym);
    rlc_tr_stamp("Factorizing A_TR");

    // DC operating point
    cout<<"Performing simulation..."<<endl;
    sim_time = time_points[time_points.size()-1];
    int time_index = 0;
 

    double * u      = (double*)calloc(NODES+INDUCTORS, sizeof(double));
    double * B_TR_x = (double*)calloc(NODES+INDUCTORS, sizeof(double));
    double * lfs    = (double*)calloc(NODES+INDUCTORS, sizeof(double));

    FILE * fout;
    fout = fopen("wave.m", "w");

    int CHECK_NODE = NODES-1;
    fprintf(fout, "%% Node: %u\n", CHECK_NODE);
    fprintf(fout, "ve = [");
//         fprintf(fout, "%e %f\n", i*TIME_STEP, fts_currents[i][0]);

//         cout<<voltage_drop_vectors[i][5389]<<endl;

//     int node = NODES - INDUCTORS/2;

    max_voltage_drops.resize(NODES);
    min_voltage_drops.resize(NODES);
    int idx = 0; 
    NUM_TIME_POINTS = ceil(sim_time/TIME_STEP);
    cout<<"Number of time points = "<<NUM_TIME_POINTS<<endl;
 
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
            current_cv[csi[i]] = current_vectors[time_index][i] 
            +( ( current_vectors[time_index+1][i] - current_vectors[time_index][i])
                    /(time_points[time_index+1]-time_points[time_index])
                * (time - time_points[time_index]) );
        }
        for (int i = 0; i<NODES+INDUCTORS; i++)
            u[i] = 0.5*prev_cv[i] + 0.5*current_cv[i];

        for (int j = 0; j<NODES+INDUCTORS; j++)
            B_TR_x[j] = 0;
        for (int j = 0; j<NODES+INDUCTORS; j++)
        {
            for (int p = ((int*)B_TR->p)[j]; p<((int*)B_TR->p)[j+1]; p++)
            {
                B_TR_x[((int*)B_TR->i)[p]] += ((double*)B_TR->x)[p] * prev_vd[j];
            }
        }
        for (int j = 0; j<NODES+INDUCTORS; j++)
            lfs[j] = u[j] + B_TR_x[j];

        (void) umfpack_di_solve(UMFPACK_A, (int*)A_TR->p, (int*)A_TR->i, (double*)A_TR->x, current_vd, lfs, A_TR_num, NULL, NULL);
        
        fprintf(fout, "%e %f\n", time, current_vd[CHECK_NODE]);
        for (int j = 0; j<NODES; j++)
        {
            if (max_voltage_drops[j] < current_vd[j])
                max_voltage_drops[j] = current_vd[j];
            if (min_voltage_drops[j] > current_vd[j])
                min_voltage_drops[j] = current_vd[j];
        }
        //prev_cv = current_cv;
        //prev_vd = current_vd;
        for (int i = 0; i<NODES+INDUCTORS; i++)
        {
            prev_cv[i] = current_cv[i];
            prev_vd[i] = current_vd[i];
        }

        idx++;
        if (idx % (int)(NUM_TIME_POINTS/100) == 0)
        {
            printf("\r%2d%% complete...", (int)(100.0*idx/NUM_TIME_POINTS));
            fflush(stdout);
            printf("\r");
        }
    } 
    fprintf(fout, "];\n\n");
    fclose(fout);

    free(prev_cv);
    free(prev_vd);
    free(current_cv);
    free(current_vd);

    rlc_tr_stamp("Peforming simulation");

    // Memory management
    free(u); free(B_TR_x); free(lfs);
    cholmod_free_sparse(&B_TR, &c);
    cholmod_free_sparse(&A_TR, &c);
    umfpack_di_free_numeric(&GG_num);
    umfpack_di_free_numeric(&A_TR_num);

    for (unsigned int i = 0; i<fts_currents.size(); i++)
        free(fts_currents[i]);

    cout<<endl;
    printf("\nPrinting the profiler details [timer records]\
            \n---------------------------------------------\n");
    printf("%60s %11s\n", "Task Description", "Wall (s)");
    printf("%60s %11s\n", "----------------", "--------");
    for (unsigned int i = 1; i<rlc_tr_stamps.size(); i++)
    {
        printf("%60s %11.2f\n", rlc_tr_stamps[i].second.c_str(), rlc_tr_stamps[i].first-rlc_tr_stamps[i-1].first);
    }
}

vector<double> rlc_tr::get_max_voltage_drops()
{
    if (!RLC_SAVE_MEMORY)
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

vector<double> rlc_tr::get_min_voltage_drops()
{
    if (!RLC_SAVE_MEMORY)
    {
        vector<double> rv(NODES, 0);
        for (unsigned int i = 0; i<voltage_drop_vectors.size(); i++)
        {
            double * vd = (double*)voltage_drop_vectors[i]->x;
            for (int q = 0; q<NODES; q++)
                if (rv[q] > vd[q])
                    rv[q] = vd[q];
        }
        return rv;
    }
    else
        return min_voltage_drops;
}



void rlc_tr::compute_voltage_drops()
{
    ;
}

void rlc_tr::get_current_points(vector<double> & time_points, vector<vector<double> > & current_vectors)
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



cholmod_sparse * rlc_tr::smatrix_smatrix_add(double a1, cholmod_sparse *A1,
                                                    double a2, cholmod_sparse *A2)
{
    double alpha[1] = {a1}; 
    double beta [1] = {a2}; 
    cholmod_sparse * A3 = cholmod_add(A1, A2, alpha, beta, true, false, &c);
    return A3;
}

cholmod_dense  * rlc_tr::smatrix_dmatrix_multiply(double a1, cholmod_sparse *A1,
                                                         cholmod_dense  *A2)       
{
    cholmod_dense *A3 = cholmod_zeros(A1->nrow, A2->ncol, CHOLMOD_REAL, &c);
    double alpha[1] = {a1};
    double beta [1] = {0};
    cholmod_sdmult(A1, 0, alpha, beta, A2, A3, &c);
    return A3;
}

cholmod_dense * rlc_tr::dmatrix_dmatrix_add(double a1, cholmod_dense *A1,
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
void rlc_tr_stamp(string op)
{
    pair<double, string> st;
    st.first  = omp_get_wtime();
    st.second = op;
    rlc_tr_stamps.push_back(st);
}
