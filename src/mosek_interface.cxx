#include "mosek_interface.h"
#include <vector>

using namespace std;

/* This function prints log output from MOSEK to the terminal. */
static void MSKAPI printstr(void *handle, MSKCONST char str[])
{
    printf("%s",str);
} /* printstr */


//Default constructor
lp_opt_problem::lp_opt_problem()
{
    env  = NULL;
    task = NULL;
}

lp_opt_problem::~lp_opt_problem()
{
    ;
}


//Set number of variables
void lp_opt_problem::set_numvar(int _numvar)
{
    numvar = _numvar;
}

//Set number of variables
void lp_opt_problem::set_numcon(int _numcon)
{
    numcon = _numcon;
}


// Create and environment
void lp_opt_problem::create_env()
{
    r = MSK_makeenv(&env, NULL);
}

//Create and optimization task
void lp_opt_problem::create_opt_task()
{
    r = MSK_maketask(env, numcon, numvar, &task);

}

//Links the mosek stream to the printing function
void lp_opt_problem::request_mosek_log()
{
    if ( r==MSK_RES_OK ) 
        r = MSK_linkfunctotaskstream(task,MSK_STREAM_LOG,NULL,printstr);
}

//Set optimization method (simplex or interior-point)
void lp_opt_problem::set_opt_method(int _opt_method)
{
    if (r == MSK_RES_OK)
    {
        opt_method = _opt_method;

        r = MSK_putintparam(task,MSK_IPAR_PRESOLVE_USE,MSK_PRESOLVE_MODE_OFF); 
        if (_opt_method == SIMPLEX)
            MSK_putintparam(task, MSK_IPAR_OPTIMIZER, MSK_OPTIMIZER_PRIMAL_SIMPLEX);
        else if (_opt_method == INTERIOR)
        {
            MSK_putintparam(task, MSK_IPAR_INTPNT_BASIS, MSK_BI_NEVER);            
            MSK_putintparam(task, MSK_IPAR_OPTIMIZER, MSK_OPTIMIZER_INTPNT);
        }
	    else
        {
            printf("Invalid optimization method.\n");
            exit(0); 
        }
    }
}

//Set optimization method (simplex or interior-point)
void lp_opt_problem::set_opt_method_vanilla(int _opt_method)
{
    if (r == MSK_RES_OK)
        opt_method = _opt_method;

    if (_opt_method == INTEGER)
        MSK_putintparam(task, MSK_IPAR_OPTIMIZER, MSK_OPTIMIZER_MIXED_INT);
//     if (_opt_method == CONINT)
//         MSK_putintparam(task, MSK_IPAR_OPTIMIZER, MSK_OPTIMIZER_MIXED_INT_CONIC);
        
}

//Read input data
void lp_opt_problem::set_A  (vector<int> &Ap,   vector<int> &Ai,   vector<double> &Ax)
{
    if (numcon > 0) 
    {
        if (Ap.size() != 0 && Ai.size() != 0 && Ax.size() != 0)
        {
            aptrb = &Ap[0];
            asub  = &Ai[0];
            aval  = &Ax[0];
        }
        else
        {   
            printf("Invalid A matrix.\n");
            exit(0);
        }   
    }
}   

//Read input data
void lp_opt_problem::set_A_2  (int *Ap, int * Ai, double * Ax)
{
    if (numcon > 0) 
    {
        aptrb = Ap;
        asub  = Ai;
        aval  = Ax;
    }
}

void lp_opt_problem::set_bounds_on_constraints(vector<double> &b_0, vector<double> &b_1)
{
    if (numcon > 0)
    {
        if ((MSKint32t)b_0.size() == numcon && (MSKint32t)b_1.size() == numcon)
        {   
            blc = &b_0[0];
            buc = &b_1[0];
            bkc = new MSKboundkeye[numcon];
            for (int i = 0; i<numcon; i++)
            {
                if(blc[i] == -MSK_INFINITY &&  buc[i] == MSK_INFINITY)
			        bkc[i] = MSK_BK_FR;
		        else if(blc[i] == -MSK_INFINITY)
			        bkc[i] = MSK_BK_UP;
		        else if(buc[i] ==  MSK_INFINITY)
			        bkc[i] = MSK_BK_LO;
		        else if(blc[i] == buc[i])
		    	    bkc[i] = MSK_BK_FX;
		        else
			        bkc[i] = MSK_BK_RA;
            }
        }
        else
        {
            printf("Inconsistent bounds on constraints.\n");
            exit(0);
        }
    }
}

void lp_opt_problem::set_bounds_on_variables(vector<double> &x_0, vector<double> &x_1)
{
    if (numvar > 0)
    {
        if ((MSKint32t)x_0.size() == numvar && (MSKint32t)x_1.size() == numvar)
        {   
            blx = &x_0[0];
            bux = &x_1[0];
            bkx = new MSKboundkeye[numvar];
            for (int i = 0; i<numvar; i++)
            {
                if(blx[i] == -MSK_INFINITY &&  bux[i] == MSK_INFINITY)
			        bkx[i] = MSK_BK_FR;
		        else if(blx[i] == -MSK_INFINITY)
			        bkx[i] = MSK_BK_UP;
		        else if(bux[i] ==  MSK_INFINITY)
			        bkx[i] = MSK_BK_LO;
		        else if(blx[i] == bux[i])
		    	    bkx[i] = MSK_BK_FX;
		        else
			        bkx[i] = MSK_BK_RA;
            }
        }
        else
        {
            printf("Inconsistent bounds on variables.\n");
            exit(0);
        }
    }

}

//Set objective function
void lp_opt_problem::set_obj(vector<double> &obj)
{
    if (obj.size() > 0)
    {
        if ((MSKint32t)obj.size() == numvar)
            c = &obj[0]; 
        else
        {
            printf("Invalid objective function size.\n");
            exit(0);
        }
    }
    else
    {
        printf("Invalid objective function.\n");
        exit(0);
    }
}



//Load a problem into the task object
void lp_opt_problem::load_constraints_into_task()
{
    //Append 'numcon' empty constraints
    //The constraints will initially have no bounds 
    if (r == MSK_RES_OK)
        r = MSK_appendcons(task, numcon);

    //Append 'numvar' empty constraints 
    //The variables will initially be fixed at zero (x=0)
    if (r == MSK_RES_OK)
        r = MSK_appendvars(task, numvar);

    //Loop over all the variables
    for (int j = 0; j<numvar && r == MSK_RES_OK; ++j)
    {
        //Set the bounds on variable j
        //blx[j] <= x_j <= bux[j]
        if (r == MSK_RES_OK)
            r = MSK_putvarbound(task, j, bkx[j], blx[j], bux[j]);

        //Input column j of A
        if (r == MSK_RES_OK)
            r = MSK_putacol    (task, j,
                                aptrb[j+1]-aptrb[j], //Number of non-zeros in column j
                                asub+aptrb[j],       //Pointer to row indexes of column j
                                aval+aptrb[j]);      //Poiter to Values of column j
    }
    delete [] bkx;

    //Set the bounds on constraints
    //for i=1,...,numcon: blc[j] <= constraint i <= buc[i]
    for (int i = 0; i<numcon && r == MSK_RES_OK; ++i)
        r = MSK_putconbound(task, i, bkc[i], blc[i], buc[i]);
    delete [] bkc;
}

void lp_opt_problem::load_objective_into_task()
{
    for (int j = 0; j<numvar && r == MSK_RES_OK; ++j)
    {
        //Set the linear term c_j in the objective
        if (r == MSK_RES_OK)
            r = MSK_putcj(task, j, c[j]);
    }
}

//Set optimization sense
void lp_opt_problem::set_opt_sense(int _opt_sense)
{
    if (r == MSK_RES_OK)
    {
        opt_sense = _opt_sense;
        if      (opt_sense == MINIMIZE)
            r = MSK_putobjsense(task, MSK_OBJECTIVE_SENSE_MINIMIZE);
        else if (opt_sense == MAXIMIZE)
            r = MSK_putobjsense(task, MSK_OBJECTIVE_SENSE_MAXIMIZE);
        else
        {   
            printf("Invalid optimization sense.\n");
            exit(0);
        }
    }
}


//Optimization
void lp_opt_problem::optimize()
{
    //Run optimizer
    if (r == MSK_RES_OK)
    {
        MSKrescodee trmcode;

        //Run optimzer
        r = MSK_optimizetrm(task, &trmcode);
        
        // Solution Summary
        MSK_solutionsummary (task,MSK_STREAM_LOG);
        
        //Get solution status
        if (r == MSK_RES_OK)
        {
            MSKsolstae solsta;

            if (r == MSK_RES_OK)
            {
                if (opt_method == SIMPLEX) 
                    r = MSK_getsolsta(task, MSK_SOL_BAS, &solsta);
                else if (opt_method == INTERIOR)
                    r = MSK_getsolsta(task, MSK_SOL_ITR, &solsta);
                
            }
            switch (solsta)
            {
                case MSK_SOL_STA_OPTIMAL:   
                case MSK_SOL_STA_NEAR_OPTIMAL:
                    //printf("Optimal solution found.\n");
                    problem_solved = true;
                    break;
                case MSK_SOL_STA_DUAL_INFEAS_CER:
                case MSK_SOL_STA_PRIM_INFEAS_CER:
                case MSK_SOL_STA_NEAR_DUAL_INFEAS_CER:
                case MSK_SOL_STA_NEAR_PRIM_INFEAS_CER:  
                    printf("Primal or dual infeasibility certificate found.\n");
                    break;
                case MSK_SOL_STA_UNKNOWN:
                    char symname[MSK_MAX_STR_LEN];
                    char desc[MSK_MAX_STR_LEN];
                    MSK_getcodedesc(trmcode, symname, desc);
                    printf("The solution status is unknown.\n");
                    printf("The optimizer terminitated with code: %s\n",symname);
                    break;
                default:
                    printf("Other solution status.\n");
                    break;
            }
        }
        
        if (r != MSK_RES_OK)
        {
            /* In case of an error print error code and description. */      
            char symname[MSK_MAX_STR_LEN];
            char desc[MSK_MAX_STR_LEN];
      
            printf("An error occurred while optimizing.\n");     
            MSK_getcodedesc (r, symname, desc);
            printf("Error %s - '%s'\n",symname,desc);
        }   
    }
}


 
//Extracting the solution
//Get optimized variables
vector<double> lp_opt_problem::get_optimized_variables()
{
    if (r == MSK_RES_OK && problem_solved)
    {
        optimized_variables.resize(numvar);
        double *xx = (double*) calloc(numvar, sizeof(double));
        if (opt_method == SIMPLEX)
            MSK_getxx(task, MSK_SOL_BAS, xx); 
        else if (opt_method == INTERIOR)
            MSK_getxx(task, MSK_SOL_ITR, xx);
        else if (opt_method == INTEGER || opt_method == CONINT)
            MSK_getxx(task, MSK_SOL_ITG, xx);
        for (int i = 0; i<numvar; i++)
            optimized_variables[i] = xx[i];
        return optimized_variables;
    }
    else
    {
        vector<double> temp(numvar, MSK_INFINITY);
        printf("No solution available.\n");
	    return temp ;
    }
}

//Get optimized objective
double lp_opt_problem::get_optimized_objective()
{
    if (r == MSK_RES_OK && problem_solved)
    {
        if (opt_method == SIMPLEX)
            MSK_getprimalobj(task, MSK_SOL_BAS, &optimized_objective);
        else if (opt_method == INTERIOR)
            MSK_getprimalobj(task, MSK_SOL_ITR, &optimized_objective);
        else if (opt_method == INTEGER || opt_method == CONINT)
            MSK_getprimalobj(task, MSK_SOL_ITG, &optimized_objective);
        return optimized_objective;
    }
    else
    {
        printf("No solution available.\n");
        return MSK_INFINITY;
    }
}


void lp_opt_problem::delete_task_and_env()
{
    MSK_deletetask(&task);
    MSK_deleteenv(&env);
}

void cq_opt_problem::append_cone(bool type, vector<int> cs)
{
    if (r == MSK_RES_OK)
    {
        int num_var_in_cone = cs.size();
        MSKint32t * csub = new MSKint32t[num_var_in_cone];

        for (int i = 0; i<num_var_in_cone; i++)
            csub[i] = cs[i];
        
        if (type == QUADRATIC)
            r = MSK_appendcone(task, MSK_CT_QUAD,  0.0, num_var_in_cone, csub);
        else
            r = MSK_appendcone(task, MSK_CT_RQUAD, 0.0, num_var_in_cone, csub);
    }            
}

void ilp_opt_problem::set_int_var(vector<int> intsub)
{
    for (unsigned int j = 0; j<intsub.size() && r == MSK_RES_OK; j++)
        r = MSK_putvartype(task, intsub[j], MSK_VAR_TYPE_INT);
}


void ilp_opt_problem::append_cone(bool type, vector<int> cs)
{
    if (r == MSK_RES_OK)
    {
        int num_var_in_cone = cs.size();
        MSKint32t * csub = new MSKint32t[num_var_in_cone];

        for (int i = 0; i<num_var_in_cone; i++)
            csub[i] = cs[i];
        
        if (type == QUADRATIC)
            r = MSK_appendcone(task, MSK_CT_QUAD,  0.0, num_var_in_cone, csub);
        else
            r = MSK_appendcone(task, MSK_CT_RQUAD, 0.0, num_var_in_cone, csub);
    }            
}

void ilp_opt_problem::set_initial_feasible_solution(vector<double> init_sol, int first, int last)
{
    if (r == MSK_RES_OK) 
        r = MSK_putintparam(task,MSK_IPAR_MIO_CONSTRUCT_SOL,MSK_ON); 
 
    if (r == MSK_RES_OK) 
    { 
        double *xx = new double[init_sol.size()];
        for (unsigned int i = 0; i<init_sol.size(); i++)
            xx[i] = init_sol[i];
    
        /* Assign values 0,2,0 to integer variables */ 
        r = MSK_putxxslice(task,MSK_SOL_ITG,first,last,xx); 
    } 
}




