#include <iostream>

// Core
#include "sparse_matrix.h"
#include "container.h"
#include "mosek_interface.h"
#include "rc_grid.h"
#include "rlc_grid.h"

// Projects
#include "rc_vectorless.h"
#include "rlc_vectorless.h"
#include "rc_simulation.h"
#include "rlc_simulation.h"
#include "rc_tran_constraints.h"

using namespace std;

int project = 7;

int main(int argc, char** argv)
{
    // Sanity check
    if (argc < 2)
    {
        cout<<"Please specify a grid options file"<<endl;
        exit(0);
    }
    if (argc < 3)
    {
        cout<<"Please specify an options file for the problem"<<endl;
        exit(0);
    }
    
    switch(project)
    {
        case 3:
        {
            cout<<"Vectorless RC Verification"<<endl;
            rc_grid g(argv[1]);
            rc_vectorless * rc = new rc_vectorless(g, argv[2]);
            delete rc;
            break;
        }
        case 4:
        {
            cout<<"Vectorless RLC Verification"<<endl;
            rlc_grid g(argv[1]);
            rlc_vectorless * rlc = new rlc_vectorless(g, argv[2]);
            delete rlc;
            break;
        }
        case 5:
        {
            cout<<"Vector-based RC Verification"<<endl;
            rc_grid g(argv[1]);
            rc_simulation * rc = new rc_simulation(g, argv[2]);
            delete rc;
            break;
        }
        case 6:
        {
            cout<<"Vector-based RLC Verification"<<endl;
            rlc_grid g(argv[1]);
            rlc_simulation * rlc = new rlc_simulation(g, argv[2]);
            delete rlc;
            break;
        }
        case 7:
        {
            cout<<"Vectorless RC Verification - Transient Constraints"<<endl;
            rc_grid g(argv[1]);
            rc_tran_constraints * rc = new rc_tran_constraints(g, argv[2]);
            delete rc;
            break;
        }
        default:
        {
            cout<<"Invalid project selected."<<endl;
            break;
        }
    }
    return 0;
}
