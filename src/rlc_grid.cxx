#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <string.h>
#include "rlc_grid.h"

using namespace std;

rlc_grid::rlc_grid()
{
    cholmod_start(&c);
}

rlc_grid::rlc_grid(string options_file)
{
    cholmod_start(&c);

    load_grid(options_file);
}

// Load a grid given the name of its config file
// The config file will tell all the parameters of the grid
// These will be used to determine the destination of the grid matrices/vectors
// The function will then load all those into cholmod_sparse structures and vectors
void rlc_grid::load_grid(string options_file)
{
    ifstream infile;
    infile.open(options_file.c_str());
    if (!infile)
    {
        cout<<"ERROR: could not open grid options file."<<endl;
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

                 if (strcmp(param_name, "number_of_layers")                == 0) number_of_layers               = param_value_double;
            else if (strcmp(param_name, "number_of_layers_global_grid")    == 0) number_of_layers_global_grid   = param_value_double;
            else if (strcmp(param_name, "number_of_sub_grids_X_DIM")       == 0) number_of_sub_grids_X_DIM      = param_value_double;
            else if (strcmp(param_name, "number_of_sub_grids_Y_DIM")       == 0) number_of_sub_grids_Y_DIM      = param_value_double;
            else if (strcmp(param_name, "GRID_DIMENSION_X_UM")             == 0) GRID_DIMENSION_X_UM            = param_value_double;
            else if (strcmp(param_name, "GRID_DIMENSION_Y_UM")             == 0) GRID_DIMENSION_Y_UM            = param_value_double;
        }
    }

    cout<<"Grid parameters as read from the grid options file"<<endl;
    cout<<"--------------------------------------------------"<<endl;
    cout<<"number_of_layers             = "<<number_of_layers              <<endl;
    cout<<"number_of_layers_global_grid = "<<number_of_layers_global_grid  <<endl;
    cout<<"number_of_sub_grids_X_DIM    = "<<number_of_sub_grids_X_DIM     <<endl;
    cout<<"number_of_sub_grids_Y_DIM    = "<<number_of_sub_grids_Y_DIM     <<endl;
    cout<<"GRID_DIMENSION_X_UM          = "<<GRID_DIMENSION_X_UM           <<endl;
    cout<<"GRID_DIMENSION_Y_UM          = "<<GRID_DIMENSION_Y_UM           <<endl;

    stringstream gridname, ss_g, ss_c, ss_l, ss_m, ss_csi;

    // Destination and grid name
    string destination = "//autofs//fs1.ece//fs1.eecg.najm//r//r0//fawazmoh//Desktop//rlc_power_grids//RLC_";
    gridname<<GRID_DIMENSION_X_UM<<"_"
            <<GRID_DIMENSION_Y_UM<<"_"
            <<number_of_layers<<"_"
            <<number_of_layers_global_grid<<"_"
            <<number_of_sub_grids_X_DIM<<"_"
            <<number_of_sub_grids_Y_DIM;

    // Read grid files
    FILE *f_g, *f_c, *f_l, *f_m, *f_csi;
    ss_g  <<destination<<gridname.str()<<"//g_RLC_"  <<gridname.str();
    ss_c  <<destination<<gridname.str()<<"//c_RLC_"  <<gridname.str();
    ss_l  <<destination<<gridname.str()<<"//l_RLC_"  <<gridname.str();
    ss_m  <<destination<<gridname.str()<<"//m_RLC_"  <<gridname.str();
    ss_csi<<destination<<gridname.str()<<"//csi_RLC_"<<gridname.str();
    
    f_g   = fopen(ss_g.str().  c_str(), "r");
    if (!f_g)
    {
        cout<<"ERROR: could not open the G file."<<endl;
        exit(1);
    }
    G     = cholmod_read_sparse(f_g, &c); G->stype = 1;
    fclose(f_g);
    
    f_c   = fopen(ss_c.str().  c_str(), "r");
    if (!f_c)
    {
        cout<<"ERROR: could not open the C file."<<endl;
        exit(1);
    }
    C     = cholmod_read_sparse(f_c, &c); C->stype = 1;
    fclose(f_c);

    f_l   = fopen(ss_l.str().  c_str(), "r");
    if (!f_l)
    {
        cout<<"ERROR: could not open the L file."<<endl;
        exit(1);
    }
    L     = cholmod_read_sparse(f_l, &c); L->stype = 1;
    fclose(f_l);

    f_m   = fopen(ss_m.str().  c_str(), "r");
    if (!f_m)
    {
        cout<<"ERROR: could not open the M file."<<endl;
        exit(1);
    }
    M     = cholmod_read_sparse(f_m, &c); 
    fclose(f_m);

    f_csi = fopen(ss_csi.str().c_str(), "r");
    if (!f_csi)
    {
        cout<<"ERROR: could not open the CSI file."<<endl;
        exit(1);
    }
    int idx = 0;
    while(fscanf(f_csi, "%d", &idx) != EOF)
        csi.push_back(idx);        
    fclose(f_csi);

    number_of_nodes     = G->nrow;
    number_of_sources   = csi.size();
    number_of_inductors = L->nrow;
    
    cout<<endl<<"Resulting grid dimensions"<<endl;
    cout      <<"-------------------------"<<endl;
    cout<<"Number of nodes              = "<<number_of_nodes<<endl;
    cout<<"Number of sources            = "<<number_of_sources<<endl;
    cout<<"Number of inductors          = "<<number_of_inductors<<endl;
    cout<<endl;
}
