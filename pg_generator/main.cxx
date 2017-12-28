/**
* \file main.cxx \brief Testing code
*
* This file is to run some simulations for the purpose of testing
*/
#include <getopt.h>
#include "power_grid.h"
#include "sparse_matrix.h"
#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include "profiler.h"
#include "common.h"

using namespace std;
using namespace matrix;

//! Function to parse the input arguments
void parse_argument_list (int argc, char** argv);

void print_matrix_matlab(FILE* matlab, matrix::sparse_matrix M, string M_name);
void print_vector_matlab(FILE* matlab, vector<double> v, string v_name);
void print_vector_matlab(FILE* matlab, vector<int> v, string v_name);

void print_triplet_to_file(sparse_matrix &A, string filename);
void print_nonsym_triplet_to_file(sparse_matrix &A, string filename);

bool RLC = false;

int main(int argc, char** argv)
{
    if (!RLC)
    {
        cout<<endl;
        cout<<"RC Grid Generator"<<endl;
        cout<<"================="<<endl;
        parse_argument_list( argc, argv );
    }
    else
    {
        cout<<endl;
        cout<<"RLC Grid Generator"<<endl;
        cout<<"=================="<<endl;
        parse_argument_list( argc, argv );
    }
	
	using namespace global_vars;
	// check if all inputs are provided
	if ( config_file.empty() )
	{
		fprintf(stderr, "Cannot proceed without config/options file. Exiting!!\n");
		return 0;
	}
	if ( sp_file.empty() && !print_grid_info )
	{
		print_grid_info = 1;
	}

	cout<<endl;
	power_grid pg( config_file );
	
	if ( print_grid_info )
		pg.print_power_grid_info();
	
	fprintf(stdout, "\n");

    // Create G and C
	pg.make_G_matrix();
	pg.make_C_matrix();
	matrix:: sparse_matrix &G = pg.get_G_matrix_ref();
	matrix:: sparse_matrix &C = pg.get_C_matrix_ref();

	// calculate voltage drops
	vector<double> Is( pg.get_number_of_nodes() );
	vector<current_source*> pg_cs = pg.get_list_of_current_sources();
	for( unsigned int i = 0; i < pg_cs.size(); i++ )
	{
		Is[ pg_cs[i]->get_n0()->get_index()-1 ] = pg_cs[i]->get_value();
	}

	vector<double> v0;
	{
		profiler record1("Calculate v0 using cholmod");
		matrix::chol_factorize(G);
		v0 = matrix::chol_solve(G, Is);
	}
	
	FILE* mat = fopen("v0.m", "w");
	print_vector_matlab(mat, v0, "v");
	fclose(mat);
	
	// find result_file statistical details about intial voltage drops
	double max_vdrop = 0.0, min_vdrop = pg.get_vdd(), avg_vdrop = 0.0, stddev_vdrop = 0.0;
	for ( unsigned int i = 0; i < v0.size(); i++ )
	{
		if ( v0[i] > max_vdrop )			max_vdrop = v0[i];
		else if ( v0[i] < min_vdrop )		min_vdrop = v0[i];
		avg_vdrop    += v0[i];
		stddev_vdrop += v0[i]*v0[i];
	}
	if ( v0.size() > 0 )
	{
		avg_vdrop    /= v0.size();
		stddev_vdrop /= v0.size();
		stddev_vdrop = sqrt( stddev_vdrop - avg_vdrop*avg_vdrop );
	}
	fprintf(stdout, "Voltage Drop:\n  - Max = %g (%g %%),\n  - Min = %g (%g %%),"
		"\n  - Avg = %.5f,\n  - Std = %g\n", max_vdrop, 100*max_vdrop/pg.get_vdd(),
		min_vdrop, 100*min_vdrop/pg.get_vdd(), avg_vdrop, stddev_vdrop );

	// calculate total grid (logic) power dissipation
	double tpd = 0;
	for( unsigned int i = 0; i < pg_cs.size(); i++ )
	{
		int cs_idx = pg_cs[i]->get_n0()->get_index()-1;
		tpd += ( pg.get_vdd() - v0[cs_idx] )*Is[ cs_idx ];
	}
	fprintf(stdout, "Total Logic Power Dissipation = %g W\n", tpd);
	fprintf(stdout, "Overall Power density = %g W/cm^2\n",
		tpd/(pg.get_grid_dimension_x()*pg.get_grid_dimension_y()*1e-8));

	if ( !sp_file.empty() )
	{
		cerr << endl << "Generating spice file ...." << endl;
		pg.generate_spice_file_for_resistive_grid( sp_file );
		cerr << "Write spile file " << sp_file <<endl;
	}

	FILE* out = stdout;
	profiler::print_all_task_timer_records( out );
	profiler::clear_all_task_timer_records();

	// get and print resource usage
	char buffer[2000];
	profiler::get_resource_usage( buffer, RUSAGE_SELF );
	printf("\n%s", buffer);

    // Get current sources indices (zero-based)
    vector<int> current_sources_indices;
    vector<current_source* > cs_vector = pg.get_list_of_current_sources();
    for (unsigned int i = 0; i<pg.get_number_of_current_sources(); i++)
        current_sources_indices.push_back(cs_vector[i]->get_n0()->get_index() - 1);
    sort(current_sources_indices.begin(), current_sources_indices.end());

    //Number of nodes and sources
    int _number_of_nodes   = G.get_number_of_rows();
    int _number_of_sources = current_sources_indices.size();
    double _number_of_c4   = pg.get_list_of_c4_pads().size();

    stringstream gridname, createfolder, ss_g, ss_c, ss_l, ss_m, ss_csi, copyoptions;

    if (!RLC)
    {
        string destination = "//autofs//fs1.ece//fs1.eecg.najm//r//r0//fawazmoh//Desktop//rc_power_grids//";
        double xdim   = pg.get_grid_dimension_x();
        double ydim   = pg.get_grid_dimension_y();
        int    layers = pg.get_list_of_layers().size();
        gridname<<"RC"<<"_"
                <<xdim<<"_"
                <<ydim<<"_"
                <<layers<<"_"
                <<pg.get_number_of_layers_global_grid()<<"_"
                <<pg.get_number_of_sub_grids_X_DIM()<<"_"
                <<pg.get_number_of_sub_grids_X_DIM();

        // Create a folder for the grid and copy the config file
        createfolder<<"mkdir -p /autofs/fs1.ece/fs1.eecg.najm/r/r0/fawazmoh/Desktop/rc_power_grids/"<<gridname.str();
        copyoptions<<"cp "<<config_file<<" /autofs/fs1.ece/fs1.eecg.najm/r/r0/fawazmoh/Desktop/rc_power_grids/"<<gridname.str()<<"//config_"  <<gridname.str();
        system(createfolder.str().c_str());
        system(copyoptions.str().c_str());
        
        ss_g  <<destination<<gridname.str()<<"//g_"  <<gridname.str();
        ss_c  <<destination<<gridname.str()<<"//c_"  <<gridname.str();
        ss_csi<<destination<<gridname.str()<<"//csi_"<<gridname.str();
        print_triplet_to_file(G, ss_g.str());
        print_triplet_to_file(C, ss_c.str());

        //Print current source indices
        FILE * fptr = fopen(ss_csi.str().c_str(), "w");
        for (unsigned int i = 0; i<current_sources_indices.size(); i++)
            fprintf(fptr, "%d\n", current_sources_indices[i]);

        fclose(fptr);
    }
    else
    {
        int _number_of_inductors = _number_of_c4;

        // Rebuild G
        // ASSUME C4 RESISTANCE IS 0.05 OHMS
        G.get_triplet_format();
        G.set_number_of_columns(_number_of_nodes + _number_of_inductors);
        G.set_number_of_rows(_number_of_nodes + _number_of_inductors);
        for (int j = _number_of_nodes; j<_number_of_nodes+_number_of_inductors; j++)
        {
            int i = pg.get_list_of_c4_pads().at(j-_number_of_nodes)->get_port_node_ptr()->get_index()-1;
            G.insert_element(j,j,  1.0/pg.get_C4_resistance());
            G.insert_element(j,i, -1.0/pg.get_C4_resistance());
            G.insert_element(i,j, -1.0/pg.get_C4_resistance());
        }
        G.set_nz(G.get_nz()+3*_number_of_inductors);
        G.get_column_format();

        // Rebuild C
        C.get_triplet_format();
        C.set_number_of_columns(_number_of_nodes + _number_of_inductors);
        C.set_number_of_rows(_number_of_nodes + _number_of_inductors);
        for (int j = _number_of_nodes; j<_number_of_nodes+_number_of_inductors; j++)
        {
            C.insert_element(j,j,  1e-13 + 
               (1e-12 - 1e-13)*((double)rand()/RAND_MAX));
        }
        C.set_nz(C.get_nz()+_number_of_inductors);
        C.get_column_format();

        // build L
        sparse_matrix L(_number_of_inductors, _number_of_inductors, _number_of_inductors);
        for (int i = 0; i<_number_of_inductors; i++)
            L.insert_element(i,i,pg.get_C4_inductance());
        L.get_column_format();
    
        // build M
        sparse_matrix M(_number_of_nodes+_number_of_inductors, _number_of_inductors, _number_of_inductors);
        for (int j =0 ; j<_number_of_inductors; j++)
        {
            int i = _number_of_nodes + j; 
            M.insert_element(i, j, 1);
        }
        M.get_column_format();

        string destination = "//autofs//fs1.ece//fs1.eecg.najm//r//r0//fawazmoh//Desktop//rlc_power_grids//";
        
        double xdim   = pg.get_grid_dimension_x();
        double ydim   = pg.get_grid_dimension_y();
        int    layers = pg.get_list_of_layers().size();
        gridname<<"RLC"<<"_"
                <<xdim<<"_"
                <<ydim<<"_"
                <<layers<<"_"
                <<pg.get_number_of_layers_global_grid()<<"_"
                <<pg.get_number_of_sub_grids_X_DIM()<<"_"
                <<pg.get_number_of_sub_grids_X_DIM();
       
        // Create a folder for the grid and copy the config file
        createfolder<<"mkdir -p /autofs/fs1.ece/fs1.eecg.najm/r/r0/fawazmoh/Desktop/rlc_power_grids/"<<gridname.str();
        copyoptions<<"cp "<<config_file<<" /autofs/fs1.ece/fs1.eecg.najm/r/r0//fawazmoh/Desktop/rlc_power_grids/"<<gridname.str()<<"//config_"  <<gridname.str();
        system(createfolder.str().c_str());
        system(copyoptions.str().c_str());
        
        ss_g  <<destination<<gridname.str()<<"//g_"  <<gridname.str();
        ss_c  <<destination<<gridname.str()<<"//c_"  <<gridname.str();
        ss_l  <<destination<<gridname.str()<<"//l_"  <<gridname.str();
        ss_m  <<destination<<gridname.str()<<"//m_"  <<gridname.str();
        ss_csi<<destination<<gridname.str()<<"//csi_"<<gridname.str();
        print_triplet_to_file(G, ss_g.str());
        print_triplet_to_file(C, ss_c.str());
        print_triplet_to_file(L, ss_l.str());
        print_nonsym_triplet_to_file(M, ss_m.str());

        //Print current source indices
        FILE * fptr = fopen(ss_csi.str().c_str(), "w");
        for (unsigned int i = 0; i<current_sources_indices.size(); i++)
            fprintf(fptr, "%d\n", current_sources_indices[i]);

        fclose(fptr);

        cout<<endl<<"After alterations, an RLC grid is generated with the following (updated) paramters: "<<endl;
        cout << setw(40) << "Total number of nodes: "
             << _number_of_inductors +_number_of_nodes << endl;
        cout << setw(40) << "Total number of inductors: "
             << _number_of_inductors << endl;
    }

	pg.deallocate_memory();
	G.destroy_chol_L();

	return 0;
}

void print_matrix_matlab(FILE* matlab, sparse_matrix M, string M_name)
{   
    stringstream Ti, Tj, Tx;
    Ti << M_name << "Ti";
    Tj << M_name << "Tj";
    Tx << M_name << "Tx";
    string sTi, sTj, sTx;
    sTi = Ti.str();
    sTj = Tj.str();
    sTx = Tx.str();
    print_vector_matlab(matlab, M.get_Ti(), sTi);
    print_vector_matlab(matlab, M.get_Tj(), sTj);
    print_vector_matlab(matlab, M.get_Tx(), sTx);
    fprintf(matlab, "%s = triplet_format_to_matrix(%s, %s, %s, %d, %d);\n", 
            M_name.c_str(), sTi.c_str(), sTj.c_str(), sTx.c_str(),
            M.get_number_of_rows(), M.get_number_of_columns());
    fprintf(matlab, "clear %s %s %s \n", sTi.c_str(), sTj.c_str(),
                    sTx.c_str());
}

void print_vector_matlab(FILE* matlab, vector<double> v, string v_name)
{
	fprintf(matlab, "%s = [", v_name.c_str());
	for (unsigned int i = 0; i < v.size() - 1; i++)
	{
		fprintf(matlab, "%-15g;", v[i]);
		if ((i + 1) % 10 == 0)
			fprintf(matlab, "\n");
	}
	fprintf(matlab, "%-15g", v[v.size()-1]);
	fprintf(matlab, "];\n\n");
}
void print_vector_matlab(FILE* matlab, vector<int> v, string v_name)
{
	fprintf(matlab, "%s = [", v_name.c_str());
	for (unsigned int i = 0; i < v.size() - 1; i++)
	{
		fprintf(matlab, "%-15d;", v[i]);
		if ((i + 1) % 10 == 0)
			fprintf(matlab, "\n");
	}
	fprintf(matlab, "%-15d", v[v.size()-1]);
	fprintf(matlab, "];\n\n");
}

void parse_argument_list (int argc, char** argv)
{
	// set the default values for the global vars
	global_vars::verbose = 0;
	global_vars::print_grid_info = 0;
	
	int c;

	struct option long_options[] =
	{
	  /* These options set a flag. */
	  {"verbose", no_argument,       &global_vars::verbose, 1},
	  {"info",    no_argument,       &global_vars::print_grid_info, 1},
	  /* These options don’t set a flag.
		 We distinguish them by their indices. */
	  {"config",  required_argument,       0, 'c'},
	  {"spfile",  required_argument,       0, 's'},
	  {0, 0, 0, 0}
	};


	while (1)
	{
	  	/* getopt_long stores the option index here. */
		int option_index = 0;

		c = getopt_long (argc, argv, "vic:s:", long_options, &option_index);

		/* Detect the end of the options. */
		if (c == -1)
			break;

	switch (c)
	{
		case 0:
			/* If this option set a flag, do nothing else now. */
			if (long_options[option_index].flag != 0)
				break;
			/*printf ("option %s", long_options[option_index].name);
			if (optarg)
				printf (" with arg %s", optarg);
			printf ("\n");
			break;*/

		case 'c':
			global_vars::config_file.assign(optarg);
			break;

		case 's':
			global_vars::sp_file.assign(optarg);
			break;

		case '?':
		  /* getopt_long already printed an error message. */
		  break;

		default:
			abort ();
		}
	}

	/* Print any remaining command line arguments (not options). */
	if (optind < argc)
	{
		printf ("non-option ARGV-elements: ");
		while (optind < argc)
			printf ("%s ", argv[optind++]);
		putchar ('\n');
	}
}

void print_triplet_to_file(sparse_matrix &A, string filename)
{
    A.get_triplet_format();
    FILE * fptr = fopen(filename.c_str(), "w");
    
    vector<int>    Ti = A.get_Ti();
    vector<int>    Tj = A.get_Tj();
    vector<double> Tx = A.get_Tx();
    fprintf(fptr, "%d %d %d\n", A.get_number_of_rows(), A.get_number_of_columns(), A.get_nz());    
    for (unsigned int i = 0; i<Ti.size(); i++)
    {
        if (Ti[i] > Tj[i])
        {    
            fprintf(fptr, "%d %d %.26f\n", Ti[i], Tj[i], Tx[i]);
            fprintf(fptr, "%d %d %.26f\n", Tj[i], Ti[i], Tx[i]);
        }
        else if (Ti[i] == Tj[i])
            fprintf(fptr, "%d %d %.26f\n", Ti[i], Tj[i], Tx[i]);


    }
    fclose(fptr);
}

void print_nonsym_triplet_to_file(sparse_matrix &A, string filename)
{
    A.get_triplet_format();
    FILE * fptr = fopen(filename.c_str(), "w");
    
    vector<int>    Ti = A.get_Ti();
    vector<int>    Tj = A.get_Tj();
    vector<double> Tx = A.get_Tx();
    fprintf(fptr, "%d %d %d\n", A.get_number_of_rows(), A.get_number_of_columns(), A.get_nz());    
    for (unsigned int i = 0; i<Ti.size(); i++)
    {
        fprintf(fptr, "%d %d %.26f\n", Ti[i], Tj[i], Tx[i]);
    }
    fclose(fptr);
}
