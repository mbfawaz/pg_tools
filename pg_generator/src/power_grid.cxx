/**
* \file power_grid.cxx \brief Power Grid Implementation
*
* This file implements all the methods of the main class (power_grid) and builds the power grid
*/
#include "power_grid.h"
#include "profiler.h"
#include <algorithm>
using namespace std;
using namespace matrix;

// for facilitating profiler task descriptions
static string task_desc;

static bool comp_ycoord(node* n1, node* n2);
static bool comp_xcoord(node* n1, node* n2);

namespace common
{
	vector<bool> random_mask_vector(int size, int num_of_ones)
	{
		vector<bool> vec;
		vec.resize(size,0);
		for(int i = 0; i < num_of_ones ; i++)
		{
			vec[i] = 1;
		}
		random_shuffle(vec.begin(),vec.end());

		return vec;
	}

	vector<double> random_normalized_vector(unsigned int size)
	{
		vector<double> rnv( size );
		double sum = 0;
		for ( unsigned int i = 0; i < size; i++ )
		{
			rnv[i] = 1.0*rand();
			sum += rnv[i];
		}
		for ( unsigned int i = 0; i < size; i++ )
			rnv[i] /= sum;
		return rnv;
	}
}

namespace global_vars
{
	int verbose;
	int print_grid_info;
	string config_file;
	string sp_file;
}

namespace rand_options
{
	bool _randomize_sub_grids;
	double _ratio_of_ON_sub_grids;
	double _ratio_of_retained_trees_in_sub_grids;
	double _ratio_of_retained_trees_in_intermediate_layers;
	double _ratio_of_retained_trees_in_global_layers;
	double _max_shift_tree_spans;
	int _seed_offset = 0;
}

// Power_grid Class

//! Default Constructor
power_grid::power_grid()
{
	_number_of_nodes = 0;
	_number_of_c4s = 0;
	_number_of_sub_grids = 0;
	_number_of_branches = 0;
	_number_of_layers = 0;
	_number_of_layers_global_grid = 0;
	_number_of_current_sources = 0;
	_unit = 1e-6;
	_is_memory_deallocated = false;
	_power_scale_factor = 1;
}

//! Regular Constructor
//! @param params_file [in]: name of the file to read the parameters from
power_grid::power_grid(string params_file)
{
	{
		profiler record("Building Power Grid");

		read_config_file(params_file);
		
		_delta_t = pow(10,-12);
		_is_memory_deallocated = false;

		// Initialize the number of subgrids. This will depend on the ratio of
		// ON subgrids specified by the user
		srand(rand_options::_seed_offset+1);
		_number_of_sub_grids = floor(_number_of_sub_grids_X_DIM*_number_of_sub_grids_Y_DIM*
			rand_options::_ratio_of_ON_sub_grids);
		if ( _number_of_sub_grids == 0 )
		{
			printf("\033[1;%dmERROR: None of the sub-grids are ON. Exiting Now.\033[0m\n", 33);
			exit(1);
		}

		// Generate a random boolean vector that specifies which subgrid is ON and which is OFF
		_ON_sub_grids = common::random_mask_vector(_number_of_sub_grids_X_DIM*_number_of_sub_grids_Y_DIM,
			_number_of_sub_grids);

		// resize the list of sub grids
		_list_of_sub_grids.resize(_number_of_sub_grids);

		// Make the ground and vdd node
		_ground_node = new node(NULL, "n0", 0, -1, -1, CHIP_NODE);
		_vdd_node    = new node(NULL, "n_vdd", -1, -1, -1, CHIP_NODE);
		
		_number_of_layers_sub_grid = _number_of_layers - _number_of_layers_global_grid; 
		_lowest_layer_index = 1;
		_number_of_current_sources = 0;
		_power_scale_factor = 1;

		string msg;

		create_global_grid();
		create_internal_grid_islands();
		cerr<<endl;
		int res_idx = 0;
		task_desc.assign("Create Port Nodes & Port Vias");
		if ( create_port_nodes_and_port_vias( res_idx ) )
		{
			if ( global_vars::verbose )
				cerr<<endl<<task_desc<<string( 60-task_desc.length(), ' ' )<<"done"<<" ("
					<<profiler::get_task_timer_data_by_desc( task_desc )->elapsed.wall/1e9<<" secs)";

			/*task_desc = "Add New Nodes for C4 Layer";
			add_new_nodes_for_C4();
			if ( global_vars::verbose )
				cerr<<endl<<task_desc<<string( 60-task_desc.length(), ' ' )<<"done"<<" ("
					<<profiler::get_task_timer_data_by_desc( task_desc )->elapsed.wall/1e9<<" secs)";//*/

			task_desc.assign("Sort Grid Nodes");
			sort_nodes_and_set_indices();
			if ( global_vars::verbose )
				cerr<<endl<<task_desc<<string( 60-task_desc.length(), ' ' )<<"done"<<" ("
					<<profiler::get_task_timer_data_by_desc( task_desc )->elapsed.wall/1e9<<" secs)";

			task_desc.assign("Add Branch Resistances");
			add_branch_resistances( res_idx );
			if ( global_vars::verbose )
				cerr<<endl<<task_desc<<string( 60-task_desc.length(), ' ' )<<"done"<<" ("
					<<profiler::get_task_timer_data_by_desc( task_desc )->elapsed.wall/1e9<<" secs)";

			int cap_idx = 0;
			task_desc.assign("Add Capacitors");
			add_capacitors( cap_idx );
			if ( global_vars::verbose )
				cerr<<endl<<task_desc<<string( 60-task_desc.length(), ' ' )<<"done"<<" ("
					<<profiler::get_task_timer_data_by_desc( task_desc )->elapsed.wall/1e9<<" secs)";

			task_desc.assign("Add C4 Pads");
			add_c4_pads( res_idx, cap_idx );
			if ( global_vars::verbose )
				cerr<<endl<<task_desc<<string( 60-task_desc.length(), ' ' )<<"done"<<" ("
					<<profiler::get_task_timer_data_by_desc( task_desc )->elapsed.wall/1e9<<" secs)";//*/

			task_desc.assign("Finalize Node Names & Indices");
			finalize_node_name_and_indices();
			if ( global_vars::verbose )
				cerr<<endl<<task_desc<<string( 60-task_desc.length(), ' ' )<<"done"<<" ("
					<<profiler::get_task_timer_data_by_desc( task_desc )->elapsed.wall/1e9<<" secs)";
		}
		else
		{
			cerr<<endl<<task_desc<<string( 60-task_desc.length(), ' ' )<<"failed"<<" ("
				<<profiler::get_task_timer_data_by_desc( task_desc )->elapsed.wall/1e9<<" secs)";
			return;
		}
		
	}
}

//! Destructor of the power grid
power_grid::~power_grid()
{
	if ( !_is_memory_deallocated )
		fprintf(stderr, "WARNING!!! Function deallocate_memory() has not been called."
			"This will result in memory leaks.\n");
}

void power_grid::deallocate_memory()
{
	// delete all layer specifications
	for (unsigned int i = 0; i < _list_of_layer_specs.size(); i++)
		delete _list_of_layer_specs[i];

	// delete global grid
	delete _global_grid;

	// Delete the sub-grids. 
	for (unsigned int i = 0; i < _list_of_sub_grids.size(); i++)
		delete _list_of_sub_grids[i];

	// Deleting all layers (GLOBAL AND INTERNAL). Deleting the layers 
	// deletes the associated interconnect trees, which in turn deletes
	// its nodes and branches. 
	for (unsigned int i = 0; i < _number_of_layers; i++)
		for(unsigned int j = 0; j < _list_of_layers[i].size(); j++)
			delete _list_of_layers[i].at(j);

	// Deleting C4_pads
	for (unsigned int i = 0; i < _list_of_c4_pads.size(); i++)
		delete _list_of_c4_pads[i];

	// Deleting vias
	for (unsigned int i = 0; i < _list_of_vias.size(); i++)
		delete _list_of_vias[i];

	// Deleting current sources
	//for (unsigned int i = 0; i < _list_of_current_sources.size(); i++)
	//	delete _list_of_current_sources[i];

	// Deleting ground node
	if ( _ground_node != NULL )
	{
		delete _ground_node;
		_ground_node = NULL;
	}

	// Deleting Vdd node
	if ( _vdd_node != NULL )
	{
		delete _vdd_node;
		_vdd_node = NULL;
	}

	_is_memory_deallocated = true;
}

void power_grid::create_global_grid()
{
	// Create the Global Grid: the metal layers, the interconnect trees, the nodes and via resistors.
	// fprintf(stdout, "\033[1;%dmCreating Global Grid\033[0m\n",36);
	_global_grid = new sub_grid(this, GLOBAL, 0, 0, _GRID_DIMENSION_X_UM, 0, _GRID_DIMENSION_Y_UM);
	_number_of_nodes = _global_grid->get_number_of_nodes(); 
	
	// Check if the number of C4s is very small
	if (_list_of_layers.at(0).at(0)->get_list_of_interconnect_trees().size() < 15)
	{
		fprintf(stdout, "\033[;%dm\n   WARNING: Number of interconnect trees in the C4 layer is small!\n"
		"   This may result in small number of C4 pads.\033[0m\n",93);
	}
}

void power_grid::create_internal_grid_islands()
{
	// Create the Internal Grids: the metal layers, the interconnect trees, the nodes and
	// via resistors, and the current sources.
	double x_start, x_final, y_start, y_final;
	
	int first_layer_index = get_number_of_layers() - get_number_of_layers_global_grid() + 1;
	unsigned int index = 1;
	double pitch_vertical   =  _list_of_layer_specs.at(_number_of_layers-first_layer_index+1)->pitch;
	double pitch_horizontal =  _list_of_layer_specs.at(_number_of_layers-first_layer_index+2)->pitch;
	if(_number_of_sub_grids_X_DIM == 1 && _number_of_sub_grids_Y_DIM == 1)
	{
		pitch_vertical = 0;
		pitch_horizontal = 0;
	}
	
	// genarate normalized power scale factors for all subgrids
	srand(rand_options::_seed_offset+10);
	vector<double> psf_vec = common::random_normalized_vector( _number_of_sub_grids );

	int total_sub_grids = _number_of_sub_grids_Y_DIM*_number_of_sub_grids_X_DIM;
	int psf_vec_idx = 0;
	for (unsigned int i = 0; i < _number_of_sub_grids_Y_DIM; i++)
	{
		for (unsigned int j = 0; j < _number_of_sub_grids_X_DIM; j++)
		{
			if(_ON_sub_grids[_number_of_sub_grids_X_DIM*i+j] == 0)
				continue;

			sub_grid * sg;
			x_start = (_GRID_DIMENSION_X_UM/_number_of_sub_grids_X_DIM)*j;
			y_start = (_GRID_DIMENSION_Y_UM/_number_of_sub_grids_Y_DIM)*i; 
			
			if(_list_of_layer_specs.at(_number_of_layers - first_layer_index + 1)->
				orientation == VERTICAL)
			{
				x_final = (_GRID_DIMENSION_X_UM/_number_of_sub_grids_X_DIM)*(j+1) - pitch_vertical;
				y_final = (_GRID_DIMENSION_Y_UM/_number_of_sub_grids_Y_DIM)*(i+1) - pitch_horizontal;
			}
			else
			{
				x_final = (_GRID_DIMENSION_X_UM/_number_of_sub_grids_X_DIM)*(j+1) - pitch_horizontal;
				y_final = (_GRID_DIMENSION_Y_UM/_number_of_sub_grids_Y_DIM)*(i+1) - pitch_vertical;
			}
			sg = new sub_grid(this,INTERNAL, index, x_start, x_final, y_start, y_final,
				total_sub_grids*psf_vec[psf_vec_idx++] );
			_number_of_nodes += sg->get_number_of_nodes();
			sg->set_index(index);
			index++;
			add_sub_grid(sg);
		}
	}
}

bool power_grid::create_port_nodes_and_port_vias( int &res_idx )
{
	// Creating Nodes for the port region, i.e., between the global grid and the internal grids.
	profiler record( task_desc );
	int node_idx = 0;
	string node_name, via_name;
	node *n1, *n2;
	resistor * via;
	double xcoord, ycoord;
	bool port_via = 0;
	node_idx = _number_of_nodes;

	// get the bottom global layer
	layer* btm_gl_layer =
		get_layer(_number_of_layers - _number_of_layers_global_grid + 1).at(0);
	double via_resistance  = btm_gl_layer->get_via_resistance();
	vector<interconnect_tree*> btm_gl_layer_itree_list =
		btm_gl_layer->get_list_of_interconnect_trees();

	for (unsigned int i = 0; i < _number_of_sub_grids; i++)
	{
		port_via = 0;
		sub_grid* this_sg = _list_of_sub_grids.at(i);
		
		// get the top layer for this subgrid
		unsigned int sub_num_layers = this_sg->get_number_of_layers();
		layer* top_sg_layer = this_sg->get_layer(sub_num_layers).at(0);

		vector<interconnect_tree*> top_sg_layer_itree_list = top_sg_layer->
			get_list_of_interconnect_trees();

		if ( top_sg_layer->get_orientation() == HORIZONTAL) // PUT SOMETHING INSTEAD OF 0 !!
		{
			//cout<<"Going into if "<<top_sg_layer_itree_list.size()<<endl;
			for (unsigned int j = 0; j < top_sg_layer_itree_list.size(); j++)
			{
				//cout<<"xstart = "<<this_sg->get_x_start()<<" xfinal = "<<this_sg->get_x_final()<<endl;
				ycoord = top_sg_layer_itree_list.at(j)->get_ycoord();

				// introduce randomness by staring and ending at arbitrary trees
				int start_k = 0, end_k = btm_gl_layer_itree_list.size();
				int btm_gl_layer_itree_list_size = btm_gl_layer_itree_list.size();
				int start_k_margin = (int)( rand_options::_max_shift_tree_spans*
					btm_gl_layer_itree_list_size );
				if ( start_k_margin != 0 )
				{
					start_k = rand() % start_k_margin;
					int end_k_margin = (int)( rand_options::_max_shift_tree_spans*
						btm_gl_layer_itree_list_size );
					end_k = btm_gl_layer_itree_list_size - ( rand() % end_k_margin );
					assert( end_k - start_k >= 1 );
				}

				for (int k = start_k; k < end_k; k++)
				{
					xcoord = btm_gl_layer_itree_list.at(k)->get_xcoord();
					
					// Checking if interconnect tree of the global grid is in range.
					// The range is defined by the xstart, ystart, xfinal, yfinal of the internal grid
					if (xcoord < this_sg->get_x_start() ||
						xcoord > this_sg->get_x_final())
						continue;
					
					port_via = 1;
					res_idx++;
					stringstream ss;
					ss << "R" << res_idx << "_prt";
					via_name = ss.str();
					via = new resistor(via_name, res_idx, via_resistance, -1, -1, VIA);

					node_idx++;
					stringstream ss1;
					ss1 << "n" << node_idx;
					node_name = ss1.str();
					n1 = new node(top_sg_layer_itree_list.at(j),
						node_name, node_idx, xcoord, ycoord, CHIP_NODE);
					
					node_idx++;
					stringstream ss2;
					ss2 << "n" << node_idx;
					node_name = ss2.str();
					n2 = new node(btm_gl_layer_itree_list.at(k),
						node_name, node_idx, xcoord, ycoord, CHIP_NODE);

					// Adding Via resistors
					via->set_nodes(n1, n2);
					
					n1->add_resistor(via); // Add node n1 to the internal grids
					top_sg_layer_itree_list.at(j)->add_node(n1);
					top_sg_layer_itree_list.at(j)->add_node_to_map(xcoord, n1);
					this_sg->increment_number_of_nodes();

					n2->add_resistor(via); // Add node n2 to the global grid
					btm_gl_layer_itree_list.at(k)->add_node(n2);
					btm_gl_layer_itree_list.at(k)->add_node_to_map(ycoord, n2);
					_global_grid->increment_number_of_nodes();

					_list_of_vias.push_back(via);
					this_sg->add_via(via);
				}
			}
		}
		else
		{
			for (unsigned int j = 0; j < top_sg_layer_itree_list.size(); j++)
			{
				xcoord = top_sg_layer_itree_list.at(j)->get_xcoord();
				
				// introduce randomness by staring and ending at arbitrary trees
				int start_k = 0, end_k = btm_gl_layer_itree_list.size();
				int btm_gl_layer_itree_list_size = btm_gl_layer_itree_list.size();
				int start_k_margin = (int)( rand_options::_max_shift_tree_spans*
					btm_gl_layer_itree_list_size );
				if ( start_k_margin != 0 )
				{
					start_k = rand() % start_k_margin;
					int end_k_margin = (int)( rand_options::_max_shift_tree_spans*
						btm_gl_layer_itree_list_size );
					end_k = btm_gl_layer_itree_list_size - ( rand() % end_k_margin );
					assert( end_k - start_k >= 1 );
				}

				for (int k = start_k; k < end_k; k++)
				{
					ycoord = btm_gl_layer_itree_list.at(k)->get_ycoord();
					
					//cout<<this_sg->get_y_start()<<" "<<this_sg->get_y_final()<<endl; 
					// Checking if interconnect tree of the global grid is in range.
					// The range is defined by the xstart, ystart, xfinal, yfinal of the internal grid
					if (ycoord < this_sg->get_y_start() || 
						ycoord > this_sg->get_y_final())
						continue;
					port_via = 1;
					res_idx++;
					stringstream ss;
					ss << "R" << res_idx;
					via_name = ss.str();
					via = new resistor(via_name, res_idx, via_resistance, -1, -1, VIA);

					node_idx++;
					stringstream ss1;
					ss1 << "n" << node_idx;
					node_name = ss1.str();
					n1 = new node(top_sg_layer_itree_list.at(j),
					node_name, node_idx, xcoord, ycoord, CHIP_NODE);

					node_idx++;
					stringstream ss2;
					ss2 << "n" << node_idx;
					node_name = ss2.str();
					n2 = new node( btm_gl_layer_itree_list.at(k),
					node_name, node_idx, xcoord, ycoord, CHIP_NODE);

					// Adding Via resistors
					via->set_nodes(n1, n2);
					
					n1->add_resistor(via);
					top_sg_layer_itree_list.at(j)->add_node(n1);
					top_sg_layer_itree_list.at(j)->add_node_to_map(ycoord, n1);
					this_sg->increment_number_of_nodes();

					n2->add_resistor(via);
					btm_gl_layer_itree_list.at(k)->add_node(n2);
					btm_gl_layer_itree_list.at(k)->add_node_to_map(xcoord, n2);
					_global_grid->increment_number_of_nodes();

					_list_of_vias.push_back(via);
					this_sg->add_via(via);
				}
			}
		}
		if(port_via == 0)
		{
			cout<<"One or more internal grids are disconnected from the global grid. Exiting.."<<endl;
			return 0;
		}
	}
	_number_of_nodes = node_idx; 
	
	return 1;
}

void power_grid::add_new_nodes_for_C4()
{
	// Creating C4s
	profiler record( task_desc );

	int num_nodes_added = 0;

	layer* c4_layer = get_c4_layer();
	double c4_distance = get_C4_spacing_factor()*c4_layer->get_pitch();
	bool is_horizontal = (c4_layer->get_orientation() == HORIZONTAL);
	double stop_it = (is_horizontal)?_GRID_DIMENSION_Y_UM:_GRID_DIMENSION_X_UM;

	map<double, interconnect_tree*> c4_itree_map = c4_layer->get_it_map();
	map<double, interconnect_tree*>::iterator itr;

	int node_idx = 0;

	for (double dist_it = c4_layer->get_offset(); dist_it <= stop_it;
		dist_it += c4_distance)
	{
		itr = c4_itree_map.find( dist_it );
 		if (itr == c4_itree_map.end())		continue;

		interconnect_tree* it = itr->second;
		vector<node*> itree_node_list = it->get_list_of_nodes();
		double off = is_horizontal?
			itree_node_list.at(0)->get_xcoord()
			:itree_node_list.at(0)->get_ycoord();
		
		double stop_node = is_horizontal?
			itree_node_list.at(itree_node_list.size()-1)->get_xcoord()
			:itree_node_list.at(itree_node_list.size()-1)->get_ycoord();
		
		// This code creates new C4 nodes, if required. Problem is, that now it has to be called
		// before assigning branches.
		for (double dist_node = off; dist_node < stop_node; dist_node += c4_distance)
		{
			// check if a node already exists at this point
			node* c4_node = NULL;
			if ( is_horizontal )
			{
				for( unsigned int i = 0; i < itree_node_list.size(); i++ )
				{
					if ( abs( itree_node_list[i]->get_xcoord() - dist_node ) < 3 )
					{
						c4_node = itree_node_list[i]; break;
					}
				}
			}
			else
			{
				for( unsigned int i = 0; i < itree_node_list.size(); i++ )
					if ( abs( itree_node_list[i]->get_ycoord() - dist_node) < 3 )
					{
						c4_node = itree_node_list[i]; break;
					}
			}
			
			// if not, then create it
			if ( c4_node == NULL )
			{
				double xcoord = (is_horizontal) ? dist_it : it->get_xcoord();
				double ycoord = (is_horizontal) ? it->get_ycoord() : dist_it;
				
				node_idx++;
				stringstream ss1;
				ss1 << "n" << node_idx;
				string node_name = ss1.str();
				c4_node = new node( it, node_name, node_idx, xcoord, ycoord, CHIP_NODE);
				it->add_node(c4_node);
				if ( is_horizontal )
					it->add_node_to_map(xcoord, c4_node);
				else
					it->add_node_to_map(ycoord, c4_node);
				_global_grid->increment_number_of_nodes();
				num_nodes_added++;
			}//*/
		}
	}
	printf("\n*** Num nodes added = %d ***\n", num_nodes_added);
	_number_of_nodes += num_nodes_added;
}

void power_grid::sort_nodes_and_set_indices()
{
	profiler record( task_desc );
	unsigned int idx = 1;
	
	// Start by sorting the lowest layer in the global grid up to the highest metal layer
	for (unsigned int i = 0; i < _number_of_layers_global_grid ; i++)
	{
		layer* gl_layer = _list_of_layers.at(i).at(0);
		vector<interconnect_tree*> gl_layer_itree_list = gl_layer->get_list_of_interconnect_trees_ref();
		for (unsigned int j = 0; j < gl_layer_itree_list.size(); j++)
		{
			if ( gl_layer->get_orientation() == HORIZONTAL )
			{
				sort(gl_layer_itree_list.at(j)->get_list_of_nodes_ref().begin(),
					 gl_layer_itree_list.at(j)->get_list_of_nodes_ref().end(),
					comp_xcoord);
			}
			else
			{
				sort(gl_layer_itree_list.at(j)->get_list_of_nodes_ref().begin(),
					 gl_layer_itree_list.at(j)->get_list_of_nodes_ref().end(),
					comp_ycoord);
			}

			for (unsigned int k = 0; k < gl_layer_itree_list.at(j)->get_list_of_nodes_ref().size(); k++)
			{
				gl_layer_itree_list.at(j)->get_list_of_nodes_ref().at(k)->set_index(idx++);
			}
		}
	}

	// Sort the nodes in the internal grids. We start by sorting the nodes in
	// the highest layer going to the lowest layer in that internal grid, then
	// we move to the next internal grid,... etc
	for(unsigned int p = 0; p < _number_of_sub_grids; p++)
	{
		for (unsigned int i = 0; i < _number_of_layers_sub_grid ; i++)
		{
			layer* sg_layer = _list_of_sub_grids.at(p)->get_list_of_layers().at(i).at(0);
			vector<interconnect_tree*> sg_layer_itree_list = sg_layer->get_list_of_interconnect_trees_ref();
			for (unsigned int j = 0; j < sg_layer_itree_list.size(); j++)
			{
				if ( sg_layer->get_orientation() == HORIZONTAL)
				{
					sort(sg_layer_itree_list.at(j)->get_list_of_nodes_ref().begin(),
						sg_layer_itree_list.at(j)->get_list_of_nodes_ref().end(),
						comp_xcoord);
				}
				else
				{
					sort(sg_layer_itree_list.at(j)->get_list_of_nodes_ref().begin(),
						sg_layer_itree_list.at(j)->get_list_of_nodes_ref().end(),
						comp_ycoord);
				}
				for (unsigned int k = 0; k < sg_layer_itree_list.at(j)->get_list_of_nodes_ref().size(); k++)
				{
					sg_layer_itree_list.at(j)->get_list_of_nodes_ref().at(k)->set_index(idx++);
				}
			}
		}
	}
}

void power_grid::add_branch_resistances( int &res_idx )
{
	// Adding Branch Resistors
	profiler record( task_desc );
	double width;
	double length;
	string res_name;
	int num_branches = 0;
	double cs_area;
	double value;
	double sheet_resistance;
	resistor* r;

	for (int i = _number_of_layers_global_grid - 1; i >= 0; i--)
	{
		layer* gl_layer = _list_of_layers.at(i).at(0);
		vector<interconnect_tree*> gl_layer_itree_list_ref =
			gl_layer->get_list_of_interconnect_trees_ref();

		for (unsigned int j = 0; j < gl_layer_itree_list_ref.size(); j++)
		{
			width = gl_layer->get_width();
			cs_area = (width*width)*1.8; // aspect ratio = h/w = 1.8
			sheet_resistance = gl_layer->get_sheet_resistance();

			vector<node*> itree_node_list_ref = gl_layer_itree_list_ref.at(j)->get_list_of_nodes_ref();
			for (unsigned int k = 0; k < itree_node_list_ref.size()-1; k++)
			{
				if (_list_of_layers.at(i).at(0)->get_orientation() == HORIZONTAL)
				{
					length = itree_node_list_ref.at(k+1)->get_xcoord() -
						itree_node_list_ref.at(k)->get_xcoord();
				}
				else
				{
					length = itree_node_list_ref.at(k+1)->get_ycoord() -
						itree_node_list_ref.at(k)->get_ycoord();
				}

				value = (sheet_resistance*length)/width;
				res_idx++;
				num_branches++;
				stringstream ss;
				ss << "R" << res_idx;
				res_name = ss.str();
				r = new resistor(res_name, res_idx, value, length, cs_area, BRANCH);
				r->set_nodes( itree_node_list_ref.at(k), itree_node_list_ref.at(k+1) );
				
				itree_node_list_ref.at(k+1)->add_resistor(r);
				itree_node_list_ref.at(k)->add_resistor(r);
				
				gl_layer_itree_list_ref.at(j)->add_resistor(r);
			}
		}
	}
	_number_of_branches = num_branches;
	_global_grid->set_number_of_branches(num_branches);

	for(unsigned int p = 0; p < _number_of_sub_grids;p++)
	{
		num_branches = 0;
		for (int i = _list_of_sub_grids.at(p)->get_number_of_layers() - 1; i >= 0; i--)
		{
			layer* sg_layer = _list_of_sub_grids.at(p)->get_list_of_layers().at(i).at(0);
			vector<interconnect_tree*> sg_layer_itree_list_ref = sg_layer->get_list_of_interconnect_trees_ref();

			for (unsigned int j = 0; j < sg_layer_itree_list_ref.size(); j++)
			{
				width = sg_layer->get_width();
				cs_area = (width*width)/1.9;
				sheet_resistance = sg_layer->get_sheet_resistance();

				vector<node*> itree_node_list_ref = sg_layer_itree_list_ref.at(j)->get_list_of_nodes_ref();
				for (unsigned int k = 0; k < itree_node_list_ref.size() - 1; k++)
				{
					if (sg_layer->get_orientation() == HORIZONTAL)
					{
						length = itree_node_list_ref.at(k+1)->get_xcoord() - 
								itree_node_list_ref.at(k)->get_xcoord();

						/*if (length*_unit < 1e-7)
						{
							cout << itree_node_list_ref.at(k+1)->get_xcoord() << " " << 
								itree_node_list_ref.at(k)->get_xcoord() << endl;
							cout << itree_node_list_ref.at(k+1)->get_ycoord() << " " << 
								itree_node_list_ref.at(k)->get_ycoord() << endl;
						}//*/
					}
					else
					{
						length = itree_node_list_ref.at(k+1)->get_ycoord() - 
								itree_node_list_ref.at(k)->get_ycoord();

						/*if (length*_unit < 1e-7)
						{
							cout << itree_node_list_ref.at(k+1)->get_xcoord() << " " << 
								itree_node_list_ref.at(k)->get_xcoord() << endl;
							cout << itree_node_list_ref.at(k+1)->get_ycoord() << " " << 
								itree_node_list_ref.at(k)->get_ycoord() << endl;
						}//*/
					}
					/*if (length*_unit < 1e-7)
					{
						cout << "Stripes above each other! Exiting... " << endl;
						exit(0);
					}*/
					
					value = (sheet_resistance*length)/width;
					res_idx++;
					num_branches++;
					stringstream ss;
					ss << "R" << res_idx;
					res_name = ss.str();
					r = new resistor(res_name, res_idx, value, length, cs_area, BRANCH);
					r->set_nodes(itree_node_list_ref.at(k), itree_node_list_ref.at(k+1));

					itree_node_list_ref.at(k+1)->add_resistor(r);
					itree_node_list_ref.at(k)->add_resistor(r);
					sg_layer_itree_list_ref.at(j)->add_resistor(r);
				}
			}
		}
		_list_of_sub_grids.at(p)->set_number_of_branches(num_branches);
		_number_of_branches += num_branches;
	}
}

void power_grid::add_capacitors( int &cap_idx )
{
	srand(rand_options::_seed_offset+5); // re-init the random number generator
	
	// Adding Capacitors
	profiler record( task_desc );
	string cap_name;
	capacitor *c;
	for (int i = _number_of_layers_global_grid - 1; i >= 0; i--)
	{
		vector<interconnect_tree*> gl_layer_itree_list =
			_list_of_layers.at(i).at(0)->get_list_of_interconnect_trees();

		for (unsigned int j = 0; j < gl_layer_itree_list.size(); j++)
		{
			vector<node*> itree_node_list_ref = gl_layer_itree_list.at(j)->get_list_of_nodes_ref();
			for (unsigned int k = 0; k < itree_node_list_ref.size(); k++)
			{
				cap_idx++;
				stringstream ss;
				ss << "C" << cap_idx;
				cap_name = ss.str();
				double value = _MIN_DIE_CAPACITANCE_F + (_MAX_DIE_CAPACITANCE_F -
					_MIN_DIE_CAPACITANCE_F)*((double)rand()/RAND_MAX);
				c = new capacitor(cap_name, cap_idx, value, BRANCH);
				c->set_nodes(itree_node_list_ref.at(k), _ground_node);
				
				itree_node_list_ref.at(k)->add_capacitor(c);
				_ground_node->add_capacitor(c);
			}
		}
	}
		
	for(unsigned int p = 0; p < _number_of_sub_grids ; p++)
	{
		for (int i = _list_of_sub_grids.at(p)->get_number_of_layers() - 1; i >= 0; i--)
		{
			vector<interconnect_tree*> sg_layer_itree_list =
				_list_of_sub_grids.at(p)->get_list_of_layers().at(i).at(0)->get_list_of_interconnect_trees();
			for (unsigned int j = 0; j < sg_layer_itree_list.size(); j++)
			{
				vector<node*> itree_node_list_ref = sg_layer_itree_list.at(j)->get_list_of_nodes_ref();
				for (unsigned int k = 0; k < itree_node_list_ref.size(); k++)
				{
					cap_idx++;
					stringstream ss;
					ss << "C" << cap_idx;
					cap_name = ss.str();
					double value = _MIN_DIE_CAPACITANCE_F + (_MAX_DIE_CAPACITANCE_F -
						_MIN_DIE_CAPACITANCE_F)*((double)rand()/RAND_MAX);
					c = new capacitor(cap_name, cap_idx, value, BRANCH);
					c->set_nodes( itree_node_list_ref.at(k), _ground_node);
					
					itree_node_list_ref.at(k)->add_capacitor(c);
					_ground_node->add_capacitor(c);
				}
			}
		}
	}
}

void power_grid::add_c4_pads( int &res_idx, int &cap_idx )
{
	// Creating C4s
	profiler record( task_desc );
	int num_c4s = 0, ind_idx = 0;
	c4_pad * c4;
	double c4_distance = _C4_SPACING_FACTOR*get_c4_layer()->get_pitch();
	bool is_horizontal = (get_c4_layer()->get_orientation() == HORIZONTAL);
	double stop_it = is_horizontal?_GRID_DIMENSION_Y_UM:_GRID_DIMENSION_X_UM;

	map<double, interconnect_tree*> c4_itree_map = get_c4_layer()->get_it_map();
	map<double, interconnect_tree*>::iterator itr;

	for (double dist_it = get_c4_layer()->get_offset(); dist_it <= stop_it;
		dist_it += c4_distance)
	{
		itr = c4_itree_map.find( dist_it );
 		if (itr == c4_itree_map.end())		continue;

		interconnect_tree* it = itr->second;
		vector<node*> itree_node_list = it->get_list_of_nodes();
		double off = is_horizontal?
			itree_node_list.at(0)->get_xcoord()
			:itree_node_list.at(0)->get_ycoord();
		
		double stop_node = is_horizontal?
			itree_node_list.at(itree_node_list.size()-1)->get_xcoord()
			:itree_node_list.at(itree_node_list.size()-1)->get_ycoord();
		
		// This code finds the nearest node around the desired position of C4. No new nodes
		// are created in the grid.
		for (double dist_node = off; dist_node < stop_node; dist_node += c4_distance)
		{
			node* c4_node = NULL;
			double min_diff = 1e40;
			if ( is_horizontal )
			{
				for( unsigned int i = 0; i < itree_node_list.size(); i++ )
				{
					if ( abs( itree_node_list[i]->get_xcoord() - dist_node ) < min_diff )
					{
						c4_node = itree_node_list[i];
						min_diff = abs( itree_node_list[i]->get_xcoord() - dist_node );
					}
				}
			}
			else
			{
				for( unsigned int i = 0; i < itree_node_list.size(); i++ )
					if ( abs( itree_node_list[i]->get_ycoord() - dist_node) < min_diff )
					{
						c4_node = itree_node_list[i];
						min_diff = abs( itree_node_list[i]->get_ycoord() - dist_node);
					}
			}
			//cout<<"min_diff: "<<min_diff<<endl;
			assert(c4_node!= NULL);

			num_c4s++;
			stringstream c4ss;
			c4ss << "c4" << num_c4s;
			string c4_name = c4ss.str();
			c4 = new c4_pad(this, c4_node, c4_name, num_c4s);
			c4->parse_subcircuit_file(_number_of_nodes, res_idx, cap_idx, ind_idx);
			c4_node->add_c4_pad_connection(c4);
			_list_of_c4_pads.push_back(c4);
		}//*/
	}
	_number_of_c4s = num_c4s;
}

void power_grid::finalize_node_name_and_indices()
{
	profiler record( task_desc );
	string node_name;
	int node_index = 0;
	/********* Starting with C4s **********/
	for (int i = _list_of_c4_pads.size() - 1; i >=0; i--)
	{
		vector<node*> c4_node_list = _list_of_c4_pads.at(i)->get_list_of_nodes();
		for (int j = c4_node_list.size() - 1; j >= 0; j--)
		{
			node_index++;
			stringstream ss;
			ss << "_n" << _list_of_c4_pads.at(i)->get_name() << "_" << node_index;
			node_name = ss.str();
			c4_node_list.at(j)->set_index(node_index);
			c4_node_list.at(j)->set_name(node_name);
		}
	}

	/********* Grid Nodes **********/
	// temporarily adjust the indices for each layer: required for proper node
	// naming convetion in spice file
	int net_idx = 1;
	for (int i = _number_of_layers-1; i >= 0; i--)
	{
		for(unsigned int j = 0; j < _list_of_layers[i].size(); j++)
		{
			//int net_idx = ( _number_of_layers - i )*10 + j;
			_list_of_layers[i][j]->set_index( net_idx++ );
		}
	}

	double mulf = 1/get_unit();
	for (unsigned int i = 0 ; i < _number_of_layers_global_grid ; i++)
	{
		layer* gl_layer = _list_of_layers.at(i).at(0);
		vector<interconnect_tree*> gl_layer_itree_list = gl_layer->get_list_of_interconnect_trees();
		
		for (int j =  gl_layer_itree_list.size() - 1; j >= 0; j--)
		{
			vector<node*> itree_node_list = gl_layer_itree_list.at(j)->get_list_of_nodes();
			for (int k = itree_node_list.size() - 1; k >= 0; k--)
			{
				node_index++;
				/*stringstream ss;
				ss << "n" << gl_layer->get_index() << "_" <<
					itree_node_list.at(k)->get_xcoord()*mulf << "_" <<
					itree_node_list.at(k)->get_ycoord()*mulf;// << "_" << node_index;
				node_name = ss.str();*/
				char nd_name[128];
				int rl = sprintf( nd_name, "n%d_%d_%d", gl_layer->get_index(),
					(int)(itree_node_list.at(k)->get_xcoord()*mulf),
					(int)(itree_node_list.at(k)->get_ycoord()*mulf) );
				assert( rl >= 0 );
				node_name.assign( nd_name );
				itree_node_list.at(k)->set_index(node_index);
				itree_node_list.at(k)->set_name(node_name);
			}
		}
	}

	for(unsigned int p = 0; p < _number_of_sub_grids ; p++)
	{
		for (unsigned int i = 0; i < _number_of_layers_sub_grid ; i++)
		{
			layer* sg_layer = _list_of_sub_grids.at(p)->get_list_of_layers().at(i).at(0);
			vector<interconnect_tree*> sg_layer_itree_list = sg_layer->get_list_of_interconnect_trees();
			for (int j =  sg_layer_itree_list.size() - 1; j >= 0; j--)
			{
				vector<node*> itree_node_list = sg_layer_itree_list.at(j)->get_list_of_nodes();
				for (int k = itree_node_list.size() - 1; k >= 0; k--)
				{
					node_index++;
					/*stringstream ss;
					ss << "n" << sg_layer->get_index() << "_" << 
						itree_node_list.at(k)->get_xcoord()*mulf << "_" << 
						itree_node_list.at(k)->get_ycoord()*mulf;// << "_" << node_index;
					node_name = ss.str();*/
					char nd_name[128];
					int rl = sprintf( nd_name, "n%d_%d_%d", sg_layer->get_index(),
						(int)(itree_node_list.at(k)->get_xcoord()*mulf),
						(int)(itree_node_list.at(k)->get_ycoord()*mulf) );
					assert( rl >= 0 );
					node_name.assign( nd_name );
					itree_node_list.at(k)->set_index(node_index);
					itree_node_list.at(k)->set_name(node_name);
				}
			}
		}
	}

	// restore the node indices
	for (int i = _number_of_layers-1; i >= 0; i--)
	{
		for(unsigned int j = 0; j < _list_of_layers[i].size(); j++)
		{
			_list_of_layers[i][j]->set_index( _number_of_layers - i );
		}
	}
}

// Getter Functions

//! This method returns a copy of the list of layers in the grid
//! \return vector<layer*>
vector<layer_specs* > power_grid::get_list_of_layer_specs()
{
	return _list_of_layer_specs;
}

//! This method returns a copy of the list of layers in the grid
//! \return vector<layer*>
vector<vector<layer*> > power_grid::get_list_of_layers()
{
	return _list_of_layers;
}

//! This method returns a copy of the list of internal grids in the grid
//! \return vector<layer*>
vector<sub_grid* > power_grid::get_list_of_sub_grids()
{
	return _list_of_sub_grids;
}

//! This method returns a copy of the list of c4 pads in the grid
//! \return vector<c4_pad*>
vector<c4_pad*> power_grid::get_list_of_c4_pads()
{
	return _list_of_c4_pads;
}

//! This method returns a copy of the list of vias in the grid
//! \return vector<resistor*>
vector<resistor*> power_grid::get_list_of_vias()
{
	return _list_of_vias;
}

//! This method returns a copy of the list of current_sources in the grid
//! \return vector<current_source*>
vector<current_source*> power_grid::get_list_of_current_sources()
{
	return _list_of_current_sources;
}

//! This method returns a pointer to the layer of the given index
//! (Convention: Layers with higher index are at the beginning of
//! the list of layers, and the index starts from 1)
//! @param index [in]: index of the layer wanted
//! \return layer*
vector<layer*> power_grid::get_layer(int index) 
{
	if (index < 1 || (unsigned int)index > _number_of_layers)
	{
		cout<<index<<" "<<_number_of_layers<<endl;
		fprintf(stderr, "\n\033[1;%dmIndex of layer is out of bounds. Exiting.\033[0m\n", 31);
		exit(0);
	}
	return _list_of_layers.at(_number_of_layers - index);
}

//! This method returns a pointer to the C4 layer
//! (Convention: C4 layer is the first layer in the list of layers)
//! \return layer*
layer* power_grid::get_c4_layer()
{
	return _list_of_layers.at(0).at(0);
}

//! This method returns a pointer to the current sources layer
//! (Convention: CS layer is the last layer in the list of layers)
//! \return layer*
vector<layer*> power_grid::get_cs_layer()
{
	return _list_of_layers.at(_number_of_layers-1);
}

//! This method returns a pointer to the ground node in the grid
//! \return node*
node* power_grid::get_ground_node()
{
   return _ground_node;
}

//! This method returns a pointer to the Vdd node in the grid
//! \return node*
node* power_grid::get_vdd_node()
{
   return _vdd_node;
}
//! This method returns the number of nodes in the grid
//! \return unsigned int
unsigned int power_grid::get_number_of_nodes()
{
	return _number_of_nodes;
}

//! This method returns the number of layers in the grid
//! \return unsigned int
unsigned int power_grid::get_number_of_layers()
{
	return _number_of_layers;
}

//! This method returns the number of layers in the global grid
//! \return unsigned int
unsigned int power_grid::get_number_of_layers_global_grid()
{
	return _number_of_layers_global_grid;
}

//! This method returns the number of c4 pads in the grid
//! \return unsigned int
unsigned int power_grid::get_number_of_c4s()
{
	return _number_of_c4s;
}

//! This method returns the number of sub grids in the grid in the X Dimension
//! \return unsigned int
unsigned int power_grid::get_number_of_sub_grids_X_DIM()
{
	return _number_of_sub_grids_X_DIM;
}
//! This method returns the number of sub grids in the grid in the Y Dimension
//! \return unsigned int
unsigned int power_grid::get_number_of_sub_grids_Y_DIM()
{
	return _number_of_sub_grids_Y_DIM;
}
//! This method returns the number of sub grids in the grid
//! \return unsigned int
unsigned int power_grid::get_number_of_sub_grids()
{
	return _number_of_sub_grids;
}

//! This method returns the number of branches in the grid
//! \return unsigned int
unsigned int power_grid::get_number_of_branches()
{
	return _number_of_branches;
}

//! This method returns the number of current sources in the grid
//! \return unsigned int
unsigned int power_grid::get_number_of_current_sources()
{
	return _number_of_current_sources;
}



//! This method returns the x dimension of the grid
//! \return double
double power_grid::get_grid_dimension_x()
{
	return _GRID_DIMENSION_X_UM;
}

//! This method returns the y dimension of the grid
//! \return double
double power_grid::get_grid_dimension_y()
{
	return _GRID_DIMENSION_Y_UM;
}

//! This method returns the unit used for all the parameters in the grid
//! \return double
double power_grid::get_unit()
{
	return _unit;
}

//! This method returns the vdd supplied to the grid
//! \return double
double power_grid::get_vdd()
{
	return _vdd;
}

//! This method returns the delta t of the backward euler
//! \return double
double power_grid::get_delta_t()
{
	return _delta_t;
}

//! This method returns the average power density of the grid
//! \return double
double power_grid::get_avg_power_density()
{
	return _AVG_POWER_DESNITY_WATT_PER_CM_SQ;
}

//! This method returns the Total power dissipation of the grid
//! \return double
double power_grid::get_total_power_dissipation_watts()
{
	return _power_scale_factor*_AVG_POWER_DESNITY_WATT_PER_CM_SQ*
		_GRID_DIMENSION_X_UM*_GRID_DIMENSION_Y_UM*1e-8;
}

//! This method returns the ratio of current sources to number of
//! nodes in M1
//! \return double
double power_grid::get_ratio_current_sources_in_M1()
{
	return _RATIO_CS_IN_M1;
}

//! This method returns the C4 spaceing factor in top layer
//! \return double
double power_grid::get_C4_spacing_factor()
{
	return _C4_SPACING_FACTOR;
}

//! This method returns the C4 resistance value 
//! \return double
double power_grid::get_C4_resistance()
{
	return _C4_RESISTANCE;
}

//! This method returns the C4 inductance value (if the grid is RLC)
//! \return double
double power_grid::get_C4_inductance()
{
	return _C4_INDUCTANCE;
}

//! This method returns a reference copy of the G matrix of the grid
//! \return sparse_matrix
matrix::sparse_matrix & power_grid::get_G_matrix_ref()
{
	return _G;
}

//! This method returns a copy of the G matrix of the grid
//! \return sparse_matrix
matrix::sparse_matrix power_grid::get_G_matrix()
{
	return _G;
}

//! This method returns a reference copy of the C matrix of the grid
//! \return sparse_matrix
matrix::sparse_matrix & power_grid::get_C_matrix_ref()
{
	return _C;
}

//! This method returns a copy of the C matrix of the grid
//! \return sparse_matrix
matrix::sparse_matrix power_grid::get_C_matrix()
{
	return _C;
}

//! This method returns a reference copy of the A = G+C/Delta t matrix of the grid
//! \return sparse_matrix
matrix::sparse_matrix & power_grid::get_A_matrix_ref()
{
	return _A;
}

//! This method returns a copy of the A matrix of the grid
//! \return sparse_matrix
matrix::sparse_matrix power_grid::get_A_matrix()
{
	return _G;
}

//! This method returns a reference copy of the L matrix of the grid
//! \return sparse_matrix
matrix::sparse_matrix & power_grid::get_L_matrix_ref()
{
	return _L;
}

//! This method returns a copy of the L matrix of the grid
//! \return sparse_matrix
matrix::sparse_matrix power_grid::get_L_matrix()
{
	return _L;
}

//! This method returns a reference copy of the M matrix of the grid
//! \return sparse_matrix
matrix::sparse_matrix & power_grid::get_M_matrix_ref()
{
	return _M;
}

//! This method returns a copy of the M matrix of the grid
//! \return sparse_matrix
matrix::sparse_matrix power_grid::get_M_matrix()
{
	return _M;
}

// Setter Functions

//! This method sets the number of current sources in the grid 
//! @params num_cs [in]: number of current sources
void power_grid::set_number_of_current_sources(int num_cs)
{
	_number_of_current_sources = num_cs;
}

//! This method adds a layer to the list of layers in the grid 
//! @params l [in]: layer to add
void power_grid::add_layer(layer* l)
{
	_list_of_layers.at(_number_of_layers - l->get_index() + (_lowest_layer_index - 1)).push_back(l);
}

//! This method adds a subgrid to the list of subgrids in the grid 
//! @params sg [in]: subgrid to add
void power_grid::add_sub_grid(sub_grid* sg)
{
	_list_of_sub_grids.at(sg->get_index()-1) = sg;
}

//! This method adds a c4 pad to the list of c4 pads in the grid 
//! @params c4 [in]: c4 pad to add
void power_grid::add_c4_pad(c4_pad* c4)
{
	_list_of_c4_pads.push_back(c4);
}

//! This method adds a via to the list of vias in the grid 
//! @params via [in]: via to add
void power_grid::add_via(resistor* via)
{
	_list_of_vias.push_back(via);
}

//! This method adds a current source to the list of current sources in the grid 
//! @params cs [in]: current source to add
void power_grid::add_current_source(current_source* cs)
{
	_list_of_current_sources.push_back(cs);
}

// Other Functions 

//! This method reads the parameters of the grid to generate from the given parameters file
//! @param params_file [in]: name of the file to read the parameters from
void power_grid::read_config_file(string params_file)
{
	profiler record("Reading Configuration Options File"); 

	using namespace rand_options;

	vector<layer_specs*> layers;
	layers.resize(9);
	layers[0] = new layer_specs();
	layers[1] = new layer_specs();
	layers[2] = new layer_specs();
	layers[3] = new layer_specs();
	layers[4] = new layer_specs();
	layers[5] = new layer_specs();
	layers[6] = new layer_specs();
	layers[7] = new layer_specs();
	layers[8] = new layer_specs();

	ifstream in;
	in.open(params_file.c_str());
	if (!in)
	{
		cerr << "Could not open file for reading configuration options" << endl;
		exit(-1);
	}
	while (!in.eof())
	{
		char str[255] = "";
		in.getline(str, 255);
		if (str[0] != '\0')
		{
			char * param = strtok(str, " =");
			char * str_option = strtok(NULL, " =");
			if (param[0] == '#' || isblank(param[0]))
			  continue;                             // '#' is used for comments
			/********* Input Options **********/
			else if (strcmp(param, "number_of_layers") == 0)
			{
				_number_of_layers = atof(str_option);
				if (_number_of_layers > 9)
				{
					fprintf(stdout, "\033[1;%dmERROR: Cannot Create more than 9 layers."
						" Exiting Now.\033[0m\n", 31);
					exit(1);
				}
				if (_number_of_layers < 3)
				{
					fprintf(stdout, "\033[1;%dmERROR: There should be at least 2 layers"
						" in a grid. Exiting Now.\033[0m\n", 31);
					exit(1);
				}

				_list_of_layers.resize(_number_of_layers);
				for(unsigned int i = 0 ; i < _number_of_layers ; i++)
					_list_of_layers.at(i).resize(0);
				_list_of_layer_specs.resize(_number_of_layers);
			}
			else if (strcmp(param, "number_of_layers_global_grid") == 0)
			{
				_number_of_layers_global_grid = atof(str_option);
				if ( _number_of_layers_global_grid < 2 )
				{
					fprintf(stdout, "\033[1;%dmERROR: Number of layers in Global Grid "
						"has to be greater than 2. Exiting Now.\033[0m\n", 31);
					exit(1);	
				}
				if (_number_of_layers_global_grid > _number_of_layers)
				{
					fprintf(stdout, "\033[1;%dmERROR: Number of layers in Global Grid "
						"exceeds the total number of layers. Exiting Now.\033[0m\n", 31);
					exit(1);
				}
				if ( _number_of_layers - _number_of_layers_global_grid > 4 )
				{
					fprintf(stdout, "\033[1;%dmERROR: Number of layers in sub-grids "
						"cannot be more than 4. Exiting Now.\033[0m\n", 31);
					exit(1);
				}
			}
			else if (strcmp(param, "number_of_sub_grids_X_DIM") == 0)
				_number_of_sub_grids_X_DIM = round( atof(str_option) );
			else if (strcmp(param, "number_of_sub_grids_Y_DIM") == 0)
				_number_of_sub_grids_Y_DIM = round( atof(str_option) );
			/******** Grid Dimensions *********/
			else if (strcmp(param, "GRID_DIMENSION_X_UM") == 0)
				_GRID_DIMENSION_X_UM = atof(str_option);
			else if (strcmp(param, "GRID_DIMENSION_Y_UM") == 0)
				_GRID_DIMENSION_Y_UM = atof(str_option);
			else if (strcmp(param, "UNIT") == 0)
				_unit = atof(str_option);
			/******** Technology Used ********/
			else if (strcmp(param, "VDD") == 0)
				_vdd = atof(str_option);
			/******** DIE Capacitances ***/
			else if (strcmp(param, "MIN_DIE_CAPACITANCE_F") == 0)
				_MIN_DIE_CAPACITANCE_F = atof(str_option);
			else if (strcmp(param, "MAX_DIE_CAPACITANCE_F") == 0)
				_MAX_DIE_CAPACITANCE_F = atof(str_option);
			/******** Layer Orientation *******/
			else if (strcmp(param, "ORIENTATION_M1") == 0)
			{
				if (strcmp(str_option, "HORIZONTAL") == 0)
					 layers.at(0)->orientation = HORIZONTAL;
				else if (strcmp(str_option, "VERTICAL") == 0)
					 layers.at(0)->orientation = VERTICAL;
				else
				{
					cerr << "Unknown orientation of layer" << endl;
					exit(-1);
				}
			}
			/******* Layer Pitch ********/
			else if (strcmp(param, "BASE_PITCH_M1_UM") == 0)
				layers.at(0)->pitch = atof(str_option);
			else if (strcmp(param, "BASE_PITCH_M2_UM") == 0)
				layers.at(1)->pitch = atof(str_option);
			else if (strcmp(param, "BASE_PITCH_M3_UM") == 0)
				layers.at(2)->pitch = atof(str_option);
			else if (strcmp(param, "BASE_PITCH_M4_UM") == 0)
				layers.at(3)->pitch = atof(str_option);
			else if (strcmp(param, "BASE_PITCH_M5_UM") == 0)
				layers.at(4)->pitch = atof(str_option);
			else if (strcmp(param, "BASE_PITCH_M6_UM") == 0)
				layers.at(5)->pitch = atof(str_option);
			else if (strcmp(param, "BASE_PITCH_M7_UM") == 0)
				layers.at(6)->pitch = atof(str_option);
			else if (strcmp(param, "BASE_PITCH_M8_UM") == 0)
				layers.at(7)->pitch = atof(str_option);
			else if (strcmp(param, "BASE_PITCH_M9_UM") == 0)
				layers.at(8)->pitch = atof(str_option);
			/******* Layer Width *****/
			else if (strcmp(param, "WIDTH_M1_UM") == 0)
				layers.at(0)->width = atof(str_option);
			else if (strcmp(param, "WIDTH_M2_UM") == 0)
				layers.at(1)->width = atof(str_option);
			else if (strcmp(param, "WIDTH_M3_UM") == 0)
				layers.at(2)->width = atof(str_option);
			else if (strcmp(param, "WIDTH_M4_UM") == 0)
				layers.at(3)->width = atof(str_option);
			else if (strcmp(param, "WIDTH_M5_UM") == 0)
				layers.at(4)->width = atof(str_option);
			else if (strcmp(param, "WIDTH_M6_UM") == 0)
				layers.at(5)->width = atof(str_option);
			else if (strcmp(param, "WIDTH_M7_UM") == 0)
				layers.at(6)->width = atof(str_option);
			else if (strcmp(param, "WIDTH_M8_UM") == 0)
				layers.at(7)->width = atof(str_option);
			else if (strcmp(param, "WIDTH_M9_UM") == 0)
				layers.at(8)->width = atof(str_option);
			/******* Layer Offset *******/
			else if (strcmp(param, "OFFSET_M1_UM") == 0)
				layers.at(0)->offset = atof(str_option);
			else if (strcmp(param, "OFFSET_M2_UM") == 0)
				layers.at(1)->offset = atof(str_option);
			else if (strcmp(param, "OFFSET_M3_UM") == 0)
				layers.at(2)->offset = atof(str_option);
			else if (strcmp(param, "OFFSET_M4_UM") == 0)
				layers.at(3)->offset = atof(str_option);
			else if (strcmp(param, "OFFSET_M5_UM") == 0)
				layers.at(4)->offset = atof(str_option);
			else if (strcmp(param, "OFFSET_M6_UM") == 0)
				layers.at(5)->offset = atof(str_option);
			else if (strcmp(param, "OFFSET_M7_UM") == 0)
				layers.at(6)->offset = atof(str_option);
			else if (strcmp(param, "OFFSET_M8_UM") == 0)
				layers.at(7)->offset = atof(str_option);
			else if (strcmp(param, "OFFSET_M9_UM") == 0)
				layers.at(8)->offset = atof(str_option);
			/****** Layer Sheet Resistance *****/
			else if (strcmp(param, "SHEET_RESISTANCE_M1_OHM_SQ") == 0)
				layers.at(0)->sheet_resistance = atof(str_option);
			else if (strcmp(param, "SHEET_RESISTANCE_M2_OHM_SQ") == 0)
				layers.at(1)->sheet_resistance = atof(str_option);
			else if (strcmp(param, "SHEET_RESISTANCE_M3_OHM_SQ") == 0)
				layers.at(2)->sheet_resistance = atof(str_option);
			else if (strcmp(param, "SHEET_RESISTANCE_M4_OHM_SQ") == 0)
				layers.at(3)->sheet_resistance = atof(str_option);
			else if (strcmp(param, "SHEET_RESISTANCE_M5_OHM_SQ") == 0)
				layers.at(4)->sheet_resistance = atof(str_option);
			else if (strcmp(param, "SHEET_RESISTANCE_M6_OHM_SQ") == 0)
				layers.at(5)->sheet_resistance = atof(str_option);
			else if (strcmp(param, "SHEET_RESISTANCE_M7_OHM_SQ") == 0)
				layers.at(6)->sheet_resistance = atof(str_option);
			else if (strcmp(param, "SHEET_RESISTANCE_M8_OHM_SQ") == 0)
				layers.at(7)->sheet_resistance = atof(str_option);
			else if (strcmp(param, "SHEET_RESISTANCE_M9_OHM_SQ") == 0)
				layers.at(8)->sheet_resistance = atof(str_option);
			/***** Layer Via Resistance ******/
			else if (strcmp(param, "VIA_RESISTANCE_M2_OHM") == 0)
				layers.at(1)->via_resistance = atof(str_option);
			else if (strcmp(param, "VIA_RESISTANCE_M3_OHM") == 0)
				layers.at(2)->via_resistance = atof(str_option);
			else if (strcmp(param, "VIA_RESISTANCE_M4_OHM") == 0)
				layers.at(3)->via_resistance = atof(str_option);
			else if (strcmp(param, "VIA_RESISTANCE_M5_OHM") == 0)
				layers.at(4)->via_resistance = atof(str_option);
			else if (strcmp(param, "VIA_RESISTANCE_M6_OHM") == 0)
				layers.at(5)->via_resistance = atof(str_option);
			else if (strcmp(param, "VIA_RESISTANCE_M7_OHM") == 0)
				layers.at(6)->via_resistance = atof(str_option);
			else if (strcmp(param, "VIA_RESISTANCE_M8_OHM") == 0)
				layers.at(7)->via_resistance = atof(str_option);
			else if (strcmp(param, "VIA_RESISTANCE_M9_OHM") == 0)
				layers.at(8)->via_resistance = atof(str_option);
			// C4 cpacing
			else if (strcmp(param, "C4_SPACING_FACTOR") == 0)
				_C4_SPACING_FACTOR = atof(str_option);
            // C4 resistance and C4 inductance
            else if (strcmp(param, "C4_RESISTANCE") == 0)
                _C4_RESISTANCE = atof(str_option);
            else if (strcmp(param, "C4_INDUCTANCE") == 0)
                _C4_INDUCTANCE = atof(str_option);
			// Power density
			else if (strcmp(param, "AVG_POWER_DESNITY_WATT_PER_CM_SQ") == 0)
				_AVG_POWER_DESNITY_WATT_PER_CM_SQ = atof(str_option);
			// ratio of current sources in M1
			else if (strcmp(param, "RATIO_CS_IN_M1") == 0)
			{
				_RATIO_CS_IN_M1 = atof(str_option);
				if ( _RATIO_CS_IN_M1 > 1 )	_RATIO_CS_IN_M1 = 1;
			}
			else if (strcmp(param, "randomize_sub_grids") == 0)
				_randomize_sub_grids = atoi(str_option);
			else if (strcmp(param, "ratio_of_ON_sub_grids") == 0)
			{
				_ratio_of_ON_sub_grids = atof(str_option);
				if ( _ratio_of_ON_sub_grids > 1 )	_ratio_of_ON_sub_grids = 1;
			}
			else if (strcmp(param, "ratio_of_retained_trees_in_sub_grids") == 0)
			{
				_ratio_of_retained_trees_in_sub_grids = atof(str_option);
				if ( _ratio_of_retained_trees_in_sub_grids > 1 )
					_ratio_of_retained_trees_in_sub_grids = 1;
			}
			else if (strcmp(param, "ratio_of_retained_trees_in_intermediate_layers") == 0)
			{
				_ratio_of_retained_trees_in_intermediate_layers = atof(str_option);
				if ( _ratio_of_retained_trees_in_intermediate_layers > 1 )
					_ratio_of_retained_trees_in_intermediate_layers = 1;
			}
			else if (strcmp(param, "ratio_of_retained_trees_in_global_layers") == 0)
			{
				_ratio_of_retained_trees_in_global_layers = atof(str_option);
				if ( _ratio_of_retained_trees_in_global_layers > 1 )
					_ratio_of_retained_trees_in_global_layers = 1;
			}
			else if (strcmp(param, "max_shift_tree_spans") == 0)
			{
				_max_shift_tree_spans = atof(str_option);
				if ( _max_shift_tree_spans > 0.4 )	_max_shift_tree_spans = 0.4;
			}
			else if (strcmp(param, "seed_offset") == 0)
			{
				_seed_offset = atof(str_option);
			}
		}
	}
	in.close();
	if (layers.at(0)->orientation == HORIZONTAL)
	{
		layers.at(1)->orientation = VERTICAL;
		layers.at(2)->orientation = HORIZONTAL;
		layers.at(3)->orientation = VERTICAL;
		layers.at(4)->orientation = HORIZONTAL;
		layers.at(5)->orientation = VERTICAL;
		layers.at(6)->orientation = HORIZONTAL;
		layers.at(7)->orientation = VERTICAL;
		layers.at(8)->orientation = HORIZONTAL;
	}
	else
	{
		layers.at(1)->orientation = HORIZONTAL;
		layers.at(2)->orientation = VERTICAL;
		layers.at(3)->orientation = HORIZONTAL;
		layers.at(4)->orientation = VERTICAL;
		layers.at(5)->orientation = HORIZONTAL;
		layers.at(6)->orientation = VERTICAL;
		layers.at(7)->orientation = HORIZONTAL;
		layers.at(8)->orientation = VERTICAL;
	}

	/*************/
	layers.at(0)->via_resistance = -1;
	/*************/

	for (unsigned int i = 0; i < _number_of_layers; i++)
	{
		_list_of_layer_specs.at(_number_of_layers - i - 1) = layers.at(i);
	}
	// free memory from layer_spec which is ot needed in future. For example, if
	// user requested only 5 layers, then layer_spec 5, 6, 7 and 8 need to be
	// freed
	for ( unsigned int i = _number_of_layers; i < layers.size(); i++ )
		delete layers[i];

	// error message template
	int to_exit = 0;
	char err_msg[] = "\033[1;31mNOTE: Input '%s' should be greater than 0.\033[0m\n";

	// finally, validate values
	if ( _number_of_sub_grids_X_DIM <= 0 )
		{ printf(err_msg, "number_of_sub_grids_X_DIM" ); to_exit++; }
	if ( _number_of_sub_grids_Y_DIM <= 0 )
		{ printf(err_msg, "number_of_sub_grids_Y_DIM" ); to_exit++; }
	if ( _GRID_DIMENSION_X_UM <= 0 )
		{ printf(err_msg, "GRID_DIMENSION_X_UM" ); to_exit++; }
	if ( _GRID_DIMENSION_Y_UM <= 0 )
		{ printf(err_msg, "GRID_DIMENSION_Y_UM" ); to_exit++; }
	if ( _unit <= 0 )
		{ printf(err_msg, "unit" ); to_exit++; }
	if ( _vdd <= 0 )
		{ printf(err_msg, "vdd" ); to_exit++; } 
	if ( _MIN_DIE_CAPACITANCE_F <= 0 )
		{ printf(err_msg, "MIN_DIE_CAPACITANCE_F" ); to_exit++; }
	if ( _MAX_DIE_CAPACITANCE_F <= 0 )
		{ printf(err_msg, "MIN_DIE_CAPACITANCE_F" ); to_exit++; }
	if ( _C4_SPACING_FACTOR <= 0 )
		{ printf(err_msg, "C4_SPACING_FACTOR" ); to_exit++; }
	if ( _AVG_POWER_DESNITY_WATT_PER_CM_SQ <= 0 )
		{ printf(err_msg, "AVG_POWER_DESNITY_WATT_PER_CM_SQ" ); to_exit++; }
	if ( _ratio_of_ON_sub_grids <= 0 )
		{ printf(err_msg, "ratio_of_ON_sub_grids" ); to_exit++; }
	if ( _ratio_of_retained_trees_in_sub_grids  <= 0 )
		{ printf(err_msg, "ratio_of_retained_trees_in_sub_grids " ); to_exit++; }
	if ( _ratio_of_retained_trees_in_intermediate_layers <= 0 )
		{ printf(err_msg, "ratio_of_retained_trees_in_intermediate_layers" ); to_exit++; }
	if ( _ratio_of_retained_trees_in_global_layers <= 0 )
		{ printf(err_msg, "ratio_of_retained_trees_in_global_layers" ); to_exit++; }
	//if ( _max_shift_tree_spans <= 0 )
	//	{ printf(err_msg, "max_shift_tree_spans" ); to_exit++; }

	for( unsigned int i = 0; i < _list_of_layer_specs.size(); i++ )
	{
		if ( !_list_of_layer_specs[i]->validate_values() )
		{
			fprintf(stderr, "\033[1;%dmValidation of layer specifications for layer M%d"
				" failed.\033[0m\n", 31, _number_of_layers-i );
			to_exit++;
		}
	}

	if (to_exit)
	{
		printf("\n******** EXITING!!! %d invalid inputs. ********\n", to_exit);
		exit(1);
	}
}

//! This method builds the G matrix of the grid 
void power_grid::make_G_matrix()
{
	//fprintf(stdout, "\033[1;%dmBuilding G Matrix\033[0m\n",36);
	profiler record("Building G matrix");
	_G.set_number_of_rows(_number_of_nodes);
	_G.set_number_of_columns(_number_of_nodes);

	int num_res = 0;
	// Going over all branches in the grid (excluding C4s)
	vector<resistor*> branches;
	for (int i = _number_of_layers_global_grid-1; i >= 0; i--)
	{
		branches = _list_of_layers.at(i).at(0)->get_list_of_resistors();
		num_res += branches.size();
		for (unsigned int j = 0; j < branches.size(); j++)
		{
			_G.get_Ti().push_back(branches[j]->get_n0()->get_index()-1);
			_G.get_Tj().push_back(branches[j]->get_n0()->get_index()-1);
			_G.get_Tx().push_back(1/branches[j]->get_value());

			_G.get_Ti().push_back(branches[j]->get_n1()->get_index()-1);
			_G.get_Tj().push_back(branches[j]->get_n1()->get_index()-1);
			_G.get_Tx().push_back(1/branches[j]->get_value());

			_G.get_Ti().push_back(branches[j]->get_n0()->get_index()-1);
			_G.get_Tj().push_back(branches[j]->get_n1()->get_index()-1);
			_G.get_Tx().push_back(-1/branches[j]->get_value());

			_G.get_Ti().push_back(branches[j]->get_n1()->get_index()-1);
			_G.get_Tj().push_back(branches[j]->get_n0()->get_index()-1);
			_G.get_Tx().push_back(-1/branches[j]->get_value());
		}
	}
	branches.clear();

	for(unsigned int p = 0; p < _number_of_sub_grids ; p++)
	{
		for (int i =_list_of_sub_grids.at(p)->get_number_of_layers()-1; i >= 0; i--)
		{
			branches = _list_of_sub_grids.at(p)->get_list_of_layers().at(i).at(0)->get_list_of_resistors();
			num_res += branches.size();
			for (unsigned int j = 0; j < branches.size(); j++)
			{
				_G.get_Ti().push_back(branches[j]->get_n0()->get_index()-1);
				_G.get_Tj().push_back(branches[j]->get_n0()->get_index()-1);
				_G.get_Tx().push_back(1/branches[j]->get_value());

				_G.get_Ti().push_back(branches[j]->get_n1()->get_index()-1);
				_G.get_Tj().push_back(branches[j]->get_n1()->get_index()-1);
				_G.get_Tx().push_back(1/branches[j]->get_value());

				_G.get_Ti().push_back(branches[j]->get_n0()->get_index()-1);
				_G.get_Tj().push_back(branches[j]->get_n1()->get_index()-1);
				_G.get_Tx().push_back(-1/branches[j]->get_value());

				_G.get_Ti().push_back(branches[j]->get_n1()->get_index()-1);
				_G.get_Tj().push_back(branches[j]->get_n0()->get_index()-1);
				_G.get_Tx().push_back(-1/branches[j]->get_value());
			}
		}
	}

	// Going over all the vias
	num_res += _list_of_vias.size();
	for (unsigned int i = 0; i < _list_of_vias.size(); i++)
	{
			_G.get_Ti().push_back(_list_of_vias[i]->get_n0()->get_index()-1);
			_G.get_Tj().push_back(_list_of_vias[i]->get_n0()->get_index()-1);
			_G.get_Tx().push_back(1/_list_of_vias[i]->get_value());
					   
			_G.get_Ti().push_back(_list_of_vias[i]->get_n1()->get_index()-1);
			_G.get_Tj().push_back(_list_of_vias[i]->get_n1()->get_index()-1);
			_G.get_Tx().push_back(1/_list_of_vias[i]->get_value());
					   
			_G.get_Ti().push_back(_list_of_vias[i]->get_n0()->get_index()-1);
			_G.get_Tj().push_back(_list_of_vias[i]->get_n1()->get_index()-1);
			_G.get_Tx().push_back(-1/_list_of_vias[i]->get_value());
					   
			_G.get_Ti().push_back(_list_of_vias[i]->get_n1()->get_index()-1);
			_G.get_Tj().push_back(_list_of_vias[i]->get_n0()->get_index()-1);
			_G.get_Tx().push_back(-1/_list_of_vias[i]->get_value());
	}

	// Going over all the C4s
	for (unsigned int i = 0; i < _list_of_c4_pads.size(); i++)
	{
		branches = _list_of_c4_pads[i]->get_list_of_resistors();
		for (unsigned int j = 0; j < branches.size(); j++)
		{
			if (branches[j]->get_n0()->get_index() != 0 && 
				branches[j]->get_n0()->get_index() != -1 &&
				branches[j]->get_n1()->get_index() != 0 &&
				branches[j]->get_n1()->get_index() != -1)
			{
				num_res++;
				_G.get_Ti().push_back(branches[j]->get_n0()->get_index()-1);
				_G.get_Tj().push_back(branches[j]->get_n0()->get_index()-1);
				_G.get_Tx().push_back(1/branches[j]->get_value());
						   
				_G.get_Ti().push_back(branches[j]->get_n1()->get_index()-1);
				_G.get_Tj().push_back(branches[j]->get_n1()->get_index()-1);
				_G.get_Tx().push_back(1/branches[j]->get_value());
						   
				_G.get_Ti().push_back(branches[j]->get_n0()->get_index()-1);
				_G.get_Tj().push_back(branches[j]->get_n1()->get_index()-1);
				_G.get_Tx().push_back(-1/branches[j]->get_value());
						   
				_G.get_Ti().push_back(branches[j]->get_n1()->get_index()-1);
				_G.get_Tj().push_back(branches[j]->get_n0()->get_index()-1);
				_G.get_Tx().push_back(-1/branches[j]->get_value());
			}
			else
			{
				_G.get_Ti().push_back(branches[j]->get_n1()->get_index()-1);
				_G.get_Tj().push_back(branches[j]->get_n1()->get_index()-1);
				_G.get_Tx().push_back(1/branches[j]->get_value());
			}
		}
	}
	_G.set_nz(_number_of_nodes + 2*num_res);
	_G.get_column_format();
	branches.clear();
}

//! This method builds the C matrix of the grid 
void power_grid::make_C_matrix()
{
	fprintf(stdout, "\033[1;%dmBuilding C Matrix\033[0m\n",36);
	profiler record("Building C matrix");
	_C.set_number_of_rows(_number_of_nodes);
	_C.set_number_of_columns(_number_of_nodes);
	_C.set_nz(_number_of_nodes);

	// Going over all capacitors in the grid (excluding C4s)
	vector<node*> nodes;
	for (unsigned int i = 0; i < _number_of_layers_global_grid; i++)
	{
		for (unsigned int j = 0; j < _list_of_layers.at(i).at(0)->get_list_of_interconnect_trees().size(); j++)
		{
			nodes  = _list_of_layers.at(i).at(0)->get_list_of_interconnect_trees().at(j)->get_list_of_nodes_ref();
			for (unsigned int k = 0; k < nodes.size(); k++)
			{
				for (unsigned int m = 0; m < nodes.at(k)->get_c_connections().size(); m++)
				{
					_C.get_Ti().push_back(nodes.at(k)->get_index() - 1);
					_C.get_Tj().push_back(nodes.at(k)->get_index() - 1);
					_C.get_Tx().push_back(nodes.at(k)->get_c_connections().at(m)->get_value());
				}
			}
		}
	}
	//if(_number_of_sub_grids>=2)
	//{
		for(unsigned int p = 0 ; p < _number_of_sub_grids ; p++)
		{
			for (unsigned int i = 0; i < _list_of_sub_grids.at(p)->get_number_of_layers(); i++)
			{
				for (unsigned int j = 0; j < _list_of_sub_grids.at(p)->get_list_of_layers().at(i).at(0)->get_list_of_interconnect_trees().size(); j++)
				{
					nodes  = _list_of_sub_grids.at(p)->get_list_of_layers().at(i).at(0)->get_list_of_interconnect_trees().at(j)->get_list_of_nodes_ref();
					for (unsigned int k = 0; k < nodes.size(); k++)
					{
						for (unsigned int m = 0; m < nodes.at(k)->get_c_connections().size(); m++)
						{
							_C.get_Ti().push_back(nodes.at(k)->get_index() - 1);
							_C.get_Tj().push_back(nodes.at(k)->get_index() - 1);
							_C.get_Tx().push_back(nodes.at(k)->get_c_connections().at(m)->get_value());
						}
					}
				}
			}
		}
	//}
	// Going over C4s
	vector<capacitor*> caps;
	for (unsigned int i = 0; i < _list_of_c4_pads.size(); i++)
	{
		caps = _list_of_c4_pads[i]->get_list_of_capacitors();
		if (caps.size() != 0)
		{
			for (unsigned int j = 0; j < caps.size(); j++)
			{
				_C.get_Ti().push_back(caps[j]->get_n1()->get_index()-1);
				_C.get_Tj().push_back(caps[j]->get_n1()->get_index()-1);
				_C.get_Tx().push_back(caps[j]->get_value());
			}
		}
	}
	_C.get_column_format();
}

//! This method builds the A matrix of the grid 
void power_grid::make_A_matrix()
{
	double dt = get_delta_t();
	sparse_matrix_add(1, get_G_matrix_ref(), 1/dt, get_C_matrix_ref(), _A);
	_A.get_triplet_format();
}

//! This method prints some grid information
void power_grid::print_power_grid_info()
{
	cout << endl << endl << "POWER GRID INFO" << endl << "------------------" <<endl;
	cout << setw(40) << "Grid dimension X: "
		 << get_grid_dimension_x() << " um" << endl;
	cout << setw(40) << "Grid dimension Y: "
		 << get_grid_dimension_y() << " um" << endl;
	cout << setw(40) << "Grid Area: "
		 <<  get_grid_dimension_x()*get_grid_dimension_y() << " um^2" <<endl;
	cout << setw(40) << "Number of layers in the grid: "
		 << _number_of_layers << endl;
	cout << setw(40) << "Number of layers in the global grid: "
		 << _number_of_layers_global_grid << endl;
	cout << setw(40) << "Number of layers in the internal grid: " 
		 << _number_of_layers - _number_of_layers_global_grid << endl;
	cout << setw(40) << "Number of sub_grids in internal grid(s): "
		 << _number_of_sub_grids << endl;
	cout << setw(40) << "Total number of nodes: "
		 << _number_of_nodes << endl;
	cout << setw(40) << "Total number of branches: "
		 << _number_of_branches << endl;
	cout << setw(40) << "Total number of vias: "
		 <<_list_of_vias.size() << endl;
	cout << setw(40) << "Number of c4s in the grid: "
		 << _number_of_c4s << endl;
	cout << setw(40) << "Number of current sources in the grid: "
		 << _number_of_current_sources << endl;
	cout << endl;

	vector<vector<layer*> > layers = _list_of_layers;
	cout << "\nGlobal Grid Layers" << endl;
	cout << "Layer Name" << setw(10) << "Index" << setw(10) << "Width" << setw(10) << "Pitch" << setw(10) << "Offset" << setw(20);
	cout << "Sheet Resistance" << setw(20) << "Via Resistance" << setw(10) << "C4?" << setw(10) << "CS?" << setw(15) << "Orientation";
	cout << setw(10) << "#ITs" << endl;
	string orientation;
	for (unsigned int i = 0; i < _number_of_layers_global_grid; i++)
	{
		orientation = layers[i].at(0)->get_orientation()==0?"H":"V"; 
		cout << layers[i].at(0)->get_name() << setw(15) << layers[i].at(0)->get_index() << setw(13) << layers[i].at(0)->get_width() << setw(10);
		cout << layers[i].at(0)->get_pitch() << setw(10) << layers[i].at(0)->get_offset() << setw(20) << layers[i].at(0)->get_sheet_resistance();
		cout << setw(20) << layers[i].at(0)->get_via_resistance() << setw(10) << layers[i].at(0)->is_c4_layer() << setw(10);
		cout << layers[i].at(0)->is_current_source_layer() << setw(10) << orientation<< setw(15);
		cout << layers[i].at(0)->get_list_of_interconnect_trees().size() << endl;
	}
	for(unsigned int p = 0; p < _number_of_sub_grids; p++)
	{
		vector<vector<layer*> > layers = _list_of_sub_grids[p]->get_list_of_layers();
		cout << "\nInternal Grid #" << _list_of_sub_grids[p]->get_index() << endl;
		cout <<  "Layer Name" << setw(10) << "Index" << setw(10) << "Width" << setw(10) << "Pitch" << setw(10) << "Offset" << setw(20);
		cout << "Sheet Resistance" << setw(20) << "Via Resistance" << setw(10) << "C4?" << setw(10) << "CS?" << setw(15) << "Orientation";
		cout << setw(10) << "#ITs" << endl;
		for (unsigned int i = 0; i < _list_of_sub_grids[p]->get_number_of_layers(); i++)
		{
			orientation = layers[i].at(0)->get_orientation()==0?"H":"V"; 
			cout << layers[i].at(0)->get_name() << setw(15) << layers[i].at(0)->get_index() << setw(13) << layers[i].at(0)->get_width() << setw(10);
			cout << layers[i].at(0)->get_pitch() << setw(10) << layers[i].at(0)->get_offset() << setw(20) << layers[i].at(0)->get_sheet_resistance();
			cout << setw(20) << layers[i].at(0)->get_via_resistance() << setw(10) << layers[i].at(0)->is_c4_layer() << setw(10);
			cout << layers[i].at(0)->is_current_source_layer() << setw(10) << orientation << setw(15);
			cout << layers[i].at(0)->get_list_of_interconnect_trees().size() << endl;
		}
	}
	cout << endl;
	
	for(unsigned int j = 0 ; j < _number_of_layers_sub_grid ; j++)
	{
		int m_nodes = 0, m_branches = 0;
		for(unsigned int p = 0 ; p < _number_of_sub_grids ; p++)
		{
			layer* m = _list_of_sub_grids.at(p)->get_layer(j+1).at(0);
			for (unsigned int i = 0; i < m->get_list_of_interconnect_trees().size(); i++)
			{
				m_nodes += m->get_list_of_interconnect_trees()[i]->
					get_list_of_nodes_ref().size();
				m_branches += m->get_list_of_interconnect_trees()[i]->
					get_list_of_resistors_ref().size();
			}
		}
		cout << "Layer M" << j+1 << ", #nodes: " << setw(6) << m_nodes << ", #branches: "
			 << setw(6) << m_branches << endl;
	}

	for(unsigned int j = 0 ; j < _number_of_layers_global_grid ; j++)
	{ 
		layer* m_last = get_layer(_number_of_layers_sub_grid + j + 1).at(0);
		int m_last_nodes = 0, m_last_branches = 0;
		for (unsigned int i = 0; i < m_last->get_list_of_interconnect_trees().size(); i++)
		{
			m_last_nodes += m_last->get_list_of_interconnect_trees()[i]->
				get_list_of_nodes_ref().size();
			m_last_branches += m_last->get_list_of_interconnect_trees()[i]->
				get_list_of_resistors_ref().size();
		}
		cout << "Layer M" << _number_of_layers_sub_grid + j + 1 << ", #nodes: "<< setw(6)
			 << m_last_nodes << ", #branches "<< setw(6) << m_last_branches << endl;
	}
	
	/*vector<resistor*> links;
	double avg_length = 0;
	for (unsigned int i = 0; i < _number_of_layers_global_grid; i++)
	{
		links = _list_of_layers.at(i).at(0)->get_list_of_resistors();
		for (unsigned j = 0; j < links.size(); j++)
		{
			if (links[j]->get_type() == BRANCH)
				avg_length += links[j]->get_length();
		}
	}
	for(unsigned int p = 0 ; p < _number_of_sub_grids ; p++)
	{
		for (unsigned int i = 0; i < _list_of_sub_grids.at(p)->get_number_of_layers(); i++)
		{
			links = _list_of_sub_grids.at(p)->get_list_of_layers().at(i).at(0)->get_list_of_resistors();
			for (unsigned j = 0; j < links.size(); j++)
			{
				if (links[j]->get_type() == BRANCH)
					avg_length += links[j]->get_length();
			}
		}
	}
	avg_length = avg_length/_number_of_branches;
	cout << "Average length = " << avg_length << endl << endl;*/

	// code to count number of very short branhces, so that ratio of lengths of
	// neighboring branches do not exceed a certain thershold
	vector<resistor*> branches;
	double th1 = 100, th2 = 0.01;
	int num_short_br;
	double maxratio, minratio;
	double maxlen, minlen;
	for (int i = _number_of_layers - 1; i >= 0; i--)
	{
		num_short_br = 0;
		maxratio = -1; minratio = 1e12;
		for(unsigned int p = 0; p < _list_of_layers[i].size(); p++)
		{
			branches = _list_of_layers[i][p]->get_list_of_resistors();
			maxlen   = branches[0]->get_length();
			minlen   = branches[0]->get_length();
			// check
			for (unsigned int j = 1; j < branches.size(); j++)
			{
				double len0 = branches[j-1]->get_length();//*_unit;
				double len1 = branches[j]  ->get_length();//*_unit;
				double ratio = len1/len0;

				if ( ratio < th2 || ratio > th1 )
					num_short_br++;

				if ( ratio < minratio )	minratio = ratio;
				if ( ratio > maxratio )	maxratio = ratio;

				if ( branches[j]->get_length() < minlen )	minlen = branches[j]->get_length();
				if ( branches[j]->get_length() > maxlen )	maxlen = branches[j]->get_length();
			}

		}

		// print status
		if ( num_short_br != 0 )
		{
			fprintf(stdout, "\033[;%dm\n   WARNING: There are some very short branches that may"
				" lead to badly scaled system matrices for EKM model.\033[0m",93);
			fprintf(stdout, "\033[;%dm\n   **** Layer %d (%s)****\n     - Num short branhces: %d,"
				"\n     - min. ratio: %g, max. ratio: %g,\n     - min. len: %g um, max. len: %g um"
				"\033[0m\n", 93,i+1,
				(char*)_list_of_layers[i][0]->get_name().c_str(), num_short_br, minratio,
				maxratio, minlen, maxlen);	
		}
	}
}

//! This method generates a spice file of the power grid
void power_grid::generate_spice_file_for_resistive_grid( string filename )
{
	FILE* spice = fopen(filename.c_str(), "w");
	
	// Number of nodes -- not required  
	// fprintf(spice, "* Number of Nodes: %d\n", _number_of_nodes);
	// Printing Layer Names
	fprintf(spice, "* layers ");
	for (int i = _number_of_layers-1; i >= 0 ; i--)
	{
		stringstream strs;
		strs<<_list_of_layers[i][0]->get_name();
		fprintf(spice, "%s ",  (char*)strs.str().c_str());
	}
	fprintf(spice, "\n");

	fprintf(spice, "* vdd %g\n", _vdd);

	// we forcibly change the index of each layer to something unique for printing
	// at the end, the actual indices are restored
	fprintf(spice, "* VDD nets ");
	int net_idx = 1;
	for (int i = _number_of_layers-1; i >= 0; i--)
	{
		//fprintf(stdout, "\n%d. %s: ", i, (char*)_list_of_layers[i][0]->get_name().c_str() );
		for(unsigned int j = 0; j < _list_of_layers[i].size(); j++)
		{
			_list_of_layers[i][j]->set_index( net_idx++ );
			//fprintf(stdout, "%d:%d, ", j, _list_of_layers[i][j]->get_index() );
			fprintf(spice, "%d ", _list_of_layers[i][j]->get_index() );
		}
	}
	fprintf(spice, "\n");

	fprintf(spice, "* GND nets ");
	fprintf(spice, "\n");

	fprintf(spice, "* via unshorted\n");

	// Printing Vias
	vector<bool> found(_list_of_vias.size(), 0);
	for (int i = _number_of_layers-1; i >= 1; i--)
	{
		for(unsigned int j = 0; j < _list_of_layers[i].size(); j++)
		{
			int via_from_net_idx = _list_of_layers[i][j]->get_index();
			for(unsigned int k = 0; k < _list_of_layers[i-1].size(); k++)
			{
				int via_to_net_idx = _list_of_layers[i-1][k]->get_index();
				bool vias_found = 0;
				for (unsigned int l = 0; l < _list_of_vias.size(); l++)
				{
					if (!found[l] && 
					   ( (_list_of_vias[l]->get_n0()->get_interconnect_tree_ptr()->
							get_layer_ptr()->get_index() == via_from_net_idx 
						 && _list_of_vias[l]->get_n1()->get_interconnect_tree_ptr()->
						 	get_layer_ptr()->get_index() == via_to_net_idx )
						 || (_list_of_vias[l]->get_n1()->get_interconnect_tree_ptr()->
						 	get_layer_ptr()->get_index() == via_from_net_idx 
						 && _list_of_vias[l]->get_n0()->get_interconnect_tree_ptr()->
						 	get_layer_ptr()->get_index() == via_to_net_idx ) ) )
					{
						if ( !vias_found )
						{
							vias_found = 1;
							fprintf(spice, "* vias from: %d to %d\n", via_from_net_idx, via_to_net_idx);
						}

						fprintf(spice, "%s %s %s %g\n", _list_of_vias[l]->get_name().c_str(),
							_list_of_vias[l]->get_n0()->get_name().c_str(),
							_list_of_vias[l]->get_n1()->get_name().c_str(),
							_list_of_vias[l]->get_value());
						found[l] = true;
					}
				}
				if (vias_found)
					fprintf(spice, "* done with vias from: %d to %d\n",
						via_from_net_idx, via_to_net_idx);
			}
		}
	}

	//sanity check
	for( unsigned int i = 0; i < _list_of_vias.size(); i++)
	{
		if ( !found[i] )
		{
			fprintf(stdout, "%s %s %s %g\n", _list_of_vias[i]->get_name().c_str(),
				_list_of_vias[i]->get_n0()->get_name().c_str(),
				_list_of_vias[i]->get_n1()->get_name().c_str(),
				_list_of_vias[i]->get_value());
		}
	}
	
	// Printing C4s
	vector<resistor*>  c4_res;
	//vector<capacitor*> c4_cap;
	//vector<inductor*>  c4_ind;
	for (unsigned int i = 0; i < _number_of_c4s; i++)
	{
		c4_res = _list_of_c4_pads.at(i)->get_list_of_resistors();
		/*// original code follows
		c4_cap = _list_of_c4_pads.at(i)->get_list_of_capacitors();
		c4_ind = _list_of_c4_pads.at(i)->get_list_of_inductors();
		for (unsigned int j = 0; j < c4_res.size(); j++)
		{
			fprintf(spice, "%s %s %s %g\n", c4_res[j]->get_name().c_str(),
					c4_res[j]->get_n0()->get_name().c_str(), c4_res[j]->get_n1()->get_name().c_str(),
					c4_res[j]->get_value());
		}
		for (unsigned int j = 0; j < c4_cap.size(); j++)
		{
			fprintf(spice, "%s %s %s %g\n", c4_cap[j]->get_name().c_str(),
					c4_cap[j]->get_n1()->get_name().c_str(), c4_cap[j]->get_n0()->get_name().c_str(),
					c4_cap[j]->get_value());
		}
		for (unsigned int j = 0; j < c4_ind.size(); j++)
		{
			fprintf(spice, "%s %s %s %g\n", c4_ind[j]->get_name().c_str(),
					c4_ind[j]->get_n0()->get_name().c_str(), c4_ind[j]->get_n1()->get_name().c_str(),
					c4_ind[j]->get_value());
		}*/
		// for printing to spice file, we reverse the node order of C4 resistor
		// and ignore any capactir or inductor in subscircuit
		for (unsigned int j = 0; j < c4_res.size(); j++)
		{
			fprintf(spice, "%s %s %s %g\n", c4_res[j]->get_name().c_str(),
				c4_res[j]->get_n1()->get_name().c_str(), c4_res[j]->get_n0()->get_name().c_str(), 
				c4_res[j]->get_value());
		}

		stringstream ss;
		ss << "vc" << _list_of_c4_pads.at(i)->get_index();
		string temp = ss.str();
		fprintf(spice, "%s %s %s %g\n", temp.c_str(), _vdd_node->get_name().c_str(),
				_ground_node->get_name().c_str(), _vdd);
	}

	// Printing Branch Resistors
	vector<resistor*> branches;
	//vector<node*> nodes;
	//vector<capacitor*> caps;

	for (int i = _number_of_layers - 1; i >= 0; i--)
	{
		for(unsigned int p = 0; p < _list_of_layers[i].size(); p++)
		{
			fprintf(spice, "* layer: %s, VDD net: %d\n",
				(char*)_list_of_layers[i][p]->get_name().c_str(),
				_list_of_layers[i][p]->get_index() ); 
			branches = _list_of_layers[i][p]->get_list_of_resistors();
			for (unsigned int j = 0; j < branches.size(); j++)
			{
				fprintf(spice, "%s %s %s %g\n", branches[j]->get_name().c_str(),
					branches[j]->get_n0()->get_name().c_str(), 
					branches[j]->get_n1()->get_name().c_str(),
					branches[j]->get_value());
			}

			/*nodes = _list_of_layers[i][p]->get_list_of_nodes();
			for (unsigned int j = 0; j < nodes.size(); j++)
			{
				caps = nodes[j]->get_c_connections();
				for (unsigned int k = 0; k < caps.size(); k++)
				{
					fprintf(spice, "%s %s %s %g\n", caps[k]->get_name().c_str(),
						caps[k]->get_n1()->get_name().c_str(), 
						caps[k]->get_n0()->get_name().c_str(),
						caps[k]->get_value());
				}
			}*/
		}
	}

	// signifies start of current sources
	fprintf(spice, "*\n");

	// Printing current sources
	for (unsigned int i = 0; i < _number_of_current_sources; i++)
	{
		fprintf(spice, "%s %s n0 %g\n", _list_of_current_sources[i]->get_name().c_str(),
			_list_of_current_sources[i]->get_n0()->get_name().c_str(),
			_list_of_current_sources[i]->get_value() );
	}
	fprintf(spice, ".op\n.end");
	fclose(spice);

	// finally resore the indices for each layer
	for (int i = _number_of_layers-1; i >= 0; i--)
	{
		for(unsigned int j = 0; j < _list_of_layers[i].size(); j++)
		{
			_list_of_layers[i][j]->set_index( _number_of_layers - i );
		}
	}

	return;

	/*fprintf(spice, "* Simulation *\n");
	fprintf(spice, ".tran 5e-10 5e-8\n");

	fprintf(spice, "* MEASURE *\n");
	fprintf(spice, ".PRINT v(n3_476_165_25)\n");
	//fprintf(spice, ".MEASURE MAXVAL00001 MAX v(    1) FROM=0 TO= 5.000000e-08");
	fprintf(spice, "\n.end");*/
}

bool layer_specs::validate_values()
{
	bool width_gr_0  = (width > 0);
	bool pitch_gr_0  = (pitch > 0);
	bool offset_gr_0 = (offset > 0);
	bool shres_gr_0  = (sheet_resistance > 0);
	bool viares_gr_0 = (via_resistance > 0 || via_resistance == -1 );
	bool pitch_gr_width = ( pitch > width );
	return ( width_gr_0 && pitch_gr_0 && offset_gr_0 && shres_gr_0 &&
		viares_gr_0 && pitch_gr_width );
}

//! This method compares two nodes according to their y_coordinate
//! @param n1 [in]: first node to compare
//! @param n2 [in]: second node to compare
//! \return bool
static bool comp_ycoord(node* n1, node* n2)
{
	return (n1->get_ycoord() < n2->get_ycoord());
}

//! This method compares two nodes according to their x_coordinate
//! @param n1 [in]: first node to compare
//! @param n2 [in]: second node to compare
//! \return bool
static bool comp_xcoord(node* n1, node* n2)
{
	return (n1->get_xcoord() < n2->get_xcoord());
}

