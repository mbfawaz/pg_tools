/**
* \file sub_grid.cxx \brief Power Grid Implementation
*
* This file implements all the methods of the main class (sub_grid) and builds the sub grid
*/
#include "power_grid.h"
#include "profiler.h"
#include <algorithm>

using namespace std;
using namespace matrix;

// for facilitating profiler task descriptions
static string task_desc;
static string gtype_str;

// sub_grid Class

//! Default Constructor
sub_grid::sub_grid()
{
	_number_of_nodes = 0;
	_number_of_c4s = 0;
	_number_of_branches = 0;
	_number_of_layers = 0;
	_number_of_current_sources = 0;
	_unit = 1e-6;
	_index = 0;
	_x_start = 0;
	_y_start = 0;
	_x_final = 0;
	_y_final = 0;
}

//! Regular Constructor
//! @param params_file [in]: name of the file to read the parameters from
sub_grid::sub_grid(power_grid* pg, bool grid_type, unsigned int index, double x_start,
	double x_final, double y_start, double y_final, double power_scale_factor)
{
	// Initialize everything
	_delta_t = pow(10,-12); 
	_x_start = x_start;
	_x_final = x_final;
	_y_start = y_start;
	_y_final = y_final;
	_grid_ptr = pg;
	_grid_type = grid_type; // INTERNAL = 0, GLOBAL = 1
	_number_of_nodes = 0;
	_number_of_c4s = 0;
	_number_of_branches = 0;
	_number_of_current_sources = 0;
	_index = (_grid_type == INTERNAL) ? index: 0;
	_power_scale_factor = power_scale_factor;

	_GRID_DIMENSION_X_UM = _grid_type ? _grid_ptr->get_grid_dimension_x(): _x_final - _x_start ;
	_GRID_DIMENSION_Y_UM = _grid_type ? _grid_ptr->get_grid_dimension_y(): _y_final - _y_start ;
	_lowest_layer_index  = _grid_type ? _grid_ptr->get_number_of_layers() -
		_grid_ptr->get_number_of_layers_global_grid() + 1 : 1;
	_number_of_layers    = _grid_type ? _grid_ptr->get_number_of_layers_global_grid():
		_grid_ptr->get_number_of_layers() - _grid_ptr->get_number_of_layers_global_grid();
	_AVG_POWER_DESNITY_WATT_PER_CM_SQ = _grid_ptr->get_avg_power_density();
	_RATIO_CS_IN_M1 = _grid_type ? 0 : _grid_ptr->get_ratio_current_sources_in_M1();
	_vdd = _grid_ptr->get_vdd();

	gtype_str.assign( (_grid_type)?"Global Grid":"Internal Grid" );
	
	if ( global_vars::verbose )
		cerr<<endl<<gtype_str<<" #"<<_index<<" :";
	
	task_desc = "Building Layers";
	build_layers();
	if ( global_vars::verbose )
		cerr<<endl<<"   - "<<task_desc<<string( 55-task_desc.length(), ' ' )<<"done"<<" ("<<
			profiler::get_task_timer_data_by_desc( gtype_str + ": " + task_desc )->elapsed.wall/1e9/
			profiler::get_task_timer_data_by_desc( gtype_str + ": " + task_desc )->task_count<<" secs)";
	
	task_desc = "Build Interconnect Trees";
	build_interconnect_trees();
	if ( global_vars::verbose )
		cerr<<endl<<"   - "<<task_desc<<string( 55-task_desc.length(), ' ' )<<"done"<<" ("<<
			profiler::get_task_timer_data_by_desc( gtype_str + ": " + task_desc )->elapsed.wall/1e9/
			profiler::get_task_timer_data_by_desc( gtype_str + ": " + task_desc )->task_count<<" secs)";

	task_desc = "Create Nodes and Vias";
	create_nodes_and_vias();
	if ( global_vars::verbose )
		cerr<<endl<<"   - "<<task_desc<<string( 55-task_desc.length(), ' ' )<<"done"<<" ("<<
			profiler::get_task_timer_data_by_desc( gtype_str + ": " + task_desc )->elapsed.wall/1e9/
			profiler::get_task_timer_data_by_desc( gtype_str + ": " + task_desc )->task_count<<" secs)";
	
	if(_grid_type == INTERNAL)
	{
		task_desc = "Add Source Currents to M1";
		add_current_sources_to_M1();
		if ( global_vars::verbose )
			cerr<<endl<<"   - "<<task_desc<<string( 55-task_desc.length(), ' ' )<<"done"<<" ("<<
				profiler::get_task_timer_data_by_desc( gtype_str + ": " + task_desc )->elapsed.wall/1e9/
				profiler::get_task_timer_data_by_desc( gtype_str + ": " + task_desc )->task_count<<" secs)";
	}
}

//! Destructor of the sub grid
sub_grid::~sub_grid()
{
	// all sub-grid layers, vias and current sources are also added to
	// the corresponding lists in power_grid class. We delete it from 
	// the destrcutor of power_grid class.

	// the only thing to be deleted from here are the allocated layer_specs
	// for internal grids in build_layers(). This is because they are not
	// part of the _list_of_layer_specs in power_grid class.
	if(_grid_type == INTERNAL && rand_options::_randomize_sub_grids )
		for (unsigned int i = 0; i < _number_of_layers; i++)
			for(unsigned int j = 0; j < _list_of_layers[i].size(); j++)
				delete _list_of_layers[i].at(j)->get_layer_specs();

	// so that we dont print the warning message in base's destructor
	_is_memory_deallocated = true;
}

void sub_grid::build_layers()
{
	profiler record( gtype_str + ": " + task_desc );

	// reset the seed
	srand(rand_options::_seed_offset+7875);

	// number of layers in the power grid
	unsigned int pg_num_layers = _grid_ptr->get_number_of_layers();
	vector<layer_specs*> layers;
	if( _grid_type == INTERNAL && rand_options::_randomize_sub_grids )
	{
		double random_factor1 = 1, random_factor2 = 1;
		int index = get_index();
		if( index > 0 && index%2 == 0 )
		{
			random_factor1 = 0.9 + ((double)rand()/RAND_MAX)*0.2;
		}
		else if ( index > 0 && index%2 == 1 )
		{
			random_factor2 = 0.7 + ((double)rand()/RAND_MAX)*0.5;
		}//*/

		layers.resize(_number_of_layers);
		vector<layer_specs*> pg_layer_spec = _grid_ptr->get_list_of_layer_specs();

		//printf("%d\n", _index);
		//if (index == 1 || index == 3 || index == 4)
		//	random_factor2 = 2;

		for(unsigned int i = 0 ; i < _number_of_layers ; i++)
		{

			layer_specs* this_layer_spec = pg_layer_spec.at(pg_num_layers - i-1);
			layers[i] = new layer_specs();
			layers[i]->orientation      = this_layer_spec -> orientation;
			layers[i]->offset           = this_layer_spec -> offset;
			layers[i]->sheet_resistance = this_layer_spec -> sheet_resistance*random_factor1;
			layers[i]->via_resistance   = this_layer_spec -> via_resistance*random_factor1;
			layers[i]->pitch            = this_layer_spec -> pitch*random_factor2;
			layers[i]->width            = this_layer_spec -> width*random_factor2;
			/*layers[i]->pitch            = this_layer_spec -> pitch;
			layers[i]->width            = this_layer_spec -> width*random_factor2;*/
		}
	}

	_list_of_layers.resize(_number_of_layers);
	for(unsigned int i = 0 ; i < _number_of_layers ; i++)
		_list_of_layers.at(i).resize(0);
	
	// Building Layers
	string layer_name;
	for (unsigned int i = _lowest_layer_index; i < _number_of_layers + _lowest_layer_index; i++)
	{
		stringstream ss; 
		ss << "M" << i;
		layer_name = ss.str();
		
		bool is_c4_layer = false, is_cs_layer = false;
		if (i == 1)
			is_cs_layer = true;
		else if (i == pg_num_layers )// || i == _number_of_layers)
			is_c4_layer = true;
		
		layer* l;
		if( _grid_type == GLOBAL || !rand_options::_randomize_sub_grids )
			l = new layer(_grid_ptr, _grid_ptr->get_list_of_layer_specs().at(pg_num_layers - i),
				layer_name, i, is_c4_layer, is_cs_layer, -1);
		else
			l = new layer(_grid_ptr, layers[i-1], layer_name, i, is_c4_layer, is_cs_layer, -1);
		add_layer(l);
		_grid_ptr->add_layer(l);
	}
}

void sub_grid::build_interconnect_trees()
{
	profiler record( gtype_str + ": " + task_desc );

	// set seed for randomness
	srand(rand_options::_seed_offset+45);

	// Building interconnect trees for this subgrid
	int it_counter = 0;
	interconnect_tree *it;
	int pg_top_layer_idx = _grid_ptr->get_c4_layer()->get_index();
	
	for (unsigned int i = 0; i < _number_of_layers ; i++)
	//for (int i = _number_of_layers-1; i >= 0; i--)
	{
		unsigned int number_of_interconnect_trees;
		layer* this_layer = _list_of_layers.at(i).at(0);
		int lidx = this_layer->get_index();

		// extra code
		int clayer_idx = i-2;
		bool short_rmv = ( clayer_idx >= 0 );

		if (this_layer->get_orientation() == HORIZONTAL)
		{
			number_of_interconnect_trees = 
				( (_GRID_DIMENSION_Y_UM - this_layer->get_offset() )
				/this_layer->get_pitch() ) + 1;
		}
		else
		{
			number_of_interconnect_trees =
				( (_GRID_DIMENSION_X_UM - this_layer->get_offset() )
				/this_layer->get_pitch() ) + 1;
		}

		// introduce randomness by removing some trees
		vector<bool> rmv;
		int num_retained_trees = 0;
		if ( _grid_type == INTERNAL )
		{
			num_retained_trees = (int)round(
				rand_options::_ratio_of_retained_trees_in_sub_grids*
				number_of_interconnect_trees);
		}
		else if ( lidx == pg_top_layer_idx || lidx == pg_top_layer_idx-1 )
		{
			num_retained_trees = (int)round(
				rand_options::_ratio_of_retained_trees_in_global_layers*
				number_of_interconnect_trees);
		}
		else
		{
			num_retained_trees = (int)round(
				rand_options::_ratio_of_retained_trees_in_intermediate_layers*
				number_of_interconnect_trees);
		}
		if ( num_retained_trees == 0 )
		{
			fprintf(stdout, "ERROR: Number of retained trees for layer %s is 0. Exiting now!",
				(char*)this_layer->get_name().c_str() );
			exit(1);
		}

		rmv = common::random_mask_vector( number_of_interconnect_trees,
			num_retained_trees );
		double distance_threshold = 0.01*this_layer->get_pitch();
		//printf("This layer = %s, pitch = %g, threshold = %g\n",
		//	(char*)this_layer->get_name().c_str(), this_layer->get_pitch(), distance_threshold );

		//int count;
		if (this_layer->get_orientation() == HORIZONTAL)
		{
			vector<double> ycoord_clayer;
			//count = 0;
			if (short_rmv)
			{
				layer* clayer = _list_of_layers.at(clayer_idx).at(0);
				//printf("Compare layer = %s\n", (char*)clayer->get_name().c_str());
				vector<interconnect_tree*> clayer_itree_list =
					clayer->get_list_of_interconnect_trees();
				ycoord_clayer.resize( clayer_itree_list.size() );
				for ( unsigned int j = 0; j < clayer_itree_list.size(); j++ )
					ycoord_clayer[j] = clayer_itree_list.at(j)->get_ycoord();
			}

			double xcoord = -1;
			double ycoord = this_layer->get_offset() + _y_start;
			string it_name;
			for (unsigned int j = 0; j < number_of_interconnect_trees; j++)
			{
				if ( !rmv[j]  )
				{
					ycoord += this_layer->get_pitch();
					continue;
				}
				if ( short_rmv )
				{
					double min_diff = 1e12;
					for ( unsigned int j = 0; j < ycoord_clayer.size(); j++ )
					{
						double diff = fabs( ycoord_clayer[j]-ycoord );
						if (diff < min_diff)	min_diff = diff;
					}
					if ( min_diff < distance_threshold )
					{
						ycoord += this_layer->get_pitch();
						//count++;
						//printf("%d. md = %g, ", count, min_diff);
						continue;
					}
				}
				//*/
				it_counter++;
				stringstream ss;
				ss << "IT" << it_counter;
				it_name = ss.str();
				it = new interconnect_tree(this_layer, it_name,
					it_counter, xcoord, ycoord);
				this_layer->add_interconnect_tree(it);
				this_layer->add_it_to_map(ycoord, it);
				ycoord += this_layer->get_pitch();
			}
		}
		else
		{
			vector<double> xcoord_clayer;
			//count = 0;
			if (short_rmv)
			{
				layer* clayer = _list_of_layers.at(clayer_idx).at(0);
				//printf("Compare layer = %s\n", (char*)clayer->get_name().c_str());
				vector<interconnect_tree*> clayer_itree_list =
					clayer->get_list_of_interconnect_trees();
				xcoord_clayer.resize( clayer_itree_list.size() );
				for ( unsigned int j = 0; j < clayer_itree_list.size(); j++ )
					xcoord_clayer[j] = clayer_itree_list.at(j)->get_xcoord();
			}

			double xcoord = this_layer->get_offset() + _x_start;
			double ycoord = -1;
			string it_name;
			for (unsigned int j = 0; j < number_of_interconnect_trees; j++)
			{
				if ( !rmv[j] )
				{
					xcoord += this_layer->get_pitch();
					continue;
				}
				if ( short_rmv )
				{
					double min_diff = 1e12;
					for ( unsigned int j = 0; j < xcoord_clayer.size(); j++ )
					{
						double diff = fabs( xcoord_clayer[j]-xcoord );
						if (diff < min_diff)	min_diff = diff;
					}
					if ( min_diff < distance_threshold )
					{
						xcoord += this_layer->get_pitch();
						//count++;
						//printf("%d. md = %g, ", count, min_diff);
						continue;
					}
				}
				//*/
				it_counter++;
				stringstream ss;
				ss << "IT" << it_counter;
				it_name = ss.str();
				it = new interconnect_tree(this_layer, it_name,
					it_counter, xcoord, ycoord);
				this_layer->add_interconnect_tree(it);
				this_layer->add_it_to_map(xcoord, it);
				xcoord += this_layer->get_pitch();
			}
		}
		//printf("\nRemoved %d trees\n-----\n", count);
	}
}

void sub_grid::create_nodes_and_vias()
{
	// set seed again
	srand(rand_options::_seed_offset+52);

	// Creating Nodes
	profiler record( gtype_str + ": " + task_desc );
	int node_idx = 0, res_idx = 0;
	string node_name, via_name;
	node *n1, *n2;
	resistor * via;
	double via_resistance;
	double xcoord, ycoord;
	
	// new code
	/*if ( _number_of_layers > 2 )
		for (int i = _number_of_layers - 1; i > 0; i--)
		{
			int next_next_layer_idx = i-2;
			if ( i - 2 < 0 )	continue;

			layer* this_layer = _list_of_layers.at(i).at(0);
			vector<interconnect_tree*> this_itree_list = this_layer->get_list_of_interconnect_trees();

			layer* nnext_layer = _list_of_layers.at(next_next_layer_idx).at(0);
			vector<interconnect_tree*> nnext_itree_list = nnext_layer->get_list_of_interconnect_trees();

			assert( this_layer->get_orientation() == nnext_layer->get_orientation() );

			double distance_threshold = 0.01*_list_of_layers.at(i).at(0)->get_pitch();
			if (this_layer->get_orientation() == HORIZONTAL)
			{
				vector<bool> to_delete( this_itree_list.size(), false );
				for (unsigned int j = 0; j < this_itree_list.size(); j++)
				{
					ycoord = this_itree_list.at(j)->get_ycoord();
					for (unsigned int k = 0; k < nnext_itree_list.size(); k++)
					{
						double ycoord1 = nnext_itree_list.at(j)->get_ycoord();
						if ( fabs(ycoord1 - ycoord) < distance_threshold  )
						{
							to_delete[j] = true;
						}
					}
				}



			}
			else
			{

			}
		}//*/

	//int pg_top_layer_idx = _grid_ptr->get_c4_layer()->get_index();
	for (int i = _number_of_layers - 1; i > 0; i--)
	{
		layer* this_layer = _list_of_layers.at(i).at(0);
		vector<interconnect_tree*> this_itree_list = this_layer->get_list_of_interconnect_trees();

		layer* next_layer = _list_of_layers.at(i-1).at(0);
		vector<interconnect_tree*> next_itree_list = next_layer->get_list_of_interconnect_trees();
		via_resistance = next_layer->get_via_resistance();

		int diff = this_layer->get_index() - next_layer->get_index();
		assert( diff == 1 || diff == -1 );

		if (this_layer->get_orientation() == HORIZONTAL)
		{
			for (unsigned int j = 0; j < this_itree_list.size(); j++)
			{
				ycoord = this_itree_list.at(j)->get_ycoord();

				// introduce randomness by staring and ending at arbitrary trees
				int start_k = 0, end_k = next_itree_list.size();
				int next_itree_list_size = next_itree_list.size();
				int start_k_margin = (int)( rand_options::_max_shift_tree_spans*
					next_itree_list_size );
				if ( start_k_margin != 0 )
				{
					start_k = rand() % start_k_margin;
					int end_k_margin = (int)( rand_options::_max_shift_tree_spans*
						next_itree_list_size );
					end_k = next_itree_list_size - ( rand() % end_k_margin );
					assert( end_k - start_k >= 1 );
				}
			
				//for (unsigned int k = 0; k < next_itree_list.size(); k++)
				for (int k = start_k; k < end_k; k++)
				{
					xcoord = next_itree_list.at(k)->get_xcoord();
					
					res_idx++;
					stringstream ss;
					ss << "R" << res_idx << "_" << _index;
					via_name = ss.str();
					via = new resistor(via_name, res_idx, via_resistance, -1, -1, VIA);

					node_idx++;
					stringstream ss1;
					ss1 << "n" << node_idx;
					node_name = ss1.str();
					n1 = new node(this_itree_list.at(j),
						node_name, node_idx, xcoord, ycoord, CHIP_NODE);
					
					node_idx++;
					stringstream ss2;
					ss2 << "n" << node_idx;
					node_name = ss2.str();
					n2 = new node(next_itree_list.at(k),
						node_name, node_idx, xcoord, ycoord, CHIP_NODE);

					// Adding Via resistors
					via->set_nodes(n1, n2);
					
					n1->add_resistor(via);
					this_itree_list.at(j)->add_node(n1);
					this_itree_list.at(j)->add_node_to_map(xcoord, n1);
					
					n2->add_resistor(via);
					next_itree_list.at(k)->add_node(n2);
					next_itree_list.at(k)->add_node_to_map(ycoord, n2);
					
					add_via(via);
					_grid_ptr->add_via(via);
				}
			}
		}
		else
		{
			for (unsigned int j = 0; j < this_itree_list.size(); j++)
			{
				xcoord = this_itree_list.at(j)->get_xcoord();

				// introduce randomness by staring and ending at arbitrary trees
				int start_k = 0, end_k = next_itree_list.size();
				int next_itree_list_size = next_itree_list.size();
				int start_k_margin = (int)( rand_options::_max_shift_tree_spans*
					next_itree_list_size );
				if ( start_k_margin != 0 )
				{
					start_k = rand() % start_k_margin;
					int end_k_margin = (int)( rand_options::_max_shift_tree_spans*
						next_itree_list_size );
					end_k = next_itree_list_size - ( rand() % end_k_margin );
					assert( end_k - start_k >= 1 );
				}

				//for (unsigned int k = 0; k < next_itree_list.size(); k++)
				for (int k = start_k; k < end_k; k++)
				{
					ycoord = next_itree_list.at(k)->get_ycoord();
					
					res_idx++;
					stringstream ss;
					ss << "R" << res_idx << "_" << _index;
					via_name = ss.str();
					via = new resistor(via_name, res_idx, via_resistance, -1, -1, VIA);

					node_idx++;
					stringstream ss1;
					ss1 << "n" << node_idx;
					node_name = ss1.str();
					n1 = new node(this_itree_list.at(j),
						node_name, node_idx, xcoord, ycoord, CHIP_NODE);

					node_idx++;
					stringstream ss2;
					ss2 << "n" << node_idx;
					node_name = ss2.str();
					n2 = new node(next_itree_list.at(k),
						node_name, node_idx, xcoord, ycoord, CHIP_NODE);
					
					// Adding Via resistors
					via->set_nodes(n1, n2);
					
					n1->add_resistor(via);
					this_itree_list.at(j)->add_node(n1);
					this_itree_list.at(j)->add_node_to_map(ycoord, n1);

					n2->add_resistor(via);
					next_itree_list.at(k)->add_node(n2);
					next_itree_list.at(k)->add_node_to_map(xcoord, n2);

					add_via(via);
					_grid_ptr->add_via(via);
				}
			}
		}
	}
	_number_of_nodes = node_idx;
}

void sub_grid::add_current_sources_to_M1()
{
	// Adding Current Sources
	profiler record( gtype_str + ": " + task_desc );
	srand(rand_options::_seed_offset+6); // re-init the random number generator

	layer* M1 = _list_of_layers.at(_number_of_layers - 1).at(0);
	vector<interconnect_tree*> M1_itree_list = M1->get_list_of_interconnect_trees();
	
	/*  POLICY
			- Each power block, if possible, needs to have at least 10 trees.
			- Allow a min of 1 and a max of 10 power blocks per subgrid
	*/
	const unsigned int min_itrees_per_block = 10, max_sg_power_blocks = 10;

	//STEP 1: get total pwoer disspationfor this ubgrid in watts
	double tpd = get_total_power_dissipation_watts();
	
	// STEP 2: find the number of power blocks in the subgrid, and number of
	// trees per block 
	unsigned int num_sg_power_blocks = 1; 
	if ( M1_itree_list.size() > 2*min_itrees_per_block )
		num_sg_power_blocks = floor( M1_itree_list.size()/min_itrees_per_block );
	if ( num_sg_power_blocks > max_sg_power_blocks )
		num_sg_power_blocks = max_sg_power_blocks;
	unsigned int num_itrees_per_block = M1_itree_list.size()/num_sg_power_blocks;	

	// STEP 3: find normalzied random numbers to distribute the total power
	// among the blocks
	vector<double> psf_vec = common::random_normalized_vector( num_sg_power_blocks );

	unsigned int cs_idx = 0;
	string cs_name;
	current_source *cs;
	
	for( unsigned int i = 0; i < num_sg_power_blocks; i++ )
	{
		// STEP 4: get a list of all nodes within this block
		vector<node*> block_node_list;
		for (unsigned int j = i*num_itrees_per_block; j < (i+1)*num_itrees_per_block-1; j++)
		{
			vector<node*> itree_node_list = M1_itree_list[j]->get_list_of_nodes();
			block_node_list.insert( block_node_list.end(), itree_node_list.begin(),
				itree_node_list.end() );
		}
		// Adding any remaining trees to the last block 
		if ( i == num_sg_power_blocks-1 )
		{
			for (unsigned int j = (i+1)*num_itrees_per_block; j < M1_itree_list.size(); j++)
			{
				vector<node*> itree_node_list = M1_itree_list[j]->get_list_of_nodes();
				block_node_list.insert( block_node_list.end(), itree_node_list.begin(),
					itree_node_list.end() );
			}
		}

		// STEP 5: Find the number of current sources to be put and the value of each current
		// source with
		unsigned int num_cs_nodes = round( _RATIO_CS_IN_M1*block_node_list.size() );
		assert(num_cs_nodes != 0);
		double cs_value = (psf_vec[i]*tpd)/(_vdd*num_cs_nodes);
		vector<bool> node_has_cs = common::random_mask_vector( block_node_list.size(), num_cs_nodes);

		for(unsigned int j = 0; j < block_node_list.size(); j++ )
		{
			if ( !node_has_cs[j] )		continue;
			
			cs_idx++;
			stringstream ss;
			ss << "i" << cs_idx << "_" << _index << "_v"; //_g for ground
			cs_name = ss.str();
			cs = new current_source( M1, cs_name, cs_idx, cs_value);
			cs->set_node( block_node_list.at(j) );
			block_node_list.at(j)->add_current_source(cs);
			add_current_source(cs);
			_grid_ptr->add_current_source(cs);
		}
	}
	_number_of_current_sources = cs_idx;
	fprintf(stdout, "\n -- Added %d current sources to subgrid#%d, TPD = %g W, psf = %g\n",
		_number_of_current_sources, get_index(), get_total_power_dissipation_watts(),
		_power_scale_factor );
	_grid_ptr->set_number_of_current_sources( _grid_ptr->get_number_of_current_sources()
		+ _number_of_current_sources );
}

//! This method returns the base index of the subrid, i.e. the smallest index of all nodes in the sub grid
//! \return int
int sub_grid::get_base_index()
{
	vector<interconnect_tree*> itree_list =
		get_layer(_number_of_layers).at(0)->get_list_of_interconnect_trees();
	int num_of_nodes_in_highest_layer =	itree_list[0]->get_list_of_nodes().size();
	int num_of_IT_in_highest_layer = itree_list.size();
	
	int base_index = itree_list[num_of_IT_in_highest_layer-1]->
		get_list_of_nodes()[num_of_nodes_in_highest_layer-1]->get_index()-1;
	
	_base_index = base_index;
	return _base_index;
}

//! This method returns the index of the subrid
//! \return int
int sub_grid::get_index()
{
	return _index;
}
//! This method returns the x_start of the subrid
//! \return double
double sub_grid::get_x_start()
{
	return _x_start;
}
//! This method returns the x_final of the subrid
//! \return double
double sub_grid::get_x_final()
{
	return _x_final;
}
//! This method returns the y_start of the subrid
//! \return double
double sub_grid::get_y_start()
{
	return _y_start;
}
//! This method returns the y_final of the subrid
//! \return double
double sub_grid::get_y_final()
{
	return _y_final;
}

//! This method builds the G matrix of the sub grid 
void sub_grid::make_G_matrix()
{
	fprintf(stdout, "\033[1;%dmBuilding G Matrix for sub grid\033[0m\n",36);
	profiler record("Building G matrix");
	_G.set_number_of_rows(_number_of_nodes);
	_G.set_number_of_columns(_number_of_nodes);
	//cout<<"Num of nodes"<<_number_of_nodes<<endl;
	int num_res = 0;
	int num_port_res = 0;
	// Going over all branches in the grid (excluding C4s)
	vector<resistor*> branches;
	// The base index is the smallest index in an internal grid
	int base_index = get_base_index();
	for (int i = _number_of_layers-1; i >= 0; i--)
	{
		branches = _list_of_layers.at(i).at(0)->get_list_of_resistors();
		num_res += branches.size();
		for (unsigned int j = 0; j < branches.size(); j++)
		{
			//cout<<branches[j]->get_n0()->get_index()-1-base_index<<endl;
			//cout<<branches[j]->get_n1()->get_index()-1-base_index<<endl;
			_G.get_Ti().push_back(branches[j]->get_n0()->get_index()-1-base_index);
			_G.get_Tj().push_back(branches[j]->get_n0()->get_index()-1-base_index);
			_G.get_Tx().push_back(1/branches[j]->get_value());

			_G.get_Ti().push_back(branches[j]->get_n1()->get_index()-1-base_index);
			_G.get_Tj().push_back(branches[j]->get_n1()->get_index()-1-base_index);
			_G.get_Tx().push_back(1/branches[j]->get_value());

			_G.get_Ti().push_back(branches[j]->get_n0()->get_index()-1-base_index);
			_G.get_Tj().push_back(branches[j]->get_n1()->get_index()-1-base_index);
			_G.get_Tx().push_back(-1/branches[j]->get_value());

			_G.get_Ti().push_back(branches[j]->get_n1()->get_index()-1-base_index);
			_G.get_Tj().push_back(branches[j]->get_n0()->get_index()-1-base_index);
			_G.get_Tx().push_back(-1/branches[j]->get_value());
		}
	}
	branches.clear();
	//cout<<"------------------------"<<endl;
	// Going over all the vias
	num_res += _list_of_vias.size();
	for (unsigned int i = 0; i < _list_of_vias.size(); i++)
	{

		if(_list_of_vias[i]->get_n0()->get_interconnect_tree_ptr()->get_layer_ptr()->get_index() >
			 _number_of_layers)
		{
			//cout<<"n0 is outside subgird"<<endl;
			num_port_res++;
			//cout<<_list_of_vias[i]->get_n1()->get_index()-1-base_index<<endl;
			_G.get_Ti().push_back(_list_of_vias[i]->get_n1()->get_index()-1-base_index);
			_G.get_Tj().push_back(_list_of_vias[i]->get_n1()->get_index()-1-base_index);
			_G.get_Tx().push_back(1/_list_of_vias[i]->get_value());
			continue;
		}
		if(_list_of_vias[i]->get_n1()->get_interconnect_tree_ptr()->get_layer_ptr()->get_index() >
			_number_of_layers)
		{
			//cout<<"n1 is outside subgird"<<endl;
			num_port_res++;
			//cout<<_list_of_vias[i]->get_n0()->get_index()-1-base_index<<endl;
			_G.get_Ti().push_back(_list_of_vias[i]->get_n0()->get_index()-1-base_index);
			_G.get_Tj().push_back(_list_of_vias[i]->get_n0()->get_index()-1-base_index);
			_G.get_Tx().push_back(1/_list_of_vias[i]->get_value());
			continue;
		}

		_G.get_Ti().push_back(_list_of_vias[i]->get_n0()->get_index()-1-base_index);
		_G.get_Tj().push_back(_list_of_vias[i]->get_n0()->get_index()-1-base_index);
		_G.get_Tx().push_back(1/_list_of_vias[i]->get_value());
				   
		_G.get_Ti().push_back(_list_of_vias[i]->get_n1()->get_index()-1-base_index);
		_G.get_Tj().push_back(_list_of_vias[i]->get_n1()->get_index()-1-base_index);
		_G.get_Tx().push_back(1/_list_of_vias[i]->get_value());
				   
		_G.get_Ti().push_back(_list_of_vias[i]->get_n0()->get_index()-1-base_index);
		_G.get_Tj().push_back(_list_of_vias[i]->get_n1()->get_index()-1-base_index);
		_G.get_Tx().push_back(-1/_list_of_vias[i]->get_value());
				   
		_G.get_Ti().push_back(_list_of_vias[i]->get_n1()->get_index()-1-base_index);
		_G.get_Tj().push_back(_list_of_vias[i]->get_n0()->get_index()-1-base_index);
		_G.get_Tx().push_back(-1/_list_of_vias[i]->get_value());

		//cout<<_list_of_vias[i]->get_n0()->get_index()-1-base_index<<endl;
		//cout<<_list_of_vias[i]->get_n1()->get_index()-1-base_index<<endl;
	}
	//cout<<"Number of nonz = "<<_number_of_nodes + 2*num_res<<endl;
	_G.set_nz(_number_of_nodes + 2*num_res - 2*num_port_res);

	//cout<<_G.get_Ti()<<endl;
	//cout<<_G.get_Tj()<<endl;
	//cout<<_G.get_Tx()<<endl;
	_G.get_column_format();
	branches.clear();
}

//! This method builds the C matrix of the sub grid 
void sub_grid::make_C_matrix()
{
	fprintf(stdout, "\033[1;%dmBuilding C Matrix for sub grid\033[0m\n",36);
	profiler record("Building C matrix");
	_C.set_number_of_rows(_number_of_nodes);
	_C.set_number_of_columns(_number_of_nodes);
	_C.set_nz(_number_of_nodes);

	int num_of_nodes_in_highest_layer = get_layer(_number_of_layers).at(0)->get_list_of_interconnect_trees()[0]->get_list_of_nodes().size();
	int num_of_IT_in_highest_layer = get_layer(_number_of_layers).at(0)->get_list_of_interconnect_trees().size();
	// The base index is the smallest index in an internal grid
	int base_index = get_layer(_number_of_layers).at(0)->get_list_of_interconnect_trees()[num_of_IT_in_highest_layer-1]->get_list_of_nodes()[num_of_nodes_in_highest_layer-1]->get_index()-1;
	// Going over all capacitors in the grid (excluding C4s)
	vector<node*> nodes;
	for (unsigned int i = 0; i < _number_of_layers; i++)
	{
		for (unsigned int j = 0; j < _list_of_layers.at(i).at(0)->get_list_of_interconnect_trees().size(); j++)
		{
			nodes  = _list_of_layers.at(i).at(0)->get_list_of_interconnect_trees().at(j)->get_list_of_nodes_ref();
			for (unsigned int k = 0; k < nodes.size(); k++)
			{
				for (unsigned int m = 0; m < nodes.at(k)->get_c_connections().size(); m++)
				{
					_C.get_Ti().push_back(nodes.at(k)->get_index() - 1-base_index);
					_C.get_Tj().push_back(nodes.at(k)->get_index() - 1-base_index);
					_C.get_Tx().push_back(nodes.at(k)->get_c_connections().at(m)->get_value());
				}
			}
		}
	}

	_C.get_column_format();
}
//! This method sets the index of the subgrid
//! @param subgrid [in]: index of the subgrid
void sub_grid::set_index(int index)
{
	_index = index;
}
//! This method increments the number of nodes of subgrid
//! @param subgrid [in]: NONE
void sub_grid::increment_number_of_nodes()
{
	_number_of_nodes++;
}
//! This method sets the number of branches of subgrid
//! @param subgrid [in]: number of branches
void sub_grid::set_number_of_branches(int num_branches)
{
	_number_of_branches = num_branches;
}

//! This method sets the x_start of the subrid
//! @param subgrid [in]: x_start of the subgrid
void sub_grid::set_x_start(double x_start)
{
	_x_start = x_start;
}
//! This method sets the x_final of the subrid
//! @param subgrid [in]: x_final of the subgrid
void sub_grid::set_x_final(double x_final)
{
	_x_final = x_final;
}
//! This method sets the y_start of the subrid
//! @param subgrid [in]: y_start of the subgrid
void sub_grid::set_y_start(double y_start)
{
	_y_start = y_start;
}
//! This method sets the y_final of the subrid
//! @param subgrid [in]: y_final of the subgrid
void sub_grid::set_y_final(double y_final)
{
	_y_final = y_final;
}
