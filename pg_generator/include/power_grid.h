/**
* \file power_grid.h \brief Power Grid Header File
*
* This file defines all the classes related to the power grid and their members
*/

#ifndef _POWER_GRID_H_
#define _POWER_GRID_H_

#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <sstream>
#include <map>
#include <algorithm>
#include "sparse_matrix.h"

class power_grid;
class sub_grid;
class layer;
class interconnect_tree;
class c4_pad;
class node;
class resistor;
class capacitor;
class inductor;
class current_source;
struct layer_specs;

// Layer Orientation
const bool HORIZONTAL = 0;
const bool VERTICAL   = 1;

// Grid Type
const bool INTERNAL   = 0;
const bool GLOBAL     = 1;

// Passive Elements Types
const int BRANCH  = 0;
const int PACKAGE = 1;
const int VIA     = 2;
const int CS      = 3;

// Node Types
const bool CHIP_NODE    = 0;
const bool PACKAGE_NODE = 1;

namespace common
{
	//! This method generates a random vector of booleans with specific number of ones
	//! @param size [in]: size of vector
	//! @param num_of_ones [in]: number of 1s in the vector
	//! \return vector<bool>
	extern vector<bool> random_mask_vector(int size, int num_of_ones);

	//! This method generates a random vector of doubles such that sum of all doubles
	//! in the vector is equal to 1.
	//! @param size [in]: size of vector
	//! \return vector<bool>
	extern vector<double> random_normalized_vector(unsigned int size);
}

namespace global_vars
{
	/*! verbose output on terminal */
	extern int verbose;
	/*! print grid info after generating? */
	extern int print_grid_info;
	/*! the config file  */
	extern string config_file;
	/*! the spice file for output*/
	extern string sp_file;
}

namespace rand_options
{
	extern bool _randomize_sub_grids;
	extern double _ratio_of_ON_sub_grids;
	extern double _ratio_of_retained_trees_in_sub_grids;
	extern double _ratio_of_retained_trees_in_intermediate_layers;
	extern double _ratio_of_retained_trees_in_global_layers;
	extern double _max_shift_tree_spans;
	extern int _seed_offset;
}

class power_grid
{
public:
	// Regular Constructor
	power_grid(string params_file);

	// Destructor
	~power_grid();

	// Getter functions
	vector<layer_specs* >   get_list_of_layer_specs();
	vector<vector<layer*> > get_list_of_layers();
	vector<sub_grid*>       get_list_of_sub_grids();
	vector<c4_pad*>         get_list_of_c4_pads();
	vector<resistor*>       get_list_of_vias();
	vector<current_source*> get_list_of_current_sources();
	vector<layer*>          get_layer(int index);
	layer*                  get_c4_layer();
	vector<layer*>          get_cs_layer();
	node*                   get_ground_node();
	node*                   get_vdd_node();

	unsigned int      get_number_of_nodes();
	unsigned int      get_number_of_layers();
	unsigned int      get_number_of_layers_global_grid();
	unsigned int      get_number_of_layers_sub_grid();
	unsigned int      get_number_of_c4s();
	unsigned int      get_number_of_sub_grids();
	unsigned int      get_number_of_sub_grids_X_DIM();
	unsigned int      get_number_of_sub_grids_Y_DIM();
	unsigned int      get_number_of_branches();
	unsigned int      get_number_of_current_sources();
	
	double            get_grid_dimension_x();
	double            get_grid_dimension_y();
	double            get_unit();
	double            get_vdd();
	double            get_delta_t();
	double            get_avg_power_density();
	double            get_total_power_dissipation_watts();
	double            get_ratio_current_sources_in_M1();
	double            get_C4_spacing_factor();
	double            get_C4_resistance();
	double            get_C4_inductance();

	matrix::sparse_matrix & get_G_matrix_ref();
	matrix::sparse_matrix   get_G_matrix();
	matrix::sparse_matrix & get_C_matrix_ref();
	matrix::sparse_matrix   get_C_matrix();
	matrix::sparse_matrix & get_A_matrix_ref();
	matrix::sparse_matrix   get_A_matrix();
	matrix::sparse_matrix & get_L_matrix_ref();
	matrix::sparse_matrix   get_L_matrix();
	matrix::sparse_matrix & get_M_matrix_ref();
	matrix::sparse_matrix   get_M_matrix();
	
	// Setter Functions
	void add_layer(layer* l);
	void add_sub_grid(sub_grid* sg);
	void add_c4_pad(c4_pad* c4);
	void add_via(resistor* via);
	void add_current_source(current_source* cs);
	void set_number_of_current_sources(int num_cs);
	
	// Other Functions
	void read_config_file(string params_file);
	void make_G_matrix();
	void make_C_matrix();
	void make_A_matrix();
	void make_L_matrix();
	void make_M_matrix();
	void generate_spice_file_for_resistive_grid( string filename );
	void eleminate_resistor(resistor * r); // Not implemented yet

	// Print Grid info
	void print_power_grid_info();

	// dellocate/free memory
	void  deallocate_memory();

protected:
	// Default Constructor: to be used only bu sub_grid class
	power_grid();

	vector<layer_specs* >   _list_of_layer_specs;
	vector<vector<layer*> > _list_of_layers;
	vector<c4_pad* >        _list_of_c4_pads;
	vector<resistor*>       _list_of_vias;
	vector<current_source*> _list_of_current_sources;

	unsigned int _number_of_nodes;
	unsigned int _number_of_layers;
	unsigned int _number_of_c4s;
	unsigned int _number_of_branches;
	unsigned int _number_of_current_sources;
	unsigned int _lowest_layer_index;
	// Grid Specs
	double _GRID_DIMENSION_X_UM;
	double _GRID_DIMENSION_Y_UM;
	double _unit;
	double _vdd;
	double _delta_t;
	double _power_scale_factor;

	// Capacitors Values
	double _MIN_DIE_CAPACITANCE_F;
	double _MAX_DIE_CAPACITANCE_F;

	// C4 Spacing
	double _C4_SPACING_FACTOR;

    // C4 resistance and C4 capacitance
    double _C4_RESISTANCE;
    double _C4_INDUCTANCE;

	// Power Density and ratio of M1 nodes to have current sources
	double _AVG_POWER_DESNITY_WATT_PER_CM_SQ;
	double _RATIO_CS_IN_M1;

	// Grid Matrices
	matrix::sparse_matrix _G;
	matrix::sparse_matrix _C;
	matrix::sparse_matrix _A;
	matrix::sparse_matrix _L;
	matrix::sparse_matrix _M;  // Incidence Matrix

	// a flag to indicate if dynamic memory has been deallocated
	bool _is_memory_deallocated;

	void create_global_grid();
	void create_internal_grid_islands();
	bool create_port_nodes_and_port_vias( int &res_idx );
	void add_new_nodes_for_C4();
	void sort_nodes_and_set_indices();
	void add_branch_resistances( int &res_idx );
	void add_capacitors( int &cap_idx );
	void add_c4_pads( int &res_idx, int &cap_idx );
	void finalize_node_name_and_indices();

private:
	node* _ground_node;
	node* _vdd_node;
	unsigned int _number_of_layers_global_grid;
	unsigned int _number_of_layers_sub_grid;
	unsigned int _number_of_sub_grids;
	unsigned int _number_of_sub_grids_X_DIM;
	unsigned int _number_of_sub_grids_Y_DIM;
	vector<bool> _ON_sub_grids; 
	sub_grid* _global_grid;
	vector<sub_grid* > _list_of_sub_grids;
};

struct layer_specs
{
	bool   orientation;
	double width;
	double pitch;
	double offset;
	double sheet_resistance;
	double via_resistance;

	bool validate_values();
};

class sub_grid: public power_grid
{
	public:
		// Regular Constructor
		sub_grid(power_grid* pg, bool grid_type, unsigned int index,
			double x_start, double x_final, double y_start, double y_final,
			double power_scale_factor = 1);
		// Destructor
		~sub_grid();
		//Getter Functions:
		double  get_x_start();
		double  get_x_final();
		double  get_y_start();
		double  get_y_final();
		int     get_index();
		int     get_base_index();
		//Setter Functions:
		void    set_x_start(double x_start);
		void    set_x_final(double x_final);
		void    set_y_start(double y_start);
		void    set_y_final(double y_final);
		void    set_index(int index);
		void    set_number_of_branches(int num_branches);
		void    increment_number_of_nodes();
	
		void make_G_matrix();
		void make_C_matrix();
	protected:
		// Default Constructor
		sub_grid();

		int    _base_index;
		double _x_start;
		double _x_final;
		double _y_start;
		double _y_final;
		int    _index;
		bool   _grid_type;
		power_grid* _grid_ptr;

		void build_layers();
		void build_interconnect_trees();
		void create_nodes_and_vias();
		void add_current_sources_to_M1();
};

class layer
{
public:
	// Default Constructor
	layer();

	// Regular Constructor
	layer(power_grid*  pg, 
		  layer_specs* specs,
		  string       name,
		  int          index,
		  bool         is_c4_layer,
		  bool         is_current_source_layer,
		  double       zcoord);

	// Destructor
	~layer();
	
	// Getter functions
	power_grid * get_power_grid_ptr();
	vector<interconnect_tree* >      get_list_of_interconnect_trees();
	vector<interconnect_tree* > &    get_list_of_interconnect_trees_ref();
	vector<node* >                   get_list_of_nodes();
	vector<resistor* >               get_list_of_resistors();
	map<double, interconnect_tree*>  get_it_map();
	layer_specs*                     get_layer_specs();
	
	string     get_name();
	int        get_index();
	double     get_width();
	double     get_pitch();
	double     get_offset();
	double     get_sheet_resistance();
	double     get_via_resistance();
	bool       is_c4_layer();
	bool       is_current_source_layer();
	bool       get_orientation();
   
	// Setter functions
	void       set_index(int index);
	void       set_name(string name);
	void       add_interconnect_tree(interconnect_tree* it);
	void       add_it_to_map(double xcoord, interconnect_tree* it);
	
	// Other Functions
	void       print();

protected:
	power_grid * _grid_ptr;

	vector<interconnect_tree* > _list_of_interconnect_trees;

	map<double, interconnect_tree*> _it_map;

	layer_specs* _specs;
	
	string       _name;
	int          _index;                         // Bottom layers have lower indices
	bool         _is_c4_layer;
	bool         _is_current_source_layer;
	double       _zcoord;                        // To keep track of the relative placement of layers
												 // = -1 if not available
};


class interconnect_tree
{
public:
	// Default Constructor
	interconnect_tree();

	// Regular Constructor
	interconnect_tree(layer* layer_ptr,
					  string name,
					  int    index,
					  double xcoord,
					  double ycoord);
	// Destructor
	~interconnect_tree();
	
	// Getter Functions
	layer*               get_layer_ptr();
	vector<node* >       get_list_of_nodes();
	vector<node* >     & get_list_of_nodes_ref();
	vector<resistor* >   get_list_of_resistors();
	vector<resistor* > & get_list_of_resistors_ref();
	map<double, node*>   get_nodes_map();
	string    get_name();
	int       get_index();
	double    get_xcoord();
	double    get_ycoord();

	// Setter Functions
	void      set_index(int index);
	void      set_name(string name);
	void      add_resistor(resistor * r);
	void      add_node(node* n);
	void      add_node_to_map(double coord, node *n);

	// Other Functions 
	bool      find_node(double xcoord, double ycoord);
	void      print();

protected:
	layer* _layer_ptr;
	
	vector<node* >     _list_of_nodes;
	vector<resistor* > _list_of_resistors;

	map<double, node*> _nodes_map;
	
	string _name;                                      
	int    _index;                                     // Unique ID
	double _xcoord;                                    // -1 if the layer is horizontal
	double _ycoord;                                    // -1 if the layer is vertical
};


class c4_pad
{
public:
	// Default Constructor
	c4_pad();

	// Regular Constructor
	c4_pad(power_grid * pg,
		   node       * port_node,
		   string       name,
		   int          index);

	// Destructor
	~c4_pad();
	
	// Getter Functions
	power_grid*            get_power_grid_ptr();
	node*                  get_port_node_ptr();
	vector<node* >         get_list_of_nodes();
	vector<node* >       & get_list_of_nodes_ref();
	vector<resistor* >     get_list_of_resistors();
	vector<resistor* >   & get_list_of_resistors_ref();
	vector<capacitor* >    get_list_of_capacitors();
	vector<capacitor* >  & get_list_of_capacitors_ref();
	vector<inductor* >     get_list_of_inductors();
	vector<inductor* >   & get_list_of_inductors_ref();
	string    get_name();
	int       get_index();

	// Setter Functions
	void      set_index(int index);
	void      set_name(string name);
	
	node*     add_node(string name, unsigned int & idx);
	void      add_resistor(string res_name, string node1, string node2,
		double value, int & idx, unsigned int & node_idx);
	void      add_capacitor(string cap_name, string node1, string node2,
		double value, int & idx, unsigned int & node_idx);
	void      add_inductor(string ind_name, string node1, string node2,
		double value, int & idx, unsigned int & node_idx);

	// Other Functions
	void      parse_subcircuit_file(unsigned int & node_idx, int & res_idx, int & cap_idx, int & ind_idx);
	void      print();
protected:
	power_grid* _power_grid;
	node*       _port_node;
	
	vector<node* >      _list_of_nodes;
	vector<resistor* >  _list_of_resistors;
	vector<capacitor* > _list_of_capacitors;
	vector<inductor* >  _list_of_inductors;
	   
	string _name;
	int _index;
};


class node
{
public:
	// Default Constructor
	node();

	// Regular Constructor
	node(interconnect_tree * it,
		 string name,
		 int    index,
		 double xcoord,
		 double ycoord,
		 bool   type);

	// Destructor
	~node();
	
	// Getter Functions
	interconnect_tree *    get_interconnect_tree_ptr(); 
	c4_pad*                get_c4_pad_connection();
	c4_pad*                get_c4_pad_ptr();
	vector<resistor* >     get_r_connections();
	vector<resistor* >   & get_r_connections_ref();
	vector<capacitor* >    get_c_connections();
	vector<capacitor* >  & get_c_connections_ref();
	vector<inductor* >     get_l_connections();
	vector<inductor* >   & get_l_connections_ref();
	current_source*        get_cs_connection();
	
	string                 get_name();
	int                    get_index();
	double                 get_xcoord();
	double                 get_ycoord();
	bool                   get_type();

	// Setter Functions
	void      set_name(string name);
	void      set_index(int index);
	void      add_resistor(resistor * r);
	void      add_capacitor(capacitor * c);
	void      add_inductor(inductor * l);
	void      add_current_source(current_source * cs);
	void      add_c4_pad_connection(c4_pad * c4);
	void      add_c4_pad_ptr(c4_pad * c4);

	// Other Functions
	vector<node* > get_list_of_neighboring_nodes();
	void           print();
protected:
	interconnect_tree* _interconnect_tree_ptr;
	c4_pad*            _c4_pad_connection;         // Used if the node is connectd to a C4 pad. NULL otherwise
	c4_pad*            _c4_pad_ptr;                // Used if the node is part of a C4 pad. NULL otherwise
	
	vector<resistor* >       _r_connections;       // Empty if not valid
	vector<capacitor* >      _c_connections;       // Empty if not valid
	vector<inductor* >       _l_connections;       // Empty if not valid
	current_source*          _cs_connection;       // Empty if not valid
  
	string _name;
	int    _index;          
	double _xcoord;                                // -1 if inside a C4 pad
	double _ycoord;                                // -1 if inside a C4 pad 
	bool   _type;                                  // can be PACKAGE_NODE (if part of C4 pad) or CHIP_NODE otherwise 
};


class resistor
{
public:
	// Default Constructor
	resistor();

	// Regular Constructor
	resistor(string name,
			 int    index,
			 double value,
			 double length,
			 double cs_area,
			 int    type);

	// Destructor
	~resistor();
	
	// Getter Functions
	node*        get_n0();
	node*        get_n1();
	string       get_name();
	int          get_index();
	double       get_value();
	double       get_length();
	double       get_cross_sectional_area();
	int          get_type();

	// Setter Functions
	void      set_name(string name);
	void      set_index(int index);
	void      set_value(double v);
	void      set_length(double l);
	void      set_cross_sectional_area(double cs_area);
	void      set_type(int type);
	void      set_nodes(node* n0, node* n1);

	// Other Functions
	void      print();
protected:
	node *_n0, *_n1;               // Pointers to the nodes to which the resistor is connected
								   // node with lower index is n0
	string _name;
	int    _index;
	double _value;
	double _length;
	double _cross_sectional_area;
	int    _type;                  // Can be VIA or PACKAGE or BRANCH 
								   // or CS if the this is part of a current source model
};


class capacitor
{
public:
	// Default Constructor
	capacitor();

	// Regular Constructor
	capacitor(string name,
			  int    index,
			  double value,
			  int    type);

	// Destructor
	~capacitor();
	
	// Getter Functions
	node*        get_n0();
	node*        get_n1();
	string       get_name();
	int          get_index();
	double       get_value();
	int          get_type();

	// Setter Functions
	void      set_name(string name);
	void      set_index(int index);
	void      set_value(double v);
	void      set_type(int type);
	void      set_nodes(node* n0, node* n1);

	// Other Functions
	void      print();
protected:
	node *_n0, *_n1;               // Pointers to the nodes to which the capacitor is connected

	string _name;
	int    _index;
	double _value;
	int    _type;                  // Can be PACKAGE or BRANCH
};


class inductor
{
public:
	// default constructor
	inductor();

	// regular constructor
	inductor(string name,
			 int    index,
			 double value,
			 int    type);

	// destructor
	~inductor();
	
	// getter functions
	node*        get_n0();
	node*        get_n1();
	string       get_name();
	int          get_index();
	double       get_value();
	int          get_type();

	// setter functions
	void      set_name(string name);
	void      set_index(int index);
	void      set_value(double v);
	void      set_type(int type);
	void      set_nodes(node* n0, node* n1);

	// Other Functions
	void      print();

protected:
	node *_n0, *_n1;               // Pointers to the nodes to which the inductor is connected

	string _name;
	int    _index;
	double _value;
	int    _type;                  // Can be PACKAGE or BRANCH
};


class current_source
{
public:
	// default constructor
	current_source();

	// regular constructor
	current_source(layer * layer_ptr,
				   string  name,
				   int     index,
				   double  value);
	// destructor
	~current_source();
	
	// getter functions
	layer*    get_layer_ptr();
	node*     get_n0();
	string    get_name();
	int       get_index();
	double    get_value();

	// setter functions
	void      set_name(string name);
	void      set_index(int index);
	void      set_node(node * n0);
	void      set_value(double value);
	
	// Other Functions 
	void      print();

protected:
	node*  _n0;
	layer* _layer_ptr;

	string _name;
	int    _index;

	double _value;
};
#endif




