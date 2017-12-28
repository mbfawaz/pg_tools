/**
* \file layer.cxx \brief Layer Implementation
*
* This file implements the methods of the class layer
*/

#include "power_grid.h"

// Layer Class

// ! Default Constructor
layer::layer()
{
    _grid_ptr                = NULL;
    _specs                   = NULL;
    _name                    = "";
    _index                   = -1;
    _is_c4_layer             = false;
    _is_current_source_layer = false;
    _zcoord                  = -1;
}

// ! Regular Constructor
// ! @param pg                      [in]: pointer to the power grid that the layer is in
// ! @param specs                   [in]: specifications of the layer (orientation, width, ... etc)
// ! @param name                    [in]: name of the layer
// ! @param index                   [in]: index of the layer
// ! @param is_c4_layer             [in]: flag that indicates the layer is part of a C4 block or not
// ! @param is_current_source_layer [in]: flag that indicates the layer contains current sources or not
// ! @param zcoord                  [in]: z_coordinate of the layer
layer::layer(power_grid*  pg, 
             layer_specs* specs,
             string       name,
             int          index,
             bool         is_c4_layer,
             bool         is_current_source_layer,
             double       zcoord)
{
    _grid_ptr                = pg;
    _specs                   = specs;
    _name                    = name;
    _index                   = index;
    _is_c4_layer             = is_c4_layer;
    _is_current_source_layer = is_current_source_layer;
    _zcoord                  = zcoord;
}

//! Destructor of class layer
layer::~layer()
{
    // delete only inyerconnec ttrees. layer specs are deleted in pwoer_grid
    // and sub_grid classes
    for (unsigned int i = 0; i < _list_of_interconnect_trees.size(); i++)
        delete _list_of_interconnect_trees[i];
}

// Getter Functions

//! This method returns a pointer to the power grid that the layer is in
//! \return power_grid*
power_grid* layer::get_power_grid_ptr()
{
    return _grid_ptr;
}

//! This method returns a copy of the list of interconnect trees that are in the layer
//! \return vector<interconnect_tree*>
vector<interconnect_tree*> layer::get_list_of_interconnect_trees()
{
    return _list_of_interconnect_trees;
}

//! This method returns a reference copy of the list of interconnect trees that are in the layer
//! \return vector<interconnect_tree*> &
vector<interconnect_tree*>& layer::get_list_of_interconnect_trees_ref()
{
    return _list_of_interconnect_trees;
}

//! This method returns a copy of the list of nodes that are in the layer
//! \return vector<node*>
vector<node*> layer::get_list_of_nodes()
{
    vector<node*> nodes;
    unsigned int size = 0;
    unsigned int loc  = 0;
    for (unsigned int i = 0; i < _list_of_interconnect_trees.size(); i++)
    {
        size += _list_of_interconnect_trees.at(i)->get_list_of_nodes().size();
        nodes.resize(size);
        for (unsigned int j = 0; j < _list_of_interconnect_trees.at(i)->get_list_of_nodes().size(); j++)
            nodes.at(loc + j) = _list_of_interconnect_trees.at(i)->get_list_of_nodes().at(j);
        loc = size;
    }
    return nodes;
}

//! This method returns a copy of the list of resistors that are in the layer
//! \return vector<resistor*>
vector<resistor*> layer::get_list_of_resistors()
{
    vector<resistor*> resistors;
    unsigned int size = 0;
    unsigned int loc  = 0;
    for (unsigned int i = 0; i < _list_of_interconnect_trees.size(); i++)
    {
        size += _list_of_interconnect_trees.at(i)->get_list_of_resistors().size();
        resistors.resize(size);
        for (unsigned int j = 0; j < _list_of_interconnect_trees.at(i)->get_list_of_resistors().size(); j++)
            resistors.at(loc + j) = _list_of_interconnect_trees.at(i)->get_list_of_resistors().at(j);
        loc = size;
    }
    return resistors;
}

//! This method returns a copy of the interconnect tress map in the layer
//! \return map<double, interconnect_tree*>
map<double, interconnect_tree*> layer::get_it_map()
{
    return _it_map;
}
//! This method returns the specifications of the layer
//! \return layer_specs*
layer_specs* layer::get_layer_specs()
{
    return _specs;
}

//! This method returns the name of the layer
//! \return string
string layer::get_name()
{
    return _name;
}

//! This method returns the index of the layer
//! \return int
int layer::get_index()
{
    return _index;
}

//! This method returns the width of the layer
//! \return double
double layer::get_width()
{
    return _specs->width;
}

//! This method returns the pitch of the layer
//! \return double
double layer::get_pitch()
{
    return _specs->pitch;
}

//! This method returns the offset of the layer
//! \return double
double layer::get_offset()
{
    return _specs->offset;
}
//! This method returns the sheet resistance of the layer
//! \return double
double layer::get_sheet_resistance()
{
    return _specs->sheet_resistance;
}

//! This method returns the via resistance of the layer
//! \return double
double layer::get_via_resistance()
{
    return _specs->via_resistance;
}

//! This method returns a flag that indicates whether the layer is part of a C4 block or not
//! \return bool
bool layer::is_c4_layer()
{
    return _is_c4_layer;
}

//! This method returns a flag that indicates whether the layer contains current sources or not
//! \return bool
bool layer::is_current_source_layer()
{
    return _is_current_source_layer;
}

//! This method returns a flag that indicates the orientation of the layer
//! \return bool
bool layer::get_orientation()
{
    return _specs->orientation;
}

// Setter Functions

//! This method sets the name of the layer
//! @param name [in]: name of the layer
void layer::set_name(string name)
{
    _name = name;
}

//! This method sets the index of the layer
//! @param index [in]: index of the layer
void layer::set_index(int index)
{
    _index = index;
}

//! This method adds an interconnect_tree to the list of interconnect_trees that are in the layer
//! @param it [in]: interconnect_tree to add
void layer::add_interconnect_tree(interconnect_tree * it)
{
    _list_of_interconnect_trees.push_back(it);
}

//! This method adds an interconnect tree to the interconnect trees map of the layer associated with a coordinate 
//! @param coord [in]: coordinate associated with the interconnect_tree in the map
//! @param it    [in]: interconnect_tree to add
void layer::add_it_to_map(double coord, interconnect_tree* it)
{
    //_it_map[coord] = it;
    _it_map.insert(pair<double, interconnect_tree*>(coord, it));
}

// Other Functions

//! This method prints the data of the layer
void layer::print()
{
	cout<<"Layer " << _name << endl;
	cout<<"Number of Interconnect Trees: " << _list_of_interconnect_trees.size() << endl;
    cout<<"Orientation and Type: "<<
     ( (get_orientation() == HORIZONTAL )?"Horizontal, ":"Vertical, ");
    if (_is_current_source_layer)
        cout<<" Current Source Layer"<<endl;
    else if (_is_c4_layer)
        cout<<" C4 Layer"<<endl;
    cout << endl;
}
