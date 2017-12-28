/**
* \file interconnect_tree.cxx \brief Interconnect Tree Implementation
*
* This file implements the methods of the class interconnect_tree
*/

#include "power_grid.h"

// Interconnect_tree Class

// ! Default Constructor
interconnect_tree::interconnect_tree()
{
    _layer_ptr = NULL;
    _name      = "";
    _index     = -1;
    _xcoord    = -1;
    _ycoord    = -1;
}

// ! Regular Constructor
// ! @param layer_ptr [in]: pointer to the layer that the interconnect tree is in
// ! @param name      [in]: name of the interconnect_tree
// ! @param index     [in]: index of interconnect_tree
// ! @param xcoord    [in]: x_coordinate of the interconnect_tree
// ! @param ycoord    [in]: y_coordinate of the interconnect_tree
interconnect_tree::interconnect_tree(layer* layer_ptr,
                                     string name,
                                     int    index,
                                     double xcoord,
                                     double ycoord)
{
    _layer_ptr = layer_ptr;
    _name      = name;
    _index     = index;
    _xcoord    = xcoord;
    _ycoord    = ycoord;
}

// ! Destructor of class interconnect_tree
interconnect_tree::~interconnect_tree()
{
    // Deleting nodes
    for (unsigned int i = 0; i < _list_of_nodes.size(); i++)
        delete _list_of_nodes[i];
    
    // Deleting Resistors
    for (unsigned int i = 0; i < _list_of_resistors.size(); i++)
        delete _list_of_resistors[i];
}
// Getter Functions

//! This method returns a pointer to the layer that the interconnect_tree is in
//! \return layer*
layer* interconnect_tree::get_layer_ptr()
{
    return _layer_ptr;
}

//! This method returns a copy of the list of nodes that are in the interconnect_tree
//! \return vector<node*>
vector<node*> interconnect_tree::get_list_of_nodes()
{
    return _list_of_nodes;
}

//! This method returns a reference copy of the list of nodes that are in the interconnect_tree
//! \return vector<node*> &
vector<node*>& interconnect_tree::get_list_of_nodes_ref()
{
    return _list_of_nodes;
}

//! This method returns a copy of the list of resistors that are in the interconnect_tree
//! \return vector<resistor*>
vector<resistor*> interconnect_tree::get_list_of_resistors()
{
    return _list_of_resistors;
}

//! This method returns a reference copy of the list of resistors that are in the interconnect_tree
//! \return vector<resistor*> &
vector<resistor*>& interconnect_tree::get_list_of_resistors_ref()
{
    return _list_of_resistors;
}

//! This method returns a copy of the nodes map in the interconnect tree 
//! \return map<double, node*>
map<double, node*> interconnect_tree::get_nodes_map()
{
    return _nodes_map;
}

//! This method returns the name of the interconnect_tree
//! \return string
string interconnect_tree::get_name()
{
    return _name;
}

//! This method returns the index of the interconnect_tree
//! \return int
int interconnect_tree::get_index()
{
    return _index;
}

//! This method returns the x_coordinate of the interconnect_tree
//! \return double
double interconnect_tree::get_xcoord()
{
    return _xcoord;
}

//! This method returns the y_coordinate of the interconnect_tree
//! \return double
double interconnect_tree::get_ycoord()
{
    return _ycoord;
}

// Setter Functions

//! This method sets the name of the interconnect_tree
//! @param name [in]: name of the interconnect_tree
void interconnect_tree::set_name(string name)
{
    _name = name;
}

//! This method sets the index of the interconnect_tree
//! @param index [in]: index of the interconnect_tree
void interconnect_tree::set_index(int index)
{
    _index = index;
}

//! This method adds a resistor to the list of resistors that are in the interconnect_tree
//! @param r [in]: resistor to add
void interconnect_tree::add_resistor(resistor * r)
{
    _list_of_resistors.push_back(r);
}

//! This method adds a node to the list of nodes that are in the interconnect_tree
//! @param n [in]: node to add
void interconnect_tree::add_node(node * n)
{
    _list_of_nodes.push_back(n);
}

//! This method adds a node to the nodes map of the interconnect tree associated with a coordinate 
//! @param coord [in]: coordinate associated with the node in the map
//! @param n     [in]: node to add
void interconnect_tree::add_node_to_map(double coord, node * n)
{
    //_nodes_map[coord] = n;
    _nodes_map.insert(pair<double, node*>(coord, n));
}

// Other Functions

//! This method checks if a node, with the given coordinates, already exists in the interconnect tree or not
//! @param xcoord [in]: x_coordinate of the node to find
//! @param ycoord [in]: y_coordinate of the node to find
//! \return bool
bool interconnect_tree::find_node(double xcoord, double ycoord)
{
    for (unsigned int i = 0; i < _list_of_nodes.size(); i++)
    {
        if (_list_of_nodes.at(i)->get_xcoord() == xcoord && _list_of_nodes.at(i)->get_ycoord() == ycoord)
        {
            return true;
        }
    }
    return false;
}

//! This method prints the data of the interconnect tree
void interconnect_tree::print()
{
    cout << "Interconnect Tree: " << _name << " in Layer: " << _layer_ptr->get_name() << endl;
    cout << "Number of Nodes: " << _list_of_nodes.size() << endl;
    cout << "Number of Resistors: " << _list_of_resistors.size() << endl;
    cout << endl;
}
