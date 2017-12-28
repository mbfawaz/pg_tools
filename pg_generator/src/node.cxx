/**
* \file node.cxx \brief node Implementation
*
* This file implements the methods of the class node
*/

#include "power_grid.h"

// Node Class

// ! Default Constructor
node::node()
{
    _interconnect_tree_ptr = NULL;
    _name                  = "";
    _index                 = -1;
    _xcoord                = -1;
    _ycoord                = -1;
    _type                  = CHIP_NODE;
    _c4_pad_connection     = NULL;
    _c4_pad_ptr            = NULL;
    _cs_connection         = NULL;
}

// ! Regular Constructor
// ! @param it       [in]: pointer to the interconnect that the node is in
// ! @param name     [in]: name of the node
// ! @param index    [in]: index of node
// ! @param xcoord   [in]: x_coordinate of the node
// ! @param ycoord   [in]: y_coordinate of the node
// ! @param type     [in]: type of the node
node::node(interconnect_tree * it,
           string name,
           int    index,
           double xcoord,
           double ycoord,
           bool   type)
{
    _interconnect_tree_ptr = it;
    _name                  = name;
    _index                 = index;
    _xcoord                = xcoord;
    _ycoord                = ycoord;
    _type                  = type;
    _c4_pad_connection     = NULL;
    _c4_pad_ptr            = NULL;
    _cs_connection         = NULL;
}

// ! Destructor of class node
node::~node()
{
    // Deleting C connections, only if this is the ground node. this check is need
    // so that a capacitor is not deleted twice (once when deleting the ground node
    // and once when deleting the power grid node)
    if ( _name.compare("n0") == 0 )
       for (unsigned int i = 0; i < _c_connections.size(); i++)
            delete _c_connections[i];

    // Deleting the current source attached
    if (_cs_connection != NULL)
        delete _cs_connection;
}

// Getter Functions

//! This method returns a pointer to the interconnect_tree that the node is in
//! \return interconnect_tree*
interconnect_tree* node::get_interconnect_tree_ptr()
{
    return _interconnect_tree_ptr;
}

//! This method returns a pointer to the C4 pad that the node is connected to
//! \return c4_pad*
c4_pad* node::get_c4_pad_connection()
{
    return _c4_pad_connection;
}

//! This method returns a pointer to the C4 pad that the node is in
//! \return c4_pad*
c4_pad* node::get_c4_pad_ptr()
{
    return _c4_pad_ptr;
}

//! This method returns a copy of the list of resistors attached to the node
//! \return vector<resistor*>
vector<resistor*> node::get_r_connections()
{
    return _r_connections;
}

//! This method returns a reference copy of the list of resistors attached to the node
//! \return vector<resistor*> &
vector<resistor*>& node::get_r_connections_ref()
{
    return _r_connections;
}

//! This method returns a copy of the list of capacitors attached to the node
//! \return vector<capacitor*>
vector<capacitor*> node::get_c_connections()
{
    return _c_connections;
}

//! This method returns a reference copy of the list of capacitors attached to the node
//! \return vector<capacitor*> &
vector<capacitor*>& node::get_c_connections_ref()
{
    return _c_connections;
}

//! This method returns a copy of the list of inductors attached to the node
//! \return vector<inductor*>
vector<inductor*> node::get_l_connections()
{
    return _l_connections;
}

//! This method returns a reference copy of the list of inductors attached to the node
//! \return vector<inductor*> &
vector<inductor*>& node::get_l_connections_ref()
{
    return _l_connections;
}

//! This method returns a pointer to the current source attached to the node
//! \return current_source*
current_source* node::get_cs_connection()
{
   return _cs_connection;
}

//! This method returns the name of the node
//! \return string
string node::get_name()
{
    return _name;
}

//! This method returns the index of the node
//! \return int
int node::get_index()
{
    return _index;
}

//! This method returns the x_coordinate of the node
//! \return double
double node::get_xcoord()
{
    return _xcoord;
}

//! This method returns the y_coordinate of the node
//! \return double
double node::get_ycoord()
{
    return _ycoord;
}

//! This method returns the type of the node
//! \return bool
bool node::get_type()
{
    return _type;
}

// Setter Functions

//! This method sets the name of the node
//! @param name [in]: name of the node
void node::set_name(string name)
{
    _name = name;
}

//! This method sets the index of the node
//! @param index [in]: index of the node
void node::set_index(int index)
{
    _index = index;
}

//! This method adds a resistor to the list of resistors attached to the node
//! @param r [in]: resistor to add
void node::add_resistor(resistor * r)
{
    _r_connections.push_back(r);
}

//! This method adds a capacitor to the list of capacitors attached to the node
//! @param c [in]: capacitor to add
void node::add_capacitor(capacitor * c)
{
    _c_connections.push_back(c);
}

//! This method adds an inductor to the list of inductors attached to the node
//! @param l [in]: inductor to add
void node::add_inductor(inductor * l)
{
    _l_connections.push_back(l);
}

//! This method sets the current source attached to the node
//! @param cs [in]: current source attached to the node
void node::add_current_source(current_source * cs)
{
    _cs_connection = cs;
}

//! This method sets the C4 pad that the node is connected to
//! @param c4 [in]: C4 pad the node is connected to 
void node::add_c4_pad_connection(c4_pad * c4)
{
    _c4_pad_connection = c4;
}

//! This method sets the C4 pad that the node is in
//! @param c4 [in]: C4 pad the node is in 
void node::add_c4_pad_ptr(c4_pad * c4)
{
    _c4_pad_ptr = c4;
}

// Other Functions

//! This method returns the list of nodes that have connections with the this node
//! \return vector<node *>
vector<node *> node::get_list_of_neighboring_nodes()
{
    vector<node *> neighbors(_r_connections.size());
    for (unsigned int i = 0; i < _r_connections.size(); i++)
    {
        if (_r_connections.at(i)->get_n0()->get_index() == _index)
            neighbors[i] = _r_connections.at(i)->get_n1();
        else
            neighbors[i] = _r_connections.at(i)->get_n0();
    }
    return neighbors;
}

//! This method is to print the data of the node
void node::print()
{
    cout << "Node " << _name<< " at (" << _xcoord << "," << _ycoord << 
        ") with aliases is connected to:" << endl;
    for (unsigned int i = 0; i<_r_connections.size(); i++)
        _r_connections[i]->print();
     for (unsigned int i = 0; i<_c_connections.size(); i++)
        _c_connections[i]->print();
     for (unsigned int i = 0; i<_l_connections.size(); i++)
        _l_connections[i]->print();
}
