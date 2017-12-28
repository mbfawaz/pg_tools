/**
* \file current_source.cxx \brief current_source Implementation
*
* This file implements the methods of the class current_source
*/

#include "power_grid.h"

// Current Source Class

// ! Default Constructor
current_source::current_source()
{
    _name       = "";
    _index      = -1;
    _n0         = NULL;
    _layer_ptr  = NULL;
    _value      = 0;
}

// ! Regular Constructor
// ! @param layer    [in]: layer that the current_source is in
// ! @param name     [in]: name of the current_source
// ! @param index    [in]: index of the current_source
// ! @param index    [in]: value of the current_source
current_source::current_source(layer* layer,
                               string name,
                               int    index,
                               double value)
{
    _layer_ptr = layer;
    _name      = name;
    _index     = index;
    _value     = value;
}

// ! Destructor of class current_source
current_source::~current_source()
{
}

// Getter Functions

//! This method returns a pointer to the layer that the current_source is in
//! \return layer*
layer* current_source::get_layer_ptr()
{
    return _layer_ptr;
}

//! This method returns a pointer to the node that the current_source is connected to
//! \return node*
node* current_source::get_n0()
{
    return _n0;
}

//! This method returns the name of the current_source
//! \return string
string current_source::get_name()
{
    return _name;
}

//! This method returns the index of the current_source
//! \return int
int current_source::get_index()
{
    return _index;
}

//! This method returns the value of the current_source
//! \return value
double current_source::get_value()
{
    return _value;
}

// Setter Functions

//! This method sets the name of the current_source
//! @param name [in]: name of the current_source
void current_source::set_name(string name)
{
    _name = name;
}

//! This method sets the index of the current_source
//! @param index [in]: index of the current_source
void current_source::set_index(int index)
{
    _index = index;
}

//! This method sets the node that the current_source is connected to
//! @param n0 [in]: node that the current_source is connected to
void current_source::set_node(node* n0)
{
    _n0 = n0;
}

//! This method sets the value (in Amps)for the current_source
//! @param n0 [in]: double value fo rthe current source
void current_source::set_value(double value)
{
    _value = value;
}

// Other Functions 

//! This method prints the data of the current source
void current_source::print()
{
    cout << "Current source: " << _name << " connected to: " << _n0->get_name() << endl;
}
