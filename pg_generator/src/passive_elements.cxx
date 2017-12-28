/**
* \file passive_elements.cxx \brief Resistor, Capacitor and Inductor Implementation
*
* This file implements the methods of the classes resistor, capacitor and inductor.
*/

#include "power_grid.h"

// Resistor Class

// ! Default Constructor
resistor::resistor()
{
    _name                 = "";
    _index                = -1;
    _value                = -1;
    _length               = -1;
    _cross_sectional_area = -1;
    _type                 = BRANCH;
    _n0                   = NULL;
    _n1                   = NULL;
}

// ! Regular Constructor
// ! @param name    [in]: name of the resistor
// ! @param index   [in]: index of the resistor
// ! @param value   [in]: value of the resistor
// ! @param length  [in]: length of the resistor
// ! @param cs_area [in]: cross sectional area of the resistor
// ! @param type    [in]: type of the resistor
resistor::resistor(string name,
                   int    index,
                   double value,
                   double length,
                   double cs_area,
                   int    type)
{
    _name                 = name;
    _index                = index;
    _value                = value;
    _length               = length;
    _cross_sectional_area = cs_area;
    _type                 = type;
}

// Destructor of class resistor
resistor::~resistor()
{
}

// Getter Functions

//! This method returns a pointer to the node at terminal 0 of the resistor
//! \return node*
node * resistor::get_n0()
{
    return _n0;
}

//! This method returns a pointer to the node at terminal 1 of the resistor
//! \return node*
node * resistor::get_n1()
{
    return _n1;
}

//! This method returns the name of the resistor
//! \return string
string resistor::get_name()
{
    return _name;
}

//! This method returns the index of the resistor
//! \return int
int resistor::get_index()
{
    return _index;
}

//! This method returns the value of the resistor
//! \return double
double resistor::get_value()
{
    return _value;
}

//! This method returns the length of the resistor
//! \return double
double resistor::get_length()
{
    return _length;
}

//! This method returns the cross sectional area of the resistor
//! \return double
double resistor::get_cross_sectional_area()
{
    return _cross_sectional_area;
}

//! This method returns the type of the resistor
//! \return int
int resistor::get_type()
{
    return _type;
}

// Setter Functions

//! This method sets the name of the resistor
//! @param name [in]: name of the resistor
void resistor::set_name(string name)
{
    _name = name;
}

//! This method sets the index of the resistor
//! @param index [in]: index of the resistor
void resistor::set_index(int index)
{
    _index = index;
}

//! This method sets the value of the resistor
//! @param v [in]: value of the resistor
void resistor::set_value(double v)
{
    _value = v;
}

//! This method sets the length of the resistor
//! @param l [in]: length of the resistor
void resistor::set_length(double l)
{
    _length = l;
}

//! This method sets the cross sectional area of the resistor
//! @param cs_area [in]: cross sectional area of the resistor
void resistor::set_cross_sectional_area(double cs_area)
{
    _cross_sectional_area = cs_area;
}

//! This method sets the type of the resistor
//! @param type [in]: type of the resistor
void resistor::set_type(int type)
{
    _type = type;
}

//! This method sets the terminal nodes of the resistor according the convention that the node with lower index will be in n0
//! @param n0 [in]: first terminal node of the resistor
//! @param n1 [in]: second terminal node of the resistor
void resistor::set_nodes(node* n0, node* n1)
{
    if (n0->get_index() == n1->get_index())
    {
        std::cerr << "Cannot create a resistor with same terminal nodes!" << std::endl;
        exit(-1);
    }
    else if (n0->get_index() < n1->get_index())
    {
        _n0 = n0;
        _n1 = n1;
    }
    else
    {
        _n0 = n1;
        _n1 = n0;
    }
}

// Other Functions

//! This method prints the data of the resistor
void resistor::print()
{
    if (_n0->get_index() == 0)
		cout<<"Resistor " << _name << " is connected between "
		    << _n1->get_name()<<" and ground and has value = "<<_value<<endl;
	else
		cout<<"Resistor "<<_name<<" is connected between "
		    <<_n0->get_name()<<" and "  << _n1->get_name() << " and has value = "
		    <<_value<<endl;
}

// Capacitor Class

// ! Default Constructor
capacitor::capacitor()
{
    _name                 = "";
    _index                = -1;
    _value                = -1;
    _type                 = BRANCH;
    _n0                   = NULL;
    _n1                   = NULL;
}

// ! Regular Constructor
// ! @param name    [in]: name of the capacitor
// ! @param index   [in]: index of the capacitor
// ! @param value   [in]: value of the capacitor
// ! @param type    [in]: type of the capacitor
capacitor::capacitor(string name,
                   int    index,
                   double value,
                   int    type)
{
    _name                 = name;
    _index                = index;
    _value                = value;
    _type                 = type;
}

// Destructor of class capacitor
capacitor::~capacitor()
{
}

// Getter Functions

//! This method returns a pointer to the node at terminal 0 of the capacitor
//! \return node*
node * capacitor::get_n0()
{
    return _n0;
}

//! This method returns a pointer to the node at terminal 1 of the capacitor
//! \return node*
node * capacitor::get_n1()
{
    return _n1;
}

//! This method returns the name of the capacitor
//! \return string
string capacitor::get_name()
{
    return _name;
}

//! This method returns the index of the capacitor
//! \return int
int capacitor::get_index()
{
    return _index;
}

//! This method returns the value of the capacitor
//! \return double
double capacitor::get_value()
{
    return _value;
}

//! This method returns the type of the capacitor
//! \return int
int capacitor::get_type()
{
    return _type;
}

// Setter Functions

//! This method sets the name of the capacitor
//! @param name [in]: name of the capacitor
void capacitor::set_name(string name)
{
    _name = name;
}

//! This method sets the index of the capacitor
//! @param index [in]: index of the capacitor
void capacitor::set_index(int index)
{
    _index = index;
}

//! This method sets the value of the capacitor
//! @param v [in]: value of the capacitor
void capacitor::set_value(double v)
{
    _value = v;
}

//! This method sets the type of the capacitor
//! @param type [in]: type of the capacitor
void capacitor::set_type(int type)
{
    _type = type;
}

//! This method sets the terminal nodes of the capacitor according the convention that the node with lower index will be in n0
//! @param n0 [in]: first terminal node of the capacitor
//! @param n1 [in]: second terminal node of the capacitor
void capacitor::set_nodes(node* n0, node* n1)
{
    if (n1 != NULL)
    {
        if (n0->get_index() == n1->get_index())
        {
            std::cerr << "Cannot create a capacitor with same terminal nodes!" << std::endl;
            exit(-1);
        }
        else if (n0->get_index() < n1->get_index())
        {
            _n0 = n0;
            _n1 = n1;
        }
        else
        {
            _n0 = n1;
            _n1 = n0;
        }
    }
    else
    {
        _n0 = n0;
        _n1 = NULL;
    }
}

// Other Functions

//! This method prints the data of the capacitor
void capacitor::print()
{
    if (_n0->get_index() == 0)
		cout<<"Capacitor " << _name << " is connected between "
		    << _n1->get_name()<<" and ground and has value = "<<_value<<endl;
	else
		cout<<"Capacitor "<<_name<<" is connected between "
		    <<_n0->get_name()<<" and " <<  _n1->get_name() <<  " and has value = "
		    <<_value<<endl;
}

// Inductor Class

// ! Default Constructor
inductor::inductor()
{
    _name                 = "";
    _index                = -1;
    _value                = -1;
    _type                 = BRANCH;
    _n0                   = NULL;
    _n1                   = NULL;
}

// ! Regular Constructor
// ! @param name    [in]: name of the inductor
// ! @param index   [in]: index of the inductor
// ! @param value   [in]: value of the inductor
// ! @param type    [in]: type of the inductor
inductor::inductor(string name,
                   int    index,
                   double value,
                   int    type)
{
    _name                 = name;
    _index                = index;
    _value                = value;
    _type                 = type;
}

// Destructor of class inductor
inductor::~inductor()
{
}

// Getter Functions

//! This method returns a pointer to the node at terminal 0 of the inductor
//! \return node*
node * inductor::get_n0()
{
    return _n0;
}

//! This method returns a pointer to the node at terminal 1 of the inductor
//! \return node*
node * inductor::get_n1()
{
    return _n1;
}

//! This method returns the name of the inductor
//! \return string
string inductor::get_name()
{
    return _name;
}

//! This method returns the index of the inductor
//! \return int
int inductor::get_index()
{
    return _index;
}

//! This method returns the value of the inductor
//! \return double
double inductor::get_value()
{
    return _value;
}

//! This method returns the type of the inductor
//! \return int
int inductor::get_type()
{
    return _type;
}

// Setter Functions

//! This method sets the name of the inductor
//! @param name [in]: name of the inductor
void inductor::set_name(string name)
{
    _name = name;
}

//! This method sets the index of the inductor
//! @param index [in]: index of the inductor
void inductor::set_index(int index)
{
    _index = index;
}

//! This method sets the value of the inductor
//! @param v [in]: value of the inductor
void inductor::set_value(double v)
{
    _value = v;
}

//! This method sets the type of the inductor
//! @param type [in]: type of the inductor
void inductor::set_type(int type)
{
    _type = type;
}

//! This method sets the terminal nodes of the inductor according the convention that the node with lower index will be in n0
//! @param n0 [in]: first terminal node of the inductor
//! @param n1 [in]: second terminal node of the inductor
void inductor::set_nodes(node* n0, node* n1)
{
    if (n0->get_index() == n1->get_index())
    {
        std::cerr << "Cannot create a inductor with same terminal nodes!" << std::endl;
        exit(-1);
    }
    else if (n0->get_index() < n1->get_index())
    {
        _n0 = n0;
        _n1 = n1;
    }
    else
    {
        _n0 = n1;
        _n1 = n0;
    }
}

// Other Functions

//! This method prints the data of the inductor
void inductor::print()
{
    if (_n0->get_index() == 0)
		cout<<"Inductor " << _name << " is connected between "
		    << _n1->get_name()<<" and ground and has value = "<<_value<<endl;
	else
		cout<<"Inductor "<<_name<<" is connected between "
		    <<_n0->get_name()<<" and " << _n1->get_name() << " and has value = "
		    <<_value<<endl;
}
