/**
* \file c4_pad.cxx \brief C4_Pad Implementation
*
* This file implements the methods of the class c4_pad
*/

#include "power_grid.h"
#include <string.h>
using namespace std;

// C4 Pad Class

// ! Default Constructor
c4_pad::c4_pad()
{
    _name         = "";
    _index        = -1;
    _port_node    = NULL;
    _power_grid   = NULL;
}

// ! Regular Constructor
// ! @param pg           [in]: pointer to the power grid that the c4 is connected to
// ! @param port_node    [in]: pointer to the port node that the c4 pad is connected to
// ! @param name         [in]: name of the c4 pad
// ! @param index        [in]: index of the c4 pad
c4_pad::c4_pad(power_grid* pg,
               node*       port_node, 
               string      name,
               int         index)
{
    _power_grid    = pg;
    _port_node     = port_node;
    _name          = name;
    _index         = index;
}

// ! Destructor of class c4_pad
c4_pad::~c4_pad()
{
    // Deleting resistors
    for (unsigned int i = 0; i < _list_of_resistors.size(); i++)
        delete _list_of_resistors[i];

    // Deleting nodes
    for (unsigned int i = 0; i < _list_of_nodes.size(); i++)
        delete _list_of_nodes[i];
    
    // Deleting inductors
    for (unsigned int i = 0; i < _list_of_inductors.size(); i++)
        delete _list_of_inductors[i];
}

// Getter Functions

//! This method returns a pointer to the power grid that the c4_pad is connected to
//! \return power_grid*
power_grid* c4_pad::get_power_grid_ptr()
{
    return _power_grid;
}

//! This method returns a pointer to the node that the c4_pad is connected to
//! \return node*
node* c4_pad::get_port_node_ptr()
{
    return _port_node;
}

//! This method returns a copy of the list of nodes that are in the c4 pad
//! \return vector<node*>
vector<node*> c4_pad::get_list_of_nodes()
{
    return _list_of_nodes;
}

//! This method returns a reference copy of the list of nodes that are in the c4 pad
//! \return vector<node*> &
vector<node*>& c4_pad::get_list_of_nodes_ref()
{
    return _list_of_nodes;
}

//! This method returns a copy of the list of resistors that are in the c4 pad
//! \return vector<resistor*>
vector<resistor*> c4_pad::get_list_of_resistors()
{
    return _list_of_resistors;
}

//! This method returns a reference copy of the list of resistors that are in the c4 pad
//! \return vector<resistor*> &
vector<resistor*>& c4_pad::get_list_of_resistors_ref()
{
    return _list_of_resistors;
}

//! This method returns a copy of the list of capacitors that are in the c4 pad
//! \return vector<capacitor*>
vector<capacitor*> c4_pad::get_list_of_capacitors()
{
    return _list_of_capacitors;
}

//! This method returns a reference copy of the list of capacitors that are in the c4 pad
//! \return vector<capacitor*> &
vector<capacitor*>& c4_pad::get_list_of_capacitors_ref()
{
    return _list_of_capacitors;
}

//! This method returns a copy of the list of inductors that are in the c4 pad
//! \return vector<intductor*>
vector<inductor*> c4_pad::get_list_of_inductors()
{
    return _list_of_inductors;
}

//! This method returns a reference copy of the list of inductors that are in the c4 pad
//! \return vector<inductor*> &
vector<inductor*>& c4_pad::get_list_of_inductors_ref()
{
    return _list_of_inductors;
}

//! This method returns the name of the c4 pad
//! \return string
string c4_pad::get_name()
{
    return _name;
}

//! This method returns the index of the c4 pad
//! \return int
int c4_pad::get_index()
{
    return _index;
}

// Setter Functions

//! This method sets the name of the c4 pad
//! @param name [in]: name of the c4 pad
void c4_pad::set_name(string name)
{
    _name = name;
}

//! This method sets the index of the c4 pad
//! @param index [in]: index of the c4 pad
void c4_pad::set_index(int index)
{
    _index = index;
}

//! This method checks if the node to add already exists. If not, it adds a
//! node to the list of nodes that are in the c4 pad and returns it
//! @param name [in]: name of the node to add
//! @param idx  [in, out]: index of the last node added
node* c4_pad::add_node(string name, unsigned int & idx)
{
    stringstream ss;
    ss << name << "_n" << _name;
    string node_name = ss.str();
    if (strcmp(name.c_str(), "0") == 0)
        return _power_grid->get_ground_node();
    else
        for (unsigned int i = 0; i < _list_of_nodes.size(); i++)
            if (strcmp(_list_of_nodes[i]->get_name().c_str(), node_name.c_str()) == 0)
                return _list_of_nodes[i];
    idx++;
    node* new_node;
    new_node = new node(NULL, node_name, idx, -1, -1, PACKAGE_NODE);
    new_node->add_c4_pad_ptr(this);
    _list_of_nodes.push_back(new_node);
    return new_node;
}

//! This method adds a resistor to the list of resistors that are in the c4 pad
//! @param name     [in]: name of resistor to add
//! @param node1    [in]: first terminal of the resistor
//! @param node2    [in]: second terminal of the resistor
//! @param value    [in]: value of the resistor
//! @param idx      [in, out]: index of the last resistor added
//! @param node_idx [in, out]: index of the last node added
void c4_pad::add_resistor(string name, string node1, string node2, double value,
    int & idx, unsigned int & node_idx)
{
    resistor* new_resistor;
    stringstream ss;
    transform(name.begin(), name.end(), name.begin(), ::tolower);
    ss << name << "_" << _name;
    string res_name = ss.str();
    idx++;
    new_resistor = new resistor(res_name, idx, value, -1, -1, PACKAGE);
    
    stringstream ss1, ss2;
    ss1 << _name << "_n" << node1;
    ss2 << _name << "_n" << node2;
    string n1, n2;
    n1 = ss1.str();
    n2 = ss2.str();

    node *node_1, *node_2;
    if (strcmp(node1.c_str(), "n_prt") == 0)
    {
        node_1 = _port_node;
    }
    else if (strcmp(node1.c_str(), "n_vdd") == 0)
    {
        node_1 = _power_grid->get_vdd_node();
    }
    else if (strcmp(node1.c_str(), "0") == 0)
    {
        node_1 = _power_grid->get_ground_node();
    }
    else
    {
        node_1 = add_node(node1, node_idx);
    }

    if (strcmp(node2.c_str(), "n_prt") == 0)
    {
        node_2 = _port_node;
    }
    else if (strcmp(node2.c_str(), "n_vdd") == 0)
    {
        node_2 = _power_grid->get_vdd_node();
    }
    else if (strcmp(node2.c_str(), "0") == 0)
    {
        node_2 = _power_grid->get_ground_node();
    }
    else
    {
        node_2 = add_node(node2, node_idx);
    }
    new_resistor->set_nodes(node_1, node_2);
    node_1->add_resistor(new_resistor);
    node_2->add_resistor(new_resistor);
    _list_of_resistors.push_back(new_resistor);
}

//! This method adds a capacitor to the list of capacitors that are in the c4 pad
//! @param name     [in]: name of capacitor to add
//! @param node1    [in]: first terminal of the capacitor
//! @param node2    [in]: second terminal of the capacitor
//! @param value    [in]: value of the capacitor
//! @param idx      [in, out]: index of the last capacitor added
//! @param node_idx [in, out]: index of the last node added
void c4_pad::add_capacitor(string name, string node1, string node2, double value, int & idx, unsigned int & node_idx)
{
    capacitor* new_capacitance;
    stringstream ss;
    transform(name.begin(), name.end(), name.begin(), ::tolower);
    ss << name << "_" << _name;
    string cap_name = ss.str();
    idx++;
    new_capacitance = new capacitor(cap_name, idx, value, PACKAGE);

    stringstream ss1, ss2;
    ss1 << _name << "_n" << node1;
    ss2 << _name << "_n" << node2;
    string n1, n2;
    n1 = ss1.str();
    n2 = ss2.str();
    
    node *node_1, *node_2;
    if (strcmp(node1.c_str(), "n_prt") == 0 || strcmp(node2.c_str(), "n_prt") == 0)
    {
        cerr << "This version does not support branch capacitance..." << endl;
        exit(1);
    }
    if (strcmp(node1.c_str(), "0") == 0) 
    {
        node_1 = _power_grid->get_ground_node();
        node_2 = add_node(node2, node_idx);
        
    }
    else if (strcmp(node2.c_str(), "0") == 0)
    {
        node_2 = _power_grid->get_ground_node();
        node_1 = add_node(node1, node_idx);
    }
    else
    {
        cerr<<"This version does not support branch capacitance..."<<endl;
        exit(1);
    }
    
    new_capacitance->set_nodes(node_1, node_2);
    node_1->add_capacitor(new_capacitance);
    node_2->add_capacitor(new_capacitance);
    _list_of_capacitors.push_back(new_capacitance);
}

//! This method adds an inductor to the list of inductors that are in the c4 pad
//! @param name     [in]: name of inductor to add
//! @param node1    [in]: first terminal of the inductor
//! @param node2    [in]: second terminal of the inductor
//! @param value    [in]: value of the inductor
//! @param idx      [in, out]: index of the last inductor added
//! @param node_idx [in, out]: index of the last node added
void c4_pad::add_inductor(string name, string node1, string node2, double value, int & idx, unsigned int & node_idx)
{
    inductor* new_inductor;
    stringstream ss;
    transform(name.begin(), name.end(), name.begin(), ::tolower);
    ss << name << "_" << _name;
    string ind_name = ss.str();
    idx++;
    new_inductor = new inductor(ind_name, idx, value, PACKAGE);
    
    stringstream ss1, ss2;
    ss1 << _name << "_n" << node1;
    ss2 << _name << "_n" << node2;
    string n1, n2;
    n1 = ss1.str();
    n2 = ss2.str();

    node *node_1, *node_2;
    if (strcmp(node1.c_str(), "n_prt") == 0)
    {
        node_1 = _port_node;
    }
    else if (strcmp(node1.c_str(), "n_vdd") == 0)
    {
        node_1 = _power_grid->get_vdd_node();
    }
    else if (strcmp(node1.c_str(), "0") == 0)
    {
        node_1 = _power_grid->get_ground_node();
    }
    else
    {
        node_1 = add_node(node1, node_idx);
    }

    if (strcmp(node2.c_str(), "n_prt") == 0)
    {
        node_2 = _port_node;
    }
    else if (strcmp(node2.c_str(), "n_vdd") == 0)
    {
        node_2 = _power_grid->get_vdd_node();
    }
    else if (strcmp(node2.c_str(), "0") == 0)
    {
        node_2 = _power_grid->get_ground_node();
    }
    else
    {
        node_2 = add_node(node2, node_idx);
    }

    new_inductor->set_nodes(node_1, node_2);
    node_1->add_inductor(new_inductor);
    node_2->add_inductor(new_inductor);
    _list_of_inductors.push_back(new_inductor);
}

// Other Functions

//! This method reads a spice file and creates the list of passive elements in the c4 pad
//! @param node_idx [in, out]: index of the last node added to the grid
//! @param res_idx  [in, out]: index of the last resistor added to the grid
//! @param cap_idx  [in, out]: index of the last capacitor added to the grid
//! @param ind_idx  [in, out]: index of the last inductor added to the grid
void c4_pad::parse_subcircuit_file(unsigned int & node_idx, int & res_idx, int & cap_idx, int & ind_idx)
{
//     string sub_ckt_file = "..//input//c4_pad.sp";
    string sub_ckt_file = "c4_pad.sp";

    ifstream infile;
    infile.open(sub_ckt_file.c_str());
    if (!infile)
    {
        cerr<<"Could not open C4 pad subcircuit file"<<endl;
        exit(0);
    }
    else
    {
        while (!infile.eof())
        {
            char str[255];
            infile.getline(str, 255);
            char *param = strtok(str, " ");
            if (param[0] == 'R')
            {
                string node1 = strtok(NULL, " ");
                string node2 = strtok(NULL, " ");
                string str_value = strtok(NULL, " ");
                add_resistor(param, node1, node2, atof(str_value.c_str()), res_idx, node_idx);
            }
            else if (param[0] == 'C')
            {   
                string node1 = strtok(NULL, " ");
                string node2 = strtok(NULL, " ");
                string str_value = strtok(NULL, " ");
                add_capacitor(param, node1, node2, atof(str_value.c_str()), cap_idx, node_idx);
            }
            else if (param[0] == 'L')
            {   
                string node1 = strtok(NULL, " ");
                string node2 = strtok(NULL, " ");
                string str_value = strtok(NULL, " ");
                add_inductor(param, node1, node2, atof(str_value.c_str()), ind_idx, node_idx);
            }
            else if (strcmp(param, ".ends") == 0)
                break; 
        }
    }
}

//! This method is to print the data of the c4 pad
void c4_pad::print()
{
    cout << "C4: " << _name << endl;
    cout<<"List of Resistors: "<<endl;
    cout<<"---------------"<<endl;
    for (unsigned int i = 0; i<_list_of_resistors.size(); i++)
        _list_of_resistors[i]->print();
    cout<<endl;
    cout<<"List of Capacitors: "<<endl;
    cout<<"---------------"<<endl;
    for (unsigned int i = 0; i<_list_of_capacitors.size(); i++)
        _list_of_capacitors[i]->print();
    cout<<endl;
    cout<<"List of Inductors: "<<endl;
    cout<<"---------------"<<endl;
    for (unsigned int i = 0; i<_list_of_inductors.size(); i++)
        _list_of_inductors[i]->print();
    cout<<endl;
    cout<<"List of Nodes: "<<endl;
    cout<<"---------------"<<endl;
    for (unsigned int i = 0; i<_list_of_nodes.size(); i++)
    {    
        _list_of_nodes[i]->print();
        cout<<endl;
    }
}
