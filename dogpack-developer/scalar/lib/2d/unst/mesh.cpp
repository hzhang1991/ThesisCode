// --------------------------------------------------------------------------
//  IMPLEMENTATION FILE (mesh.cpp)
//    2d unstructured mesh
// --------------------------------------------------------------------------

#include "mesh.h"
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <iomanip>
#include <new>
#include <sstream>
#include <string>
using namespace std;


mesh::mesh(int inNumElems, 
        int inNumPhysElems, 
        int inNumNodes, 
        int inNumPhysNodes, 
        int inNumBndNodes, 
        int inNumEdges,
        int inNumBndEdges)
// Constructor
// POST: Creates a mesh
{
    NumElems      = inNumElems;
    NumPhysElems  = inNumPhysElems;
    NumGhostElems = NumElems - NumPhysElems;
    NumNodes      = inNumNodes;
    NumPhysNodes  = inNumPhysNodes; 
    NumBndNodes   = inNumBndNodes;
    NumEdges      = inNumEdges;
    NumBndEdges   = inNumBndEdges;

    if (NumElems<1 || NumPhysElems<1 || NumGhostElems<0 || NumNodes<1 || 
            NumPhysNodes<1 || NumBndNodes<0  || NumEdges<1)
    {
        cout << endl;
        cout << " Error in mesh constructor ... "       << endl;
        cout << "         NumElems = " << NumElems      << endl;
        cout << "     NumPhysElems = " << NumPhysElems  << endl;
        cout << "    NumGhostElems = " << NumGhostElems << endl;
        cout << "         NumNodes = " << NumNodes      << endl;
        cout << "     NumPhysNodes = " << NumPhysNodes  << endl;
        cout << "      NumBndNodes = " << NumBndNodes   << endl;
        cout << "         NumEdges = " << NumEdges      << endl;
        cout << "      NumBndEdges = " << NumBndEdges   << endl;
        cout << endl;
        exit(1);
    }

    // node: list of x & y coordinates of nodes (size = NumNodes-by-2)
    node = new dTensor2(NumNodes,2);

    // edge: list of coordinates (x1,y1) and (x2,y2) that make up edge (size = NumEdges-by-4)
    edge = new dTensor2(NumEdges,4);

    // tnode: list of nodes attached to element (size = NumNodes-by-3)
    tnode = new iTensor2(NumElems,3);

    // tedge: list of edges attached to element (size = NumElems-by-3)
    tedge = new iTensor2(NumElems,3);

    // tedge_orientation: list of orientiation (either +1 or -1) for edges attached to element (size = NumElems-by-3)
    tedge_orientation = new iTensor2(NumElems,3);

    // eelem: list of elements attached to edge (size = NumEdges-by-2)
    eelem = new iTensor2(NumEdges,2);

    // enode: list of nodes attached to edge (size = NumEdges-by-2)
    enode = new iTensor2(NumEdges,2);

    // list of elements on opposite side of boundary to ghost cell (size = NumGhostElems)
    int tmp;
    if (NumGhostElems==0)
    { tmp = 1; }
    else
    { tmp = NumGhostElems; }
    ghost_link = new iTensor1(tmp);

    // list of nodes on opposite side of boundary to external node
    if ((NumNodes-NumPhysNodes)<1)
    { tmp = 1; }
    else
    { tmp = NumNodes-NumPhysNodes; }
    ext_node_link = new iTensor1(tmp);

    // list of nodes that lie on boundary (size = NumBndNodes)
    if (NumBndNodes==0)
    { tmp = 1; }
    else
    { tmp = NumBndNodes; }
    bnd_node = new iTensor1(tmp);

    // list of edges that lie on boundary (size = NumBndEdges)
    if (NumBndEdges==0)
    { tmp = 1; }
    else
    { tmp = NumBndEdges; }
    bnd_edge = new iTensor1(tmp);

    // element areas (size = NumElems)
    area_prim = new dTensor1(NumElems);

    // dual element areas (size = NumNodes)
    area_dual = new dTensor1(NumNodes);

    // Jacobian matrix (size = NumElems-by-2-by-2)
    Jmat = new dTensor3(NumElems,2,2);

    // list of number of elements per node
    NumElemsPerNode = new iTensor1(NumNodes);

    // adjacent: list of elements that are adjacent to current element
    adjacent = new iTensor2(NumElems,3);

    // Skel2Elem: matrix to convert a vector defined on the
    // mesh skeleton to an element-center quantity
    skel2elem_on = false;
    Skel2Elem = new dTensor3(NumElems,20,20);

    // KMI: matrix to convert a vector node information to
    // edge-based quantities
    kmi_on = false;
    KMI = new dTensor3(NumEdges,3,6);

    // Submesh information
    is_submesh = false;
}

mesh::mesh(const mesh& amesh)
    // Copy constructor
    // POST: New mesh created with size and contents same as amesh
{
    NumElems      = amesh.NumElems;
    NumPhysElems  = amesh.NumPhysElems;
    NumGhostElems = amesh.NumGhostElems;
    NumNodes      = amesh.NumNodes;
    NumPhysNodes  = amesh.NumPhysNodes;
    NumBndNodes   = amesh.NumBndNodes;
    NumEdges      = amesh.NumEdges;
    NumBndEdges   = amesh.NumBndEdges;

    for (int i=1; i<=NumNodes; i++)
        for (int j=1; j<=2; j++)
        {
            node->set(i,j, amesh.node->get(i,j) );
        }

    for (int i=1; i<=NumEdges; i++)
        for (int j=1; j<=4; j++)
        {
            edge->set(i,j, amesh.edge->get(i,j) );
        }

    for (int i=1; i<=NumElems; i++)
        for (int j=1; j<=3; j++)
        {
            tnode->set(i,j, amesh.tnode->get(i,j) );
        }

    for (int i=1; i<=NumElems; i++)
        for (int j=1; j<=3; j++)
        {
            tedge->set(i,j, amesh.tedge->get(i,j) );
        }

    for (int i=1; i<=NumElems; i++)
        for (int j=1; j<=3; j++)
        {
            tedge_orientation->set(i,j, amesh.tedge_orientation->get(i,j) );
        }

    for (int i=1; i<=NumEdges; i++)
        for (int j=1; j<=2; j++)
        {
            eelem->set(i,j, amesh.eelem->get(i,j) );
        }

    for (int i=1; i<=NumEdges; i++)
        for (int j=1; j<=2; j++)
        {
            enode->set(i,j, amesh.enode->get(i,j) );
        }

    for (int i=1; i<=NumGhostElems; i++)
    {
        ghost_link->set(i, amesh.ghost_link->get(i) );
    }

    for (int i=1; i<=(NumNodes-NumPhysNodes); i++)
    {
        ext_node_link->set(i, amesh.ext_node_link->get(i) );
    }

    for (int i=1; i<=NumBndNodes; i++)
    {
        bnd_node->set(i, amesh.bnd_node->get(i) );
    }

    for (int i=1; i<=NumBndEdges; i++)
    {
        bnd_edge->set(i, amesh.bnd_edge->get(i) );
    }

    for (int i=1; i<=NumElems; i++)
    {
        area_prim->set(i, amesh.area_prim->get(i) );
    }

    for (int i=1; i<=NumNodes; i++)
    {
        area_dual->set(i, amesh.area_dual->get(i) );
    } 

    for (int i=1; i<=NumElems; i++)
        for (int m1=1; m1<=2; m1++)
            for (int m2=1; m2<=2; m2++)
            {
                Jmat->set(i,m1,m2, amesh.Jmat->get(i,m1,m2) );
            }

    for (int i=1; i<=NumNodes; i++)
    {
        NumElemsPerNode->set(i, amesh.NumElemsPerNode->get(i) );
    }

    for (int i=1; i<=NumElems; i++)
        for (int m=1; m<=3; m++)
        {
            adjacent->set(i,m, amesh.adjacent->get(i,m) );
        } 

    for (int i=1; i<=NumEdges; i++)
        for (int m=1; m<=3; m++)
            for (int n=1; n<=6; n++)
            {
                KMI->set(i,m,n, amesh.KMI->get(i,m,n) );
            }

    // Submesh information
    if (amesh.is_submesh)
    {
        is_submesh      = amesh.is_submesh;
        SubFactor       = amesh.SubFactor;
        SubNumPhysElems = amesh.SubNumPhysElems;
        SubNumPhysNodes = amesh.SubNumPhysNodes;
        SubNumBndNodes  = amesh.SubNumBndNodes;

        for (int i=1; i<=SubNumPhysNodes; i++)
            for (int j=1; j<=2; j++)
            {
                sub_node->set(i,j, amesh.sub_node->get(i,j) );
            }

        for (int i=1; i<=SubNumPhysElems; i++)
            for (int j=1; j<=3; j++)
            {
                sub_tnode->set(i,j, amesh.sub_tnode->get(i,j) );
            }

        for (int i=1; i<=SubNumBndNodes; i++)
        {
            sub_bnd_node->set(i, amesh.sub_bnd_node->get(i) );
        }

        for (int i=1; i<=SubNumPhysElems; i++)
        {
            sub_area_prim->set(i, amesh.sub_area_prim->get(i) );
        }

        for (int i=1; i<=NumPhysElems; i++)
            for (int k=1; k<=(SubFactor*SubFactor); k++)
            {
                elem_subs->set(i,k, amesh.elem_subs->get(i,k) );
            }

        for (int i=1; i<=NumPhysElems; i++)
            for (int k=1; k<=(((SubFactor+1)*(SubFactor+2))/2); k++)
            {
                node_subs->set(i,k, amesh.node_subs->get(i,k) );
            }

    }

}

mesh::~mesh()
    // Destructor
    // POST: mesh no longer exists
{
}

void mesh::OutputMesh(string outputdir)
    // Ouput all mesh information
{
    string fname;

    // output: constants
    fname = outputdir + "/mesh_params.dat";
    ofstream write1(fname.c_str(), ios::out);
    write1 << setw(8) << NumElems      << "  : NumElems "      << endl;
    write1 << setw(8) << NumPhysElems  << "  : NumPhysElems "  << endl;
    write1 << setw(8) << NumGhostElems << "  : NumGhostElems " << endl;
    write1 << setw(8) << NumNodes      << "  : NumNodes "      << endl;
    write1 << setw(8) << NumPhysNodes  << "  : NumPhysNodes "  << endl;
    write1 << setw(8) << NumBndNodes   << "  : NumBndNodes "   << endl;
    write1 << setw(8) << NumEdges      << "  : NumEdges "      << endl;
    write1 << setw(8) << NumBndEdges   << "  : NumBndEdges "   << endl;
    write1 << setw(8) << is_submesh    << "  : is_submesh "    << endl;
    write1.close();
    write1.clear();

    // output: NODE
    fname = outputdir + "/mesh_node.dat";
    ofstream write2(fname.c_str(), ios::out);
    write2 << setprecision(16); 
    for (int i=1; i<=NumNodes; i++)
    {
        write2 << setw(24) << scientific << node->get(i,1);
        write2 << "  ";
        write2 << setw(24) << scientific << node->get(i,2) << endl;
    }
    write2.close();
    write2.clear();

    // output: EDGE
    fname = outputdir + "/mesh_edge.dat";
    ofstream write3(fname.c_str(), ios::out);
    write3 << setprecision(16); 
    for (int i=1; i<=NumEdges; i++)
    {
        write3 << setw(24) << scientific << edge->get(i,1);
        write3 << "  ";
        write3 << setw(24) << scientific << edge->get(i,2);
        write3 << "  ";
        write3 << setw(24) << scientific << edge->get(i,3);
        write3 << "  ";
        write3 << setw(24) << scientific << edge->get(i,4) << endl;
    }
    write3.close();
    write3.clear();

    // output: TNODE
    fname = outputdir + "/mesh_tnode.dat";
    ofstream write4(fname.c_str(), ios::out);
    for (int i=1; i<=NumElems; i++)
    {
        write4 << setw(8) << tnode->get(i,1);
        write4 << "  ";
        write4 << setw(8) << tnode->get(i,2);
        write4 << "  ";
        write4 << setw(8) << tnode->get(i,3) << endl;
    }
    write4.close();
    write4.clear();

    // output: TEDGE
    fname = outputdir + "/mesh_tedge.dat";
    ofstream write5a(fname.c_str(), ios::out);
    for (int i=1; i<=NumElems; i++)
    {
        write5a << setw(8) << tedge->get(i,1);
        write5a << "  ";
        write5a << setw(8) << tedge->get(i,2);
        write5a << "  ";
        write5a << setw(8) << tedge->get(i,3) << endl;
    }
    write5a.close();
    write5a.clear();

    // output: TEDGE_ORIENTATION
    fname = outputdir + "/mesh_tedge_orientation.dat";
    ofstream write5b(fname.c_str(), ios::out);
    for (int i=1; i<=NumElems; i++)
    {
        write5b << setw(8) << tedge_orientation->get(i,1);
        write5b << "  ";
        write5b << setw(8) << tedge_orientation->get(i,2);
        write5b << "  ";
        write5b << setw(8) << tedge_orientation->get(i,3) << endl;
    }
    write5b.close();
    write5b.clear();

    // output: EELEM
    fname = outputdir + "/mesh_eelem.dat";
    ofstream write6(fname.c_str(), ios::out);
    for (int i=1; i<=NumEdges; i++)
    {
        write6 << setw(8) << eelem->get(i,1);
        write6 << "  ";
        write6 << setw(8) << eelem->get(i,2) << endl;
    }
    write6.close();
    write6.clear();

    // output: ENODE
    fname = outputdir + "/mesh_enode.dat";
    ofstream write7(fname.c_str(), ios::out);
    for (int i=1; i<=NumEdges; i++)
    {
        write7 << setw(8) << enode->get(i,1);
        write7 << "  ";
        write7 << setw(8) << enode->get(i,2) << endl;
    }
    write7.close();
    write7.clear();

    // output: GHOST_LINK
    fname = outputdir + "/mesh_ghost_link.dat";
    ofstream write8(fname.c_str(), ios::out);
    for (int i=1; i<=NumGhostElems; i++)
    {
        write8 << setw(8) << ghost_link->get(i) << endl;
    }
    write8.close();
    write8.clear();

    // output: EXT_NODE_LINK
    fname = outputdir + "/mesh_ext_node_link.dat";
    ofstream write8b(fname.c_str(), ios::out);
    for (int i=1; i<=(NumNodes-NumPhysNodes); i++)
    {
        write8b << setw(8) << ext_node_link->get(i) << endl;
    }
    write8b.close();
    write8b.clear();

    // output: BND_NODE
    fname = outputdir + "/mesh_bnd_node.dat";
    ofstream write9(fname.c_str(), ios::out);
    for (int i=1; i<=NumBndNodes; i++)
    {
        write9 << setw(8) << bnd_node->get(i) << endl;
    }
    write9.close();
    write9.clear();

    // output: BND_NODE
    fname = outputdir + "/mesh_bnd_edge.dat";
    ofstream write9b(fname.c_str(), ios::out);
    for (int i=1; i<=NumBndEdges; i++)
    {
        write9b << setw(8) << bnd_edge->get(i) << endl;
    }
    write9b.close();
    write9b.clear();

    // output: AREA_PRIM
    fname = outputdir + "/mesh_area_prim.dat";
    ofstream write10(fname.c_str(), ios::out);
    write10 << setprecision(16); 
    for (int i=1; i<=NumElems; i++)
    {
        write10 << setw(24) << scientific << area_prim->get(i) << endl;
    }
    write10.close();
    write10.clear();

    // output: AREA_DUAL
    fname = outputdir + "/mesh_area_dual.dat";
    ofstream write11(fname.c_str(), ios::out);
    write11 << setprecision(16); 
    for (int i=1; i<=NumNodes; i++)
    {
        write11 << setw(24) << scientific << area_dual->get(i) << endl;
    }
    write11.close();
    write11.clear();

    // output: JMAT
    fname = outputdir + "/mesh_jmat.dat";
    ofstream write12(fname.c_str(), ios::out);
    write12 << setprecision(16); 
    for (int i=1; i<=NumElems; i++)
        for (int m1=1; m1<=2; m1++)
            for (int m2=1; m2<=2; m2++)
            {
                write12 << setw(24) << scientific << Jmat->get(i,m1,m2) << endl;
            }
    write12.close();
    write12.clear();

    // output: NUMELEMSPERNODE
    fname = outputdir + "/mesh_numelemspernode.dat";
    ofstream write13(fname.c_str(), ios::out);
    write13 << setprecision(16); 
    for (int i=1; i<=NumNodes; i++)
    {
        write13 << setw(8) << NumElemsPerNode->get(i) << endl;
    }
    write13.close();
    write13.clear();

    // output: ADJACENT
    fname = outputdir + "/mesh_adjacent.dat";
    ofstream write14(fname.c_str(), ios::out);
    write14 << setprecision(16); 
    for (int i=1; i<=NumElems; i++)
    {
        write14 << setw(8) << adjacent->get(i,1) << "  ";
        write14 << setw(8) << adjacent->get(i,2) << "  ";
        write14 << setw(8) << adjacent->get(i,3) << endl;
    }
    write14.close();
    write14.clear();

    // output: SKEL2ELEM
    if (skel2elem_on)
    {
        fname = outputdir + "/mesh_skel2elem.dat";
        ofstream write15(fname.c_str(), ios::out);
        write15 << setprecision(16); 
        for (int i=1; i<=NumElems; i++)
            for (int m=1; m<=20; m++)
                for (int n=1; n<=20; n++)
                {
                    write15 << setw(24) << scientific << Skel2Elem->get(i,m,n) << endl;
                }    
        write15.close();
        write15.clear();
    }

    // output: KMI
    if (kmi_on)
    {
        fname = outputdir + "/mesh_kmi.dat";
        ofstream write16(fname.c_str(), ios::out);
        write16 << setprecision(16); 
        for (int i=1; i<=NumEdges; i++)
            for (int m=1; m<=3; m++)
                for (int n=1; n<=6; n++)
                {
                    write16 << setw(24) << scientific << KMI->get(i,m,n) << endl;
                }    
        write16.close();
        write16.clear();
    }

    // sub-mesh
    if (is_submesh)
    {  OutputSubMesh(outputdir);  }
}

void mesh::InputMesh(string inputdir)
    // Input all mesh information
{
    string fname;
    int tmp_int;
    double tmp_double;
    char buffer[256];
    int NumElems_in,NumPhysElems_in,NumGhostElems_in; 
    int NumNodes_in,NumPhysNodes_in,NumBndNodes_in;
    int NumEdges_in,NumBndEdges_in;

    // input: constants
    fname = inputdir + "/mesh_params.dat";
    ifstream read1(fname.c_str(), ios::in);
    read1 >> NumElems_in;        read1.getline(buffer,256);
    read1 >> NumPhysElems_in;    read1.getline(buffer,256);
    read1 >> NumGhostElems_in;   read1.getline(buffer,256);
    read1 >> NumNodes_in;        read1.getline(buffer,256);
    read1 >> NumPhysNodes_in;    read1.getline(buffer,256);
    read1 >> NumBndNodes_in;     read1.getline(buffer,256);
    read1 >> NumEdges_in;        read1.getline(buffer,256);
    read1 >> NumBndEdges_in;     read1.getline(buffer,256);
    read1 >> is_submesh;         read1.getline(buffer,256);
    read1.close();   
    read1.clear();

    // Error check
    if (NumElems_in!=NumElems || 
            NumPhysElems_in!=NumPhysElems ||
            NumGhostElems_in!=NumGhostElems ||
            NumNodes_in!=NumNodes ||
            NumPhysNodes_in!=NumPhysNodes ||
            NumBndNodes_in!=NumBndNodes ||
            NumEdges_in!=NumEdges ||
            NumBndEdges_in!=NumBndEdges)
    {
        cout << " ERROR: trying to read-in mesh information " << endl;
        cout << "        with sizes that do not match ... " << endl;
        cout << endl;
        cout << "      NumElems_in = " << NumElems_in      << "          NumElems = " << NumElems << endl;
        cout << "  NumPhysElems_in = " << NumPhysElems_in  << "      NumPhysElems = " << NumPhysElems << endl;
        cout << " NumGhostElems_in = " << NumGhostElems_in << "     NumGhostElems = " << NumGhostElems << endl;
        cout << "      NumNodes_in = " << NumNodes_in      << "          NumNodes = " << NumNodes << endl;
        cout << "  NumPhysNodes_in = " << NumPhysNodes_in  << "      NumPhysNodes = " << NumPhysNodes << endl;
        cout << "   NumBndNodes_in = " << NumBndNodes_in   << "       NumBndNodes = " << NumBndNodes << endl;
        cout << "      NumEdges_in = " << NumEdges_in      << "          NumEdges = " << NumEdges << endl;
        cout << "   NumBndEdges_in = " << NumBndEdges_in   << "       NumBndEdges = " << NumBndEdges << endl;
        cout << endl;
        exit(1);
    }

    // input: NODE
    fname = inputdir + "/mesh_node.dat";
    ifstream read2(fname.c_str(), ios::in);
    for (int i=1; i<=NumNodes; i++)
    {
        read2 >> tmp_double;
        node->set(i,1, tmp_double );
        read2 >> tmp_double;
        node->set(i,2, tmp_double );
    }
    read2.close();
    read2.clear();

    // input: EDGE
    fname = inputdir + "/mesh_edge.dat";
    ifstream read3(fname.c_str(), ios::in);
    for (int i=1; i<=NumEdges; i++)
    {
        read3 >> tmp_double;
        edge->set(i,1, tmp_double );

        read3 >> tmp_double;
        edge->set(i,2, tmp_double );

        read3 >> tmp_double;
        edge->set(i,3, tmp_double );

        read3 >> tmp_double;
        edge->set(i,4, tmp_double );
    }
    read3.close();
    read3.clear();

    // input: TNODE
    fname = inputdir + "/mesh_tnode.dat";
    ifstream read4(fname.c_str(), ios::in);
    for (int i=1; i<=NumElems; i++)
    {
        read4 >> tmp_int;
        tnode->set(i,1, tmp_int );

        read4 >> tmp_int;
        tnode->set(i,2, tmp_int );

        read4 >> tmp_int;
        tnode->set(i,3, tmp_int );
    }
    read4.close();
    read4.clear();

    // input: TEDGE
    fname = inputdir + "/mesh_tedge.dat";
    ifstream read5a(fname.c_str(), ios::in);
    for (int i=1; i<=NumElems; i++)
    {
        read5a >> tmp_int;
        tedge->set(i,1, tmp_int );

        read5a >> tmp_int;
        tedge->set(i,2, tmp_int );

        read5a >> tmp_int;
        tedge->set(i,3, tmp_int );
    }
    read5a.close();
    read5a.clear();

    // input: TEDGE_ORIENTATION
    fname = inputdir + "/mesh_tedge_orientation.dat";
    ifstream read5b(fname.c_str(), ios::in);
    for (int i=1; i<=NumElems; i++)
    {
        read5b >> tmp_int;
        tedge_orientation->set(i,1, tmp_int );

        read5b >> tmp_int;
        tedge_orientation->set(i,2, tmp_int );

        read5b >> tmp_int;
        tedge_orientation->set(i,3, tmp_int );
    }
    read5b.close();
    read5b.clear();

    // input: EELEM
    fname = inputdir + "/mesh_eelem.dat";
    ifstream read6(fname.c_str(), ios::in);
    for (int i=1; i<=NumEdges; i++)
    {
        read6 >> tmp_int;
        eelem->set(i,1, tmp_int );

        read6 >> tmp_int;
        eelem->set(i,2, tmp_int );
    }
    read6.close();
    read6.clear();

    // input: ENODE
    fname = inputdir + "/mesh_enode.dat";
    ifstream read7(fname.c_str(), ios::in);
    for (int i=1; i<=NumEdges; i++)
    {
        read7 >> tmp_int;
        enode->set(i,1, tmp_int );

        read7 >> tmp_int;
        enode->set(i,2, tmp_int );
    }
    read7.close();
    read7.clear();

    // input: GHOST_LINK
    fname = inputdir + "/mesh_ghost_link.dat";
    ifstream read8(fname.c_str(), ios::in);
    for (int i=1; i<=NumGhostElems; i++)
    {
        read8 >> tmp_int;
        ghost_link->set(i, tmp_int );
    }
    read8.close();
    read8.clear();

    // input: EXT_NODE_LINK
    fname = inputdir + "/mesh_ext_node_link.dat";
    ifstream read8b(fname.c_str(), ios::out);
    for (int i=1; i<=(NumNodes-NumPhysNodes); i++)
    {
        read8b >> tmp_int;
        ext_node_link->set(i, tmp_int );
    }
    read8b.close();
    read8b.clear();

    // input: BND_NODE
    fname = inputdir + "/mesh_bnd_node.dat";
    ifstream read9(fname.c_str(), ios::in);
    for (int i=1; i<=NumBndNodes; i++)
    {
        read9 >> tmp_int;
        bnd_node->set(i, tmp_int );
    }
    read9.close();
    read9.clear();

    // input: BND_EDGE
    fname = inputdir + "/mesh_bnd_node.dat";
    ifstream read9b(fname.c_str(), ios::in);
    for (int i=1; i<=NumBndEdges; i++)
    {
        read9b >> tmp_int;
        bnd_edge->set(i, tmp_int );
    }
    read9b.close();
    read9b.clear();

    // input: AREA_PRIM
    fname = inputdir + "/mesh_area_prim.dat";
    ifstream read10(fname.c_str(), ios::in);
    for (int i=1; i<=NumElems; i++)
    {
        read10 >> tmp_double;
        area_prim->set(i, tmp_double );
    }
    read10.close();
    read10.clear();

    // input: AREA_DUAL
    fname = inputdir + "/mesh_area_dual.dat";
    ifstream read11(fname.c_str(), ios::in);
    for (int i=1; i<=NumNodes; i++)
    {
        read11 >> tmp_double;
        area_dual->set(i, tmp_double );
    }
    read11.close();
    read11.clear();

    // input: JMAT
    fname = inputdir + "/mesh_jmat.dat";
    ifstream read12(fname.c_str(), ios::in);
    for (int i=1; i<=NumElems; i++)
        for (int m1=1; m1<=2; m1++)
            for (int m2=1; m2<=2; m2++)
            {
                read12 >> tmp_double;
                Jmat->set(i,m1,m2, tmp_double );
            }
    read12.close();
    read12.clear();

    // input: NUMELEMSPERNODE
    fname = inputdir + "/mesh_numelemspernode.dat";
    ifstream read13(fname.c_str(), ios::in);
    for (int i=1; i<=NumNodes; i++)
    {
        read13 >> tmp_int;
        NumElemsPerNode->set(i, tmp_int );
    }
    read13.close();
    read13.clear();

    // input: ADJACENT
    fname = inputdir + "/mesh_adjacent.dat";
    ifstream read14(fname.c_str(), ios::in);
    for (int i=1; i<=NumElems; i++)
    {
        read14 >> tmp_int;
        adjacent->set(i,1, tmp_int );
        read14 >> tmp_int;
        adjacent->set(i,2, tmp_int );
        read14 >> tmp_int;
        adjacent->set(i,3, tmp_int );
    }
    read14.close();
    read14.clear();

    // input: SKEL2ELEM
    fname = inputdir + "/mesh_skel2elem.dat";
    ifstream read15(fname.c_str(), ios::in);
    skel2elem_on = false;
    if (read15.is_open())
    {      
        skel2elem_on = true;
        for (int i=1; i<=NumElems; i++)
            for (int m1=1; m1<=20; m1++)
                for (int m2=1; m2<=20; m2++)
                {
                    read15 >> tmp_double;
                    Skel2Elem->set(i,m1,m2, tmp_double );
                }
        read15.close();
        read15.clear();
    }

    // input: KMI
    fname = inputdir + "/mesh_kmi.dat";
    ifstream read16(fname.c_str(), ios::in);
    kmi_on = false;
    if (read16.is_open())
    {      
        kmi_on = true;
        for (int i=1; i<=NumEdges; i++)
            for (int m1=1; m1<=3; m1++)
                for (int m2=1; m2<=6; m2++)
                {
                    read16 >> tmp_double;
                    KMI->set(i,m1,m2, tmp_double );
                }
        read16.close();
        read16.clear();
    }

    // sub-mesh
    if (is_submesh)
    {  InputSubMesh(inputdir);  }
}

const int& mesh::get_NumElems() const
// Returns "NumElems"
{
    return NumElems;
}

const int& mesh::get_NumPhysElems() const
// Returns "NumPhysElems"
{
    return NumPhysElems;
}

const int& mesh::get_NumGhostElems() const
// Returns "NumGhostElems"
{
    return NumGhostElems;
}

const int& mesh::get_NumNodes() const
// Returns "NumNodes"
{
    return NumNodes;
}

const int& mesh::get_NumPhysNodes() const
// Returns "NumPhysNodes"
{
    return NumPhysNodes;
}

const int& mesh::get_NumBndNodes() const
// Returns "NumBndNodes"
{
    return NumBndNodes;
}

const int& mesh::get_NumBndEdges() const
// Returns "NumBndEdges"
{
    return NumBndEdges;
}

const int& mesh::get_NumEdges() const
// Returns "NumEdges"
{
    return NumEdges;
}

const double& mesh::get_node(int i,int m) const
// Returns x (m=1) or y (m=2) coordinate of ith node
{
    return node->get(i,m);
}

const double& mesh::get_edge(int i,int m) const
// Returns (x1,y1) and (x2,y2) values of edge:
//     m = 1  --->  x1
//     m = 2  --->  y1
//     m = 3  --->  x2
//     m = 4  --->  y2
{
    return edge->get(i,m);
}

const int& mesh::get_tnode(int i,int m) const
// Returns pointer to node for ith element (m=1,2, or 3)
{
    return tnode->get(i,m);
}

const int& mesh::get_tedge(int i,int m) const
// Returns pointer to edge for ith element (m=1,2, or 3)
{
    return tedge->get(i,m);
}

const int& mesh::get_tedge_orientation(int i,int m) const
// Returns +1 if normal to edge points out of ith element (m=1,2, or 3) -- outward pointing normal
// Returns -1 if normal to edge points into   ith element (m=1,2, or 3) -- inward pointing normal
{
    return tedge_orientation->get(i,m);
}

const int& mesh::get_adjacent(int i, int m) const
// Returns pointer to neighboring elements for ith element (m=1,2, or 3)
{
    return adjacent->get(i,m);
}

const int& mesh::get_eelem(int i,int m) const
// Returns pointer to element for ith edge (m=1 or 2)
{
    return eelem->get(i,m);
}

const int& mesh::get_enode(int i,int m) const
// Returns pointer to node for ith edge (m=1 or 2)
{
    return enode->get(i,m);
}

const int& mesh::get_ghost_link(int i) const
// Returns pointer to element on opposite side of boundary to ghost cell
{
    return ghost_link->get(i);
}

const int& mesh::get_ext_node_link(int i) const
// Returns pointer to node on opposite side of boundary to external node
{
    return ext_node_link->get(i);
}

const int& mesh::get_bnd_node(int i) const
// Returns pointer to node on boundary
{
    return bnd_node->get(i);
}

const double& mesh::get_area_prim(int i) const
// Returns area of element i
{
    return area_prim->get(i);
}

const double& mesh::get_area_dual(int i) const
// Returns area of dual element centered at node i
{
    return area_dual->get(i);
}

const double& mesh::get_jmat(int i, int m1, int m2) const
// Returns (m1,m2) component of Jacobian matrix in element i
{
    return Jmat->get(i,m1,m2);
}

const int& mesh::get_NumElemsPerNode(int i) const
// Returns the number of elements attached to node i
{
    return NumElemsPerNode->get(i);
}

void mesh::set_node(int i,int m, double input)
    // Sets x (m=1) or y (m=2) coordinate of ith node
{
    node->set(i,m, input );
}

void mesh::set_edge(int i,int m, double input)
    // Sets (x1,y1) and (x2,y2) values of edge:
    //     m = 1  --->  x1
    //     m = 2  --->  y1
    //     m = 3  --->  x2
    //     m = 4  --->  y2
{
    edge->set(i,m, input );
}

void mesh::set_tnode(int i,int m, int input)
    // Sets pointer to node for ith element (m=1,2, or 3)
{
    tnode->set(i,m, input );
}

void mesh::set_tedge(int i,int m, int input)
    // Sets pointer to edge for ith element (m=1,2, or 3)
{
    tedge->set(i,m, input );
}

void mesh::set_tedge_orientation(int i,int m, int input)
    // Set value to +1 if normal to edge points out of ith element (m=1,2, or 3) -- outward pointing normal
    // Set value to -1 if normal to edge points into   ith element (m=1,2, or 3) -- inward pointing normal
{
    if (input!=1 && input!=-1 && input!=0)
    {
        printf("\n");
        printf(" ERROR in set_tedge_orientation, input must be set to 0, -1 or +1. \n");
        printf(" input = %i\n",input);
        printf("\n");
        exit(1);
    }
    tedge_orientation->set(i,m, input );
}

void mesh::set_eelem(int i,int m, int input)
    // Sets pointer to element for ith edge (m=1 or 2)
{
    eelem->set(i,m, input );
}

void mesh::set_enode(int i,int m, int input)
    // Sets pointer to node for ith edge (m=1 or 2)
{
    enode->set(i,m, input );
}

void mesh::set_ghost_link(int i, int input)
    // Sets pointer to element on opposite side of boundary to ghost cell
{
    ghost_link->set(i, input );
}

void mesh::set_ext_node_link(int i, int input)
    // Sets pointer to node on opposite side of boundary to external node
{
    ext_node_link->set(i, input );
}

void mesh::set_bnd_node(int i, int input)
    // Sets pointer to node on boundary
{
    bnd_node->set(i, input );
}

void mesh::set_bnd_edge(int i, int input)
    // Sets pointer to edge on boundary
{
    bnd_edge->set(i, input );
}

void mesh::set_area_prim(int i, double input)
    // Sets area of element i
{
    area_prim->set(i, input );
}

void mesh::set_area_dual(int i, double input)
    // Sets area of dual element centered at node i
{
    area_dual->set(i, input );
}

void mesh::Compute_NumElemsPerNode()
    // Computes the number of elements attached to each node
{
    for (int i=1; i<=NumNodes; i++)
    {
        NumElemsPerNode->set(i, 0 );
    }

    for (int j=1; j<=NumElems; j++)
        for (int k=1; k<=3; k++)
        {
            int i = tnode->get(j,k);
            int num_old = NumElemsPerNode->get(i);
            NumElemsPerNode->set(i,num_old+1);
        }
}

void mesh::ComputeJacobian()
    // Compute 2x2 Jacobian transformation matrix for each element i
    // The Jacobian matrix tells one how to transform a gradient in
    //   in the canonical variables (xi,eta) to a gradient in 
    //   physical variables (x,y):
    //
    //     phi_x = Jmat(1,1)*phi_xi + Jmat(1,2)*phi_eta
    //     phi_y = Jmat(2,1)*phi_xi + Jmat(2,2)*phi_eta
{
    // loop over each element
    for (int i=1; i<=NumElems; i++)
    {
        // Find coordinates of current element
        int i1 = tnode->get(i,1);
        int i2 = tnode->get(i,2);
        int i3 = tnode->get(i,3);

        double x1 = node->get(i1,1);
        double y1 = node->get(i1,2);

        double x2 = node->get(i2,1);
        double y2 = node->get(i2,2);

        double x3 = node->get(i3,1);
        double y3 = node->get(i3,2);

        double Area = area_prim->get(i);

        Jmat->set(i,1,1, (y3-y1)/(2.0*Area) );
        Jmat->set(i,1,2, (y1-y2)/(2.0*Area) );
        Jmat->set(i,2,1, (x1-x3)/(2.0*Area) );
        Jmat->set(i,2,2, (x2-x1)/(2.0*Area) );
    }
}


void mesh::SetAdjacency()
    // Find all elements that are adjacent to current element
    // and store in "iTensor2* adjacent"
{
    for (int i=1; i<=NumElems; i++)
    {
        int te,e1,e2;

        te = tedge->get(i,1);
        if (te>0)
        {
            e1 = eelem->get(te,1);
            e2 = eelem->get(te,2);
            if (e1!=i)
            { adjacent->set(i,1, e1 ); }
            else
            { adjacent->set(i,1, e2 ); }
        }
        else
        {
            adjacent->set(i,1, -1 );
        }

        te = tedge->get(i,2);
        if (te>0)
        {
            e1 = eelem->get(te,1);
            e2 = eelem->get(te,2);
            if (e1!=i)
            { adjacent->set(i,2, e1 ); }
            else
            { adjacent->set(i,2, e2 ); }
        }
        else
        {
            adjacent->set(i,2, -1 );
        }

        te = tedge->get(i,3);
        if (te>0)
        {
            e1 = eelem->get(te,1);
            e2 = eelem->get(te,2);
            if (e1!=i)
            { adjacent->set(i,3, e1 ); }
            else
            { adjacent->set(i,3, e2 ); }
        }
        else
        {
            adjacent->set(i,3, -1 );
        }
    }
}

// get elements of matrix to convert a vector defined on the
// mesh skeleton to an element-center quantity
const double& mesh::get_Skel2Elem(int i, int m, int n) const
{
    return Skel2Elem->get(i,m,n);
}

// set elements of matrix to convert a vector defined on the
// mesh skeleton to an element-center quantity
void mesh::set_Skel2Elem(int i, int m, int n, double val)
{
    Skel2Elem->set(i,m,n, val );
}

const bool& mesh::is_skel2elem_on() const
{return skel2elem_on;}

void mesh::turn_skel2elem_on()
{skel2elem_on = true;}

void mesh::turn_skel2elem_off()
{skel2elem_on = false;}

// get elements of matrix to convert a vector node information to
// edge-based quantities
const double& mesh::get_KMI(int i, int m, int n) const
{
    return KMI->get(i,m,n);
}

// set elements of matrix to convert a vector node information to
// edge-based quantities
void mesh::set_KMI(int i, int m, int n, double val)
{
    KMI->set(i,m,n, val );
}

const bool& mesh::is_kmi_on() const
{return kmi_on;}

void mesh::turn_kmi_on()
{kmi_on = true;}

void mesh::turn_kmi_off()
{kmi_on = false;}


// Create sub-mesh
void mesh::CreateSubMesh(int num_divide, 
        double xmin, double xmax,
        double ymin, double ymax,
        double deps, double (*SignedDistance)(point))
{
    assert(is_submesh==false);

    if (num_divide>1)
    {
        is_submesh = true;
        SubFactor = num_divide;      
    }
    else
    {
        printf("\n");
        printf(" Error in Mesh.CreateSubMesh: Invalid num_divide input parameter.\n");
        printf("         num_divide = %i\n",num_divide);
        printf("\n");
        exit(1);
    }

    SubNumElems = NumElems*SubFactor*SubFactor;
    SubNumPhysElems = NumPhysElems*SubFactor*SubFactor;
    elem_subs = new iTensor2(NumPhysElems,SubFactor*SubFactor);
    node_subs = new iTensor2(NumPhysElems,((SubFactor+1)*(SubFactor+2))/2);
    sub_tnode = new iTensor2(SubNumPhysElems,3);
    sub_area_prim = new dTensor1(SubNumPhysElems);

    const int tmp_num_sub_nodes = (((SubFactor+1)*(SubFactor+2))/2)*NumPhysElems;
    int node_count = 0;
    int elem_count = 0;
    const double dF = 1.0/double(SubFactor);
    const double onethird = 1.0/3.0;
    dTensor2 tmp_sub_node(tmp_num_sub_nodes,3);
    dTensor1 node_val(tmp_num_sub_nodes);
    iTensor1 index(tmp_num_sub_nodes);

    int n = 0;
    iTensor2 nindex(SubFactor+1,SubFactor+1);
    iTensor2 eindex(2*SubFactor-1,SubFactor);

    for (int k=1; k<=(SubFactor+1); k++)
        for (int j=1; j<=(SubFactor+2-k); j++)
        {
            n = n+1;
            nindex.set(j,k, n );
        }

    int e = 0;
    for (int k=1; k<=SubFactor; k++)
        for (int j=1; j<=(2*(SubFactor-k)+1); j++)
        {
            e = e+1;
            eindex.set(j,k, e );
        }


    for (int i=1; i<=NumPhysElems; i++)
    {
        // basic information
        const int i1 = tnode->get(i,1);
        const int i2 = tnode->get(i,2);
        const int i3 = tnode->get(i,3);

        const double x1 = node->get(i1,1);
        const double y1 = node->get(i1,2);

        const double x2 = node->get(i2,1);
        const double y2 = node->get(i2,2);

        const double x3 = node->get(i3,1);
        const double y3 = node->get(i3,2);

        const double xc = onethird*(x1+x2+x3);
        const double yc = onethird*(y1+y2+y3);

        // nodes            
        for (int k=1; k<=(SubFactor+1); k++)
            for (int j=1; j<=(SubFactor+2-k); j++)
            {
                node_count = node_count + 1;
                const double xi  = -onethird + (j-1)*dF;
                const double eta = -onethird + (k-1)*dF;

                const double x = xc + xi*(x2-x1) + eta*(x3-x1);
                const double y = yc + xi*(y2-y1) + eta*(y3-y1);

                tmp_sub_node.set(node_count,1, x );
                tmp_sub_node.set(node_count,2, y );

                node_val.set(node_count, ((x-xmin)/(xmax-xmin)) + 10000.0*(1.0+(y-ymin)/(ymax-ymin)) );
                index.set(node_count, node_count );
            }
    }

    QuickSort_double(node_val,index,1,tmp_num_sub_nodes);

    dTensor1 node_val_unique(tmp_num_sub_nodes);
    iTensor1 index_unique(tmp_num_sub_nodes);
    iTensor1 rev_index(tmp_num_sub_nodes);
    for (int k=1; k<=tmp_num_sub_nodes; k++)
    {
        node_val_unique.set(k, 0.0 );
        index_unique.set(k, k );
        rev_index.set(k, k );
    }
    SubNumPhysNodes = Unique_double(node_val,index,node_val_unique,index_unique,rev_index,
            1,tmp_num_sub_nodes);

    sub_node = new dTensor2(SubNumPhysNodes,2);
    for (int i=1; i<=SubNumPhysNodes; i++)
    {
        sub_node->set(i,1, tmp_sub_node.get(index_unique.get(i),1) );
        sub_node->set(i,2, tmp_sub_node.get(index_unique.get(i),2) );
    }

    // elements
    node_count = 0;
    for (int i=1; i<=NumPhysElems; i++)
    {
        int n = 0;      
        for (int k=1; k<=SubFactor; k++)
        {
            int jtmp = 0;
            for (int j=1; j<=(2*(SubFactor-k)+1); j++)
            {
                elem_count = elem_count + 1;
                n = n+1;

                elem_subs->set(i,n, elem_count );

                if (j%2 == 1)
                {
                    jtmp = jtmp + 1;		  
                    sub_tnode->set(elem_count,1, rev_index.get(node_count + nindex.get(jtmp,k  )) );
                    sub_tnode->set(elem_count,2, rev_index.get(node_count + nindex.get(jtmp+1,k)) );
                    sub_tnode->set(elem_count,3, rev_index.get(node_count + nindex.get(jtmp,k+1)) );
                }
                else
                {
                    sub_tnode->set(elem_count,1, rev_index.get(node_count + nindex.get(jtmp+1,k))   );
                    sub_tnode->set(elem_count,2, rev_index.get(node_count + nindex.get(jtmp+1,k+1)) );
                    sub_tnode->set(elem_count,3, rev_index.get(node_count + nindex.get(jtmp,k+1))   );
                }
            }	 
        }
        node_count = node_count + ((SubFactor+1)*(SubFactor+2))/2;
    }

    // for each super-element store pointers to nodes on this super-element
    for (int i=1; i<=NumPhysElems; i++)
    {      
        iTensor1 sub_node_list(3*SubFactor*SubFactor);
        iTensor1 index(3*SubFactor*SubFactor);
        int ncount = 1;
        for (int k=1; k<=(SubFactor*SubFactor); k++)
        {  	  
            int t_tmp = elem_subs->get(i,k);

            sub_node_list.set(ncount,   sub_tnode->get(t_tmp,1) );
            sub_node_list.set(ncount+1, sub_tnode->get(t_tmp,2) );
            sub_node_list.set(ncount+2, sub_tnode->get(t_tmp,3) );

            index.set(ncount,   ncount   );
            index.set(ncount+1, ncount+1 );
            index.set(ncount+2, ncount+2 );

            ncount = ncount+3;
        }
        iTensor1 aout(3*SubFactor*SubFactor);
        iTensor1 iout(3*SubFactor*SubFactor);
        iTensor1 irout(3*SubFactor*SubFactor);
        QuickSort_int(sub_node_list,index,1,3*SubFactor*SubFactor);
        int num_nodes = Unique_int(sub_node_list,index,aout,iout,irout,1,3*SubFactor*SubFactor);

        assert( (((SubFactor+1)*(SubFactor+2))/2) == num_nodes );

        // Set pointer to all sub-nodes that live on a given element
        // First organize to make node numbering easy
        const int i1 = tnode->get(i,1);
        const int i2 = tnode->get(i,2);
        const int i3 = tnode->get(i,3);
        const double x1  = node->get(i1,1);
        const double x2  = node->get(i2,1);
        const double x3  = node->get(i3,1);
        const double y1  = node->get(i1,2);
        const double y2  = node->get(i2,2);
        const double y3  = node->get(i3,2);
        const double xc = onethird*(x1+x2+x3);
        const double yc = onethird*(y1+y2+y3);
        const double T  = area_prim->get(i);
        iTensor1 ns_ind(num_nodes);
        dTensor1 ns_val(num_nodes);
        for (int k=1; k<=num_nodes; k++)
        {
            int ind = aout.get(k);
            double xx = sub_node->get(ind,1);
            double yy = sub_node->get(ind,2);
            double xi  = ((y3-y1)*(xx-xc)+(x1-x3)*(yy-yc))/(2.0*T);
            double eta = ((y1-y2)*(xx-xc)+(x2-x1)*(yy-yc))/(2.0*T);

            ns_ind.set(k, k);
            ns_val.set(k, ((xi+onethird)) + 10000.0*(1.0+(eta+onethird)) ); 	  
        }

        QuickSort_double(ns_val,ns_ind,1,num_nodes);

        for (int k=1; k<=num_nodes; k++)
        {  node_subs->set(i,k, aout.get(ns_ind.get(k)) );  }
    }

    // find all sub-nodes that are on a boundary edge
    iTensor2 project_node(NumBndEdges,SubFactor+1);
    for (int i=1; i<=NumBndEdges; i++)
    {
        int e = bnd_edge->get(i);
        int t1 = eelem->get(e,1);
        int t2 = eelem->get(e,2);
        int n1 = enode->get(e,1);
        int n2 = enode->get(e,2);

        double x1 = node->get(n1,1);
        double y1 = node->get(n1,2);
        double x2 = node->get(n2,1);
        double y2 = node->get(n2,2);

        int t;
        if (t1<=NumPhysElems)
        {  t = t1;  }
        else
        {  t = t2;  }

        iTensor1 sub_node_list(3*SubFactor*SubFactor);
        iTensor1 index(3*SubFactor*SubFactor);
        int ncount = 1;
        for (int k=1; k<=(SubFactor*SubFactor); k++)
        {  	  
            int t_tmp = elem_subs->get(t,k);

            sub_node_list.set(ncount,   sub_tnode->get(t_tmp,1) );
            sub_node_list.set(ncount+1, sub_tnode->get(t_tmp,2) );
            sub_node_list.set(ncount+2, sub_tnode->get(t_tmp,3) );

            index.set(ncount,   ncount   );
            index.set(ncount+1, ncount+1 );
            index.set(ncount+2, ncount+2 );

            ncount = ncount+3;
        }
        iTensor1 aout(3*SubFactor*SubFactor);
        iTensor1 iout(3*SubFactor*SubFactor);
        iTensor1 irout(3*SubFactor*SubFactor);
        QuickSort_int(sub_node_list,index,1,3*SubFactor*SubFactor);
        int num_nodes = Unique_int(sub_node_list,index,aout,iout,irout,1,3*SubFactor*SubFactor);

        int num_on_edge = 0;      
        for (int k=1; k<=num_nodes; k++)
        {
            double x = sub_node->get(aout.get(k),1);
            double y = sub_node->get(aout.get(k),2);

            double f  = (y2-y1)*(x-x1)-(x2-x1)*(y-y1);

            if (fabs(f)<=1.0e-13)
            {
                double d1 = sqrt(pow(x-x1,2)+pow(y-y1,2));
                double d2 = sqrt(pow(x-x2,2)+pow(y-y2,2));
                num_on_edge = num_on_edge+1;
                project_node.set(i,num_on_edge, aout.get(k) );
            }
        }
        if (num_on_edge!=(SubFactor+1))
        {
            printf(" e = %i/%i, num_on_edge = %i,    SubFactor+1 = %i\n",e,NumBndEdges,num_on_edge,SubFactor+1);
        }
        assert(num_on_edge==(SubFactor+1));      
    }

    // find all boundary nodes
    int tmp_num_bnd_nodes = 0;
    iTensor1 tmp_bnd_nodes(NumBndEdges*(SubFactor+1));
    iTensor1 bnd_index(NumBndEdges*(SubFactor+1));
    iTensor1 bnd_iout(NumBndEdges*(SubFactor+1));
    iTensor1 bnd_riout(NumBndEdges*(SubFactor+1));
    iTensor1 bnd_aout(NumBndEdges*(SubFactor+1));
    for (int i=1; i<=NumBndEdges; i++)
        for (int k=1; k<=(SubFactor+1); k++)
        {
            tmp_num_bnd_nodes = tmp_num_bnd_nodes+1;
            tmp_bnd_nodes.set(tmp_num_bnd_nodes,project_node.get(i,k));
            bnd_index.set(tmp_num_bnd_nodes, tmp_num_bnd_nodes );
        }
    //  SubNumBndNodes = tmp_num_bnd_nodes;
    QuickSort_int(tmp_bnd_nodes,bnd_index,1,tmp_num_bnd_nodes);
    SubNumBndNodes = Unique_int(tmp_bnd_nodes,bnd_index,
            bnd_aout,bnd_iout,bnd_riout,
            1,tmp_num_bnd_nodes);
    sub_bnd_node = new iTensor1(SubNumBndNodes);
    for (int i=1; i<=SubNumBndNodes; i++)
    {
        sub_bnd_node->set(i, bnd_aout.get(i) );
    }

    // compute element areas, force CCW orientation
    for (int i=1; i<=SubNumPhysElems; i++)
    {
        const int i1 = sub_tnode->get(i,1);
        const int i2 = sub_tnode->get(i,2);
        const int i3 = sub_tnode->get(i,3);

        double area = 0.5 * ( sub_node->get(i1,1) * (sub_node->get(i2,2)-sub_node->get(i3,2)) + 
                sub_node->get(i2,1) * (sub_node->get(i3,2)-sub_node->get(i1,2)) +
                sub_node->get(i3,1) * (sub_node->get(i1,2)-sub_node->get(i2,2)) );
        if (area < 0.0)
        {
            double temp = sub_tnode->get(i,2);
            sub_tnode->set(i,2, sub_tnode->get(i,3) );
            sub_tnode->set(i,3, temp );	  
        }
        sub_area_prim->set(i, fabs(area) );    
    }

}

const bool& mesh::get_is_submesh() const
{
    return is_submesh;
}

const int& mesh::get_SubFactor() const
{
    assert(is_submesh==true);
    return SubFactor;
}

const int& mesh::get_SubNumPhysElems() const
{
    assert(is_submesh==true);
    return SubNumPhysElems;
}

const int& mesh::get_SubNumPhysNodes() const
{
    assert(is_submesh==true);
    return SubNumPhysNodes;
}

const int& mesh::get_SubNumBndNodes() const
{
    assert(is_submesh==true);
    return SubNumBndNodes;
}

const double& mesh::get_sub_node(int i, int j) const
{
    assert(is_submesh==true);
    return sub_node->get(i,j);
}

const int& mesh::get_sub_tnode(int i, int j) const
{
    assert(is_submesh==true);
    return sub_tnode->get(i,j);
}

const int& mesh::get_sub_bnd_node(int i) const
{
    assert(is_submesh==true);
    return sub_bnd_node->get(i);
}

const double& mesh::get_sub_area_prim(int i) const
{
    assert(is_submesh==true);
    return sub_area_prim->get(i);
}

const int& mesh::get_elem_subs(int i, int j) const
{
    assert(is_submesh==true);
    return elem_subs->get(i,j);
}

const int& mesh::get_node_subs(int i, int j) const
{
    assert(is_submesh==true);
    return node_subs->get(i,j);
}

void mesh::OutputSubMesh(string outputdir)
{
    string fname;

    // output: constants
    fname = outputdir + "/submesh_params.dat";
    ofstream write1(fname.c_str(), ios::out);
    write1 << setw(8) << SubFactor        << "  : SubFactor       " << endl;
    write1 << setw(8) << SubNumPhysElems  << "  : SubNumPhysElems " << endl;
    write1 << setw(8) << SubNumPhysNodes  << "  : SubNumPhysNodes " << endl;
    write1 << setw(8) << SubNumBndNodes   << "  : SubNumBndNodes  " << endl;
    write1.close();
    write1.clear();

    // output: SUBNODE
    fname = outputdir + "/submesh_node.dat";
    ofstream write2(fname.c_str(), ios::out);
    write2 << setprecision(16); 
    for (int i=1; i<=SubNumPhysNodes; i++)
    {
        write2 << setw(24) << scientific << sub_node->get(i,1);
        write2 << "  ";
        write2 << setw(24) << scientific << sub_node->get(i,2) << endl;
    }
    write2.close();
    write2.clear();

    // output: SUBTNODE
    fname = outputdir + "/submesh_tnode.dat";
    ofstream write3(fname.c_str(), ios::out);
    for (int i=1; i<=SubNumPhysElems; i++)
    {
        write3 << setw(8) << sub_tnode->get(i,1);
        write3 << "  ";
        write3 << setw(8) << sub_tnode->get(i,2);
        write3 << "  ";
        write3 << setw(8) << sub_tnode->get(i,3) << endl;
    }
    write3.close();
    write3.clear();

    // output: SUB_BND_NODE
    fname = outputdir + "/submesh_bnd_node.dat";
    ofstream write4(fname.c_str(), ios::out);
    for (int i=1; i<=SubNumBndNodes; i++)
    {
        write4 << setw(8) << sub_bnd_node->get(i) << endl;
    }
    write4.close();
    write4.clear();

    // output: SUB_AREA_PRIM
    fname = outputdir + "/submesh_area_prim.dat";
    ofstream write5(fname.c_str(), ios::out);
    write5 << setprecision(16); 
    for (int i=1; i<=SubNumPhysElems; i++)
    {
        write5 << setw(24) << scientific << sub_area_prim->get(i) << endl;
    }
    write5.close();
    write5.clear();

    // output: ELEM_SUBS
    fname = outputdir + "/submesh_elem_subs.dat";
    ofstream write6(fname.c_str(), ios::out);
    for (int i=1; i<=NumPhysElems; i++)
    {
        for (int j=1; j<=(SubFactor*SubFactor); j++)
        {
            write6 << setw(8) << elem_subs->get(i,j);
        }
        write6 << endl;
    }
    write6.close();
    write6.clear();

    // output: NODE_SUBS
    fname = outputdir + "/submesh_node_subs.dat";
    ofstream write7(fname.c_str(), ios::out);
    for (int i=1; i<=NumPhysElems; i++)
    {
        for (int j=1; j<=(((SubFactor+1)*(SubFactor+2))/2); j++)
        {
            write7 << setw(8) << node_subs->get(i,j);
        }
        write7 << endl;
    }
    write7.close();
    write7.clear();
}

void mesh::InputSubMesh(string inputdir)
{
    string fname;
    int tmp_int;
    double tmp_double;
    char buffer[256];

    // input: constants
    fname = inputdir + "/submesh_params.dat";
    ifstream read1(fname.c_str(), ios::in);
    read1 >> SubFactor;        read1.getline(buffer,256);
    read1 >> SubNumPhysElems;  read1.getline(buffer,256);
    read1 >> SubNumPhysNodes;  read1.getline(buffer,256);
    read1 >> SubNumBndNodes;   read1.getline(buffer,256);
    read1.close();   
    read1.clear();    

    // dimension arrays
    sub_node      = new dTensor2(SubNumPhysNodes,2);
    sub_tnode     = new iTensor2(SubNumPhysElems,3);
    sub_bnd_node  = new iTensor1(SubNumBndNodes);
    sub_area_prim = new dTensor1(SubNumPhysElems);
    elem_subs     = new iTensor2(NumPhysElems,SubFactor*SubFactor);
    node_subs     = new iTensor2(NumPhysElems,((SubFactor+1)*(SubFactor+2))/2);

    // input: SUBNODE
    fname = inputdir + "/submesh_node.dat";
    ifstream read2(fname.c_str(), ios::in);
    for (int i=1; i<=SubNumPhysNodes; i++)
    {
        read2 >> tmp_double;
        sub_node->set(i,1, tmp_double );
        read2 >> tmp_double;
        sub_node->set(i,2, tmp_double );
    }
    read2.close();
    read2.clear();

    // input: SUBTNODE
    fname = inputdir + "/submesh_tnode.dat";
    ifstream read3(fname.c_str(), ios::in);
    for (int i=1; i<=SubNumPhysElems; i++)
    {
        read3 >> tmp_int;
        sub_tnode->set(i,1, tmp_int );
        read3 >> tmp_int;
        sub_tnode->set(i,2, tmp_int );
        read3 >> tmp_int;
        sub_tnode->set(i,3, tmp_int );
    }
    read3.close();
    read3.clear();

    // input: SUB_BND_NODE
    fname = inputdir + "/submesh_bnd_node.dat";
    ifstream read4(fname.c_str(), ios::in);
    for (int i=1; i<=SubNumBndNodes; i++)
    {
        read4 >> tmp_int;
        sub_bnd_node->set(i, tmp_int );
    }
    read4.close();
    read4.clear();

    // input: SUB_AREA_PRIM
    fname = inputdir + "/submesh_area_prim.dat";
    ifstream read5(fname.c_str(), ios::in);
    for (int i=1; i<=SubNumPhysElems; i++)
    {
        read5 >> tmp_double;
        sub_area_prim->set(i, tmp_double );
    }
    read5.close();
    read5.clear();

    // input: ELEM_SUBS
    fname = inputdir + "/submesh_elem_subs.dat";
    ifstream read6(fname.c_str(), ios::in);
    for (int i=1; i<=NumPhysElems; i++)
        for (int j=1; j<=(SubFactor*SubFactor); j++)
        {
            read6 >> tmp_int;
            elem_subs->set(i,j, tmp_int );
        }
    read6.close();
    read6.clear(); 

    // input: NODE_SUBS
    fname = inputdir + "/submesh_node_subs.dat";
    ifstream read7(fname.c_str(), ios::in);
    for (int i=1; i<=NumPhysElems; i++)
        for (int j=1; j<=(((SubFactor+1)*(SubFactor+2))/2); j++)
        {
            read7 >> tmp_int;
            node_subs->set(i,j, tmp_int );
        }
    read7.close();
    read7.clear(); 
}


void mesh::QuickSort_int(iTensor1& a, iTensor1& index, int lo, int hi)
{
    //  lo is the lower index, hi is the upper index
    //  hi the region of array a that is to be sorted
    int i = lo;
    int j = hi;
    int x=a.get( (i+j)/2 );
    int h;
    int itmp;

    //  partition
    while(i<=j) 
    {           
        while (a.get(i)<x) 
        {i++;} 

        while (a.get(j)>x) 
        {j--;}

        if (i<=j)
        {
            h=a.get(i);
            a.set(i,a.get(j));
            a.set(j,h);

            itmp = index.get(i);
            index.set(i, index.get(j) );
            index.set(j, itmp );

            i++; j--;
        }
    }

    //  recursion
    if (lo<j) QuickSort_int(a, index, lo, j);
    if (i<hi) QuickSort_int(a, index, i, hi);
}

int mesh::Unique_int(const iTensor1& a_in, 
        const iTensor1& index_in,
        iTensor1& a_out, 
        iTensor1& index_out,
        iTensor1& rev_index,
        int lo, 
        int hi)
{
    int NumUniqueNodes = 1;


    int curr = a_in.get(lo);
    int curr_index = index_in.get(lo);
    a_out.set(NumUniqueNodes,curr);
    index_out.set(NumUniqueNodes, curr_index );
    rev_index.set(index_in.get(lo), 1 );

    for (int k=(lo+1); k<=hi; k++)
    {      
        int next = a_in.get(k);

        if (curr!=next)
        {
            NumUniqueNodes = NumUniqueNodes + 1;
            curr = next;
            a_out.set(NumUniqueNodes, curr );
            curr_index = index_in.get(k);
            index_out.set(NumUniqueNodes, curr_index );
        }
        rev_index.set(index_in.get(k), NumUniqueNodes );
    }

    return NumUniqueNodes;
}


void mesh::QuickSort_double(dTensor1& a, iTensor1& index, int lo, int hi)
{
    //  lo is the lower index, hi is the upper index
    //  hi the region of array a that is to be sorted
    int i = lo;
    int j = hi;
    double x=a.get( (i+j)/2 );
    double h;
    int itmp;

    //  partition
    while(i<=j) 
    {           
        while (a.get(i)<x) 
        {i++;} 

        while (a.get(j)>x) 
        {j--;}

        if (i<=j)
        {
            h=a.get(i);
            a.set(i,a.get(j));
            a.set(j,h);

            itmp = index.get(i);
            index.set(i, index.get(j) );
            index.set(j, itmp );

            i++; j--;
        }
    }

    //  recursion
    if (lo<j) QuickSort_double(a, index, lo, j);
    if (i<hi) QuickSort_double(a, index, i, hi);
}

int mesh::Unique_double(const dTensor1& a_in, 
        const iTensor1& index_in,
        dTensor1& a_out, 
        iTensor1& index_out,
        iTensor1& rev_index,
        int lo, 
        int hi)
{
    int NumUniqueNodes = 1;

    double curr = a_in.get(lo);
    int curr_index = index_in.get(lo);
    a_out.set(NumUniqueNodes,curr);
    index_out.set(NumUniqueNodes, curr_index );
    rev_index.set(index_in.get(lo), 1 );
    double next;

    for (int k=(lo+1); k<=hi; k++)
    {
        curr = a_in.get(k-1);
        next = a_in.get(k);

        if (fabs(next-curr)>1.0e-10)
        {
            NumUniqueNodes = NumUniqueNodes + 1;
            a_out.set(NumUniqueNodes, next );
            curr_index = index_in.get(k);
            index_out.set(NumUniqueNodes, curr_index );
        }
        rev_index.set(index_in.get(k), NumUniqueNodes );
    }

    return NumUniqueNodes;
}
