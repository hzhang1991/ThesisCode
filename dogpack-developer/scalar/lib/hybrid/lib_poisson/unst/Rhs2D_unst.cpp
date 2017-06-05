#include "tensors.h"
#include "mesh.h"

// Create right-hand side vector, when the rhs is already known in terms of
// the dg-basis elements
void Rhs2D_unst(const int space_order,
        const mesh& Mesh, const dTensor3& rhs_dg,
        dTensor1& rhs)
{

    // Constants
    const int kmax = ((space_order*(space_order+1))/2);
    const int NumElems = Mesh.get_NumElems();
    const int NumPhysElems = Mesh.get_NumPhysElems();

    assert_eq( rhs_dg.getsize(2), 1 );

    // Project DG solution onto CG basis functions
    void L2Project_DG2CG_Unst(const int istart, 
            const int iend, 
            const int QuadOrder,		    
            const int cg_order,
            const int dg_comp,
            const mesh& Mesh,
            const dTensor3& fdg,
            dTensor1& fcg);
    L2Project_DG2CG_Unst(1,NumPhysElems,space_order+1,space_order,
            1,Mesh,rhs_dg,rhs);

    // Use right-hand side to also enforce homogeneous Dirichlet boundary conditions
    switch( space_order )
    {
        case 1:
            for (int i=1; i<=Mesh.get_NumBndNodes(); i++)
            {
                int j = Mesh.get_bnd_node(i);
                rhs.set(j, 0.0 );
            }
            break;
        case 2:
            for (int i=1; i<=Mesh.get_SubNumBndNodes(); i++)
            {
                int j = Mesh.get_sub_bnd_node(i);
                rhs.set(j, 0.0 );
            }
            break;
        case 3:
            for (int i=1; i<=Mesh.get_SubNumBndNodes(); i++)
            {
                int j = Mesh.get_sub_bnd_node(i);
                rhs.set(j, 0.0 );
            }
            break;
    }

}

// Create right-hand side vector
void Rhs2D_unst(const int space_order,
        const mesh& Mesh,
        dTensor1& rhs)
{

printf("WARNING: do you really want to project RhsFunc onto basis functions?\n");

    // Constants
    const int kmax = ((space_order*(space_order+1))/2);
    const int NumElems = Mesh.get_NumElems();
    const int NumPhysElems = Mesh.get_NumPhysElems();

    // Project rhs function onto DG basis functions
    void RhsFunc(const dTensor2& xpts, 
            const dTensor2& NOT_USED_1,
            const dTensor2& NOT_USED_2, 
            dTensor2& rhs);
    void L2Project_Unst(const int istart, 
            const int iend, 
            const int QuadOrder, 
            const int BasisOrder_qin,
            const int BasisOrder_auxin,
            const int BasisOrder_fout,
            const mesh& Mesh, 
            const dTensor3* qin, 
            const dTensor3* auxin, 
            dTensor3* fout, 
            void (*Func)(const dTensor2&,const dTensor2&,
                const dTensor2&,dTensor2&));

    dTensor3 rhs_dg(NumElems,1,kmax);
    L2Project_Unst(1,NumElems,
            space_order,space_order,
            space_order,space_order,
            Mesh,&rhs_dg,&rhs_dg,
            &rhs_dg, &RhsFunc);

    // Project DG solution onto CG basis functions
    void L2Project_DG2CG_Unst(const int istart, 
            const int iend, 
            const int QuadOrder,		    
            const int cg_order,
            const int dg_comp,
            const mesh& Mesh,
            const dTensor3& fdg,
            dTensor1& fcg);
    L2Project_DG2CG_Unst(1,NumPhysElems,space_order+1,space_order,
            1,Mesh,rhs_dg,rhs);

    // Use right-hand side to also enforce homogeneous Dirichlet boundary conditions
    switch( space_order )
    {
        case 1:
            for (int i=1; i<=Mesh.get_NumBndNodes(); i++)
            {
                int j = Mesh.get_bnd_node(i);
                rhs.set(j, 0.0 );
            }
            break;
        case 2:
            for (int i=1; i<=Mesh.get_SubNumBndNodes(); i++)
            {
                int j = Mesh.get_sub_bnd_node(i);
                rhs.set(j, 0.0 );
            }
            break;
        case 3:
            for (int i=1; i<=Mesh.get_SubNumBndNodes(); i++)
            {
                int j = Mesh.get_sub_bnd_node(i);
                rhs.set(j, 0.0 );
            }
            break;
    }

}
