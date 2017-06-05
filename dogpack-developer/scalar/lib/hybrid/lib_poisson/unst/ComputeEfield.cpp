#include "mesh.h"
#include "constants.h"
#include "dogdefs.h"
#include "MonomialsToLegendre.h"
#include <cmath>

// This routine take the derivative of phi and saves the calculated values onto
// E1 = phi_x, and E2 = phi_y.
//
//   E1(1:NumPhysElems, 1:kmax_unst )
//   E2(1:NumPhysElems, 1:kmax_unst )
//
void ComputeEfield(const int space_order,
        const mesh& Mesh,		   
        const dTensor1& phi,
        dTensor2& E1,
        dTensor2& E2)
{

    const int NumPhysElems = Mesh.get_NumPhysElems();
    const int NumPhysNodes = Mesh.get_NumPhysNodes();
    const int kmax = E1.getsize(2);
    assert(kmax==E2.getsize(2));
    assert(kmax==((space_order*(space_order+1))/2));

    switch(space_order)
    {
        case 1:

            assert_eq(kmax,1);

            for (int i=1; i<=NumPhysElems; i++)
            {
                int n1 = Mesh.get_tnode(i,1);
                int n2 = Mesh.get_tnode(i,2);
                int n3 = Mesh.get_tnode(i,3);

                double phi1 = phi.get(n1);
                double phi2 = phi.get(n2);
                double phi3 = phi.get(n3);

                double phi_xi  = phi2 - phi1;
                double phi_eta = phi3 - phi1;

                double E1tmp = Mesh.get_jmat(i,1,1)*phi_xi + Mesh.get_jmat(i,1,2)*phi_eta;	  
                double E2tmp = Mesh.get_jmat(i,2,1)*phi_xi + Mesh.get_jmat(i,2,2)*phi_eta;

                E1.set(i,1, -E1tmp );
                E2.set(i,1, -E2tmp );
            }
            break;

        case 2:

            assert_eq(kmax,3);

            for (int i=1; i<=NumPhysElems; i++)
            {
                int n1 = Mesh.get_node_subs(i,1);
                int n2 = Mesh.get_node_subs(i,2);
                int n3 = Mesh.get_node_subs(i,3);
                int n4 = Mesh.get_node_subs(i,4);
                int n5 = Mesh.get_node_subs(i,5);
                int n6 = Mesh.get_node_subs(i,6);

                double phi1 = phi.get(n1);
                double phi2 = phi.get(n2);
                double phi3 = phi.get(n3);
                double phi4 = phi.get(n4);
                double phi5 = phi.get(n5);
                double phi6 = phi.get(n6);

                double A1 = onethird*(phi2+phi4+phi5);
                double A2 = (sq2/60.0)*(-3.0*phi1+4.0*phi2+6.0*phi3-8.0*phi4+4.0*phi5-3.0*phi6);
                double A3 = (sq2*sq3/60.0)*(-3.0*phi1-4.0*phi2+4.0*phi5+3.0*phi6);
                double A4 = (sq7/105.0)*(phi1-phi2-3.0*phi3-phi4+7.0*phi5-3.0*phi6);
                double A5 = (sq3*sq7/315.0)*(5.0*phi1-12.0*phi2+6.0*phi3+2.0*phi4-phi6);
                double A6 = (sq3*sq5/45.0)*(phi1-2.0*phi4+phi6);

                dTensor1 phi_xi(3);
                phi_xi.set(1,  3.0*sq2*A2 + sq2*sq3*A3 + 8.0*osq7*A4 + 2.0*osq3*sq5*A6 - 2.0*osq3*osq7*A5 );
                phi_xi.set(2, -5.0*osq7*sq2*A4 - osq2*osq3*sq5*A6 + 55.0*osq2*osq3*osq7*A5 );
                phi_xi.set(3,  5.0*osq7*sq2*sq3*A4 + 3.0*osq2*sq5*A6 + 15.0*osq2*osq7*A5 );

                dTensor1 phi_eta(3);
                phi_eta.set(1, 8.0*osq7*A4 + 4.0*osq7*sq3*A5 + 2.0*sq2*sq3*A3 );
                phi_eta.set(2, 10.0*osq7*sq2*A4 + 5.0*osq7*sq2*sq3*A5 );
                phi_eta.set(3, 3.0*sq2*sq5*A6 );

                for (int k=1; k<=3; k++)
                {
                    E1.set(i,k, -(phi_xi.get(k)*Mesh.get_jmat(i,1,1) + phi_eta.get(k)*Mesh.get_jmat(i,1,2)) );
                    E2.set(i,k, -(phi_xi.get(k)*Mesh.get_jmat(i,2,1) + phi_eta.get(k)*Mesh.get_jmat(i,2,2)) );
                }
            }
            break;

        case 3:

            assert_eq(kmax,6);      

            for (int i=1; i<=NumPhysElems; i++)
            {
                int n1  = Mesh.get_node_subs(i,1);
                int n2  = Mesh.get_node_subs(i,2);
                int n3  = Mesh.get_node_subs(i,3);
                int n4  = Mesh.get_node_subs(i,4);
                int n5  = Mesh.get_node_subs(i,5);
                int n6  = Mesh.get_node_subs(i,6);
                int n7  = Mesh.get_node_subs(i,7);
                int n8  = Mesh.get_node_subs(i,8);
                int n9  = Mesh.get_node_subs(i,9);
                int n10 = Mesh.get_node_subs(i,10);

                double phi1  = phi.get(n1);
                double phi2  = phi.get(n2);
                double phi3  = phi.get(n3);
                double phi4  = phi.get(n4);
                double phi5  = phi.get(n5);
                double phi6  = phi.get(n6);
                double phi7  = phi.get(n7);
                double phi8  = phi.get(n8);
                double phi9  = phi.get(n9);
                double phi10 = phi.get(n10);	  

                double A1  =  onethird*0.0025*( 40.0*(phi1+phi4+phi10) 
                        + 90.0*(phi2+phi3+phi5+phi7+phi8+phi9) + 540.0*phi6 );
                double A2  = -onethird*0.05*osq2*(phi1+phi10-2.0*phi4+9.0*(phi2+phi5+phi8+phi9)-18.0*(phi3+phi7) );
                double A3  = -onethird*0.05*osq2*sq3*(phi1-phi10 + 9.0*(phi2+phi5-phi8-phi9));
                double A4  =  oneseventh*0.075*osq7*(phi2+phi5+4.0*phi1+17.0*(phi7+phi9)
                        -12.0*(phi4+phi10)+30.0*phi6-23.0*(phi3+phi8));
                double A5  =  oneseventh*0.05*sq3*osq7*(2.0*(phi8-phi10)+10.0*phi1-15.0*phi2+12.0*(phi4-phi3)
                        +18.0*phi7-3.0*phi9+20.0*phi5-30.0*phi6);
                double A6  =  oneseventh*0.05*sq3*sq5*(2.0*(phi1-phi5-phi8+phi10)+3.0*(phi2+phi9)-6.0*phi6);
                double A7  = -0.0375*oneseventh*sq3*osq13*(3.0*phi1-2.0*(phi3+phi5)-7.0*(phi2+phi8)
                        +18.0*(phi4+phi6)-12.0*phi10+43.0*phi9-52.0*phi7);
                double A8  = -0.1875*oneseventh*sq3*osq5*osq13*osq19*(19.0*phi1-30.0*(phi3+phi5)-27.0*(phi2+phi8)
                        +10*phi4+28.0*phi10+114.0*phi6-57.0*phi9);
                double A9  = -0.0375*oneseventh*sq3*osq19*(19.0*phi1+3.0*phi5-3.0*phi8-20.0*phi4
                        -60.0*phi2+phi10+60.0*phi3);
                double A10 = -0.1125*osq7*(3.0*(phi8-phi5)+phi1-phi10);

                dTensor1 phi_xi(6);
                phi_xi.set(1, 3.0*sq2*A2 + sq2*sq3*A3 + 8.0*osq7*A4 - 2.0*osq3*osq7*A5 
                        + 2.0*osq3*sq5*A6 - 5.0*sq3*osq13*A7 + 3.0*sq3*sq5*osq13*osq19*A8 
                        + 13.0*sq3*osq19*A9 + sq7*A10 );
                phi_xi.set(2, -5.0*sq2*osq7*A4 + 55.0*osq2*osq3*osq7*A5 - osq2*osq3*sq5*A6 
                        + 7.0*sq2*sq3*osq13*A7 + 31.0*sq2*sq3*osq5*osq13*osq19*A8 
                        + 2.8*sq2*sq3*osq19*A9 - 0.4*sq2*sq7*A10 );
                phi_xi.set(3, 5.0*sq2*sq3*osq7*A4 + 15.0*osq2*osq7*A5 + 3.0*osq2*sq5*A6 
                        - 9.0*sq2*osq13*A7 + 79.0*sq2*osq5*osq13*osq19*A8 - 8.8*sq2*osq19*A9 
                        + 0.8*sq2*sq3*sq7*A10 );
                phi_xi.set(4, 7.0*sq3*sq7*osq13*A7 + 29.0*osq3*sq5*sq7*osq13*osq19*A8 
                        - osq3*sq7*osq19*A9 - A10 );
                phi_xi.set(5, -4.0*onethird*osq5*sq7*sq13*osq19*A8 
                        + 66.4*onethird*sq7*osq19*A9 - 0.8*osq3*A10 );
                phi_xi.set(6,  28.0*onethird*sq13*osq19*A8 + 14.0*onethird*sq5*osq19*A9 
                        + 2.0*osq3*sq5*sq7*A10 );

                dTensor1 phi_eta(6);
                phi_eta.set(1, 2.0*sq2*sq3*A3 + 8.0*osq7*A4 + 4.0*sq3*osq7*A5+5*sq3*osq13*A7 
                        - 3.0*sq3*sq5*osq13*osq19*A8 + 6.0*sq3*osq19*A9 + 2.0*sq7*A10 );
                phi_eta.set(2, 10.0*sq2*osq7*A4 + 5.0*sq2*sq3*osq7*A5 + 8.0*sq2*sq3*osq13*A7 
                        - 24.0*sq2*sq3*osq5*osq13*osq19*A8 + 9.6*sq2*sq3*osq19*A9 - 0.8*sq2*sq7*A10 );
                phi_eta.set(3, 3.0*sq2*sq5*A6 + 8.0*sq2*osq5*sq13*osq19*A8 + 4.0*sq2*osq19*A9 );
                phi_eta.set(4, -3.0*sq3*sq7*osq13*A7 + 47.0*osq3*sq5*sq7*osq13*osq19*A8 
                        - 2.8*osq3*sq7*osq19*A9 - 0.4*A10 );
                phi_eta.set(5, 6.0*sq7*osq13*A7 + 86.0*osq5*sq7*osq13*osq19*A8 
                        + 11.2*sq7*osq19*A9 + 0.8*osq3*A10 );
                phi_eta.set(6, 4.0*osq3*sq5*sq7*A10 );

                for (int k=1; k<=6; k++)
                {
                    E1.set(i,k, -(phi_xi.get(k)*Mesh.get_jmat(i,1,1) + phi_eta.get(k)*Mesh.get_jmat(i,1,2)) );
                    E2.set(i,k, -(phi_xi.get(k)*Mesh.get_jmat(i,2,1) + phi_eta.get(k)*Mesh.get_jmat(i,2,2)) );
                }
            }
            break;

        default:
            printf("\n");
            printf(" Error in ComputeEfield.cpp: space_order value is not implemented ...\n");
            printf("     space_order = %i\n",space_order);
            printf("\n");
            exit(1);
            break;
    }

}
