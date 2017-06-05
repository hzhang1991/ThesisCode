#include "tensors.h"
#include "RKinfo.h"

void SetRKinfo(int method2, RKinfo& rk)
{
    rk.mstage = 0;

    switch( method2 )
    {
        case 5:
            rk.num_stages = 8;
            break;
        case 4:
            rk.num_stages = 10;
            break;
        default:
            rk.num_stages = method2;
    }

    rk.alpha1 = new dTensor1(rk.num_stages);
    rk.alpha2 = new dTensor1(rk.num_stages);
    rk.beta   = new dTensor1(rk.num_stages);

    // 5th order stuff 
    rk.delta  = new dTensor1(rk.num_stages);
    rk.gamma  = new dTensor2(3, rk.num_stages);

    switch(method2)
    {
        case 1: // first-order

            rk.alpha1->set(1, 1.0 );
            rk.alpha2->set(1, 0.0 );
            rk.beta->set(  1, 1.0 );

            break;

        case 2: // second-order

            rk.alpha1->set(1, 1.0 );
            rk.alpha2->set(1, 0.0 );
            rk.beta->set(  1, 1.0 );

            rk.alpha1->set(2, 0.5 );
            rk.alpha2->set(2, 0.5 );
            rk.beta->set(  2, 0.5 );

            break;

        case 3: // third-order

            rk.alpha1->set(1, 1.0 );
            rk.alpha2->set(1, 0.0 );
            rk.beta->set(  1, 1.0 );

            rk.alpha1->set(2, 0.75 );
            rk.alpha2->set(2, 0.25 );
            rk.beta->set(  2, 0.25 );

            rk.alpha1->set(3, 2.0e0/3.0e0 );
            rk.alpha2->set(3, 1.0e0/3.0e0 );
            rk.beta->set(  3, 2.0e0/3.0e0 );

            break;

        case 4: // fourth-order

            for (int i=1; i<=9; i++)
            {
                rk.alpha1->set(i, 1.0 );
                rk.alpha2->set(i, 0.0 );
                rk.beta->set(  i, 1.0/6.0 );
            }

            rk.alpha1->set(10, 1.0 );
            rk.alpha2->set(10, 3.0/5.0 );
            rk.beta->set(  10, 1.0/10.0 );

            break;

        case 5: // 5th-order

            rk.delta->set(1, 1.0e0);
            rk.delta->set(2, 1.528486658778845e00);
            rk.delta->set(3, 4.720094096662784e-02);
            rk.delta->set(4, 8.801244253465348e-01);
            rk.delta->set(5, 1.019066090228480e+00);
            rk.delta->set(6, 1.049772291176110e+01);
            rk.delta->set(7, -4.254616508506826e+00);
            rk.delta->set(8, 0.0);

            rk.gamma->set(1, 1, 0.0);
            rk.gamma->set(1, 2, -1.552288007713033e+01);
            rk.gamma->set(1, 3,  4.127286635722417e-01);
            rk.gamma->set(1, 4, -1.011819196331377e+00);
            rk.gamma->set(1, 5, -2.765748383780848e-01);
            rk.gamma->set(1, 6,  5.075770311217778e-02);
            rk.gamma->set(1, 7,  6.999810478513669e+00);
            rk.gamma->set(1, 8, -1.114908881433104e+01);

            rk.gamma->set(2, 1,  1.0);
            rk.gamma->set(2, 2,  6.534691420958578e+00);
            rk.gamma->set(2, 3,  2.280056542904473e-01);
            rk.gamma->set(2, 4,  1.308684311397668e+00);
            rk.gamma->set(2, 5,  4.769419552531064e-01);
            rk.gamma->set(2, 6, -6.368809762042849e-03);
            rk.gamma->set(2, 7,  9.339446057238532e-02);
            rk.gamma->set(2, 8,  9.556626047962331e-01);

            rk.gamma->set(3, 1, 0.0);
            rk.gamma->set(3, 2, 0.0);
            rk.gamma->set(3, 3, 0.0);
            rk.gamma->set(3, 4, -2.510747784045939e+00);
            rk.gamma->set(3, 5, -8.576822794622042e-01);
            rk.gamma->set(3, 6,  1.044599944472272e+00);
            rk.gamma->set(3, 7, -7.000810861049136e+00);
            rk.gamma->set(3, 8,  1.906311811144179e+00);

            rk.beta->set(1,  8.653258038183180e-02);
            rk.beta->set(2,  9.544677980851571e-01);
            rk.beta->set(3,  2.651941386774408e-01);
            rk.beta->set(4,  2.736914413910379e-01);
            rk.beta->set(5,  5.999778649323600e-01);
            rk.beta->set(6,  4.433177471748104e-03);
            rk.beta->set(7,  5.309971130968292e-03);
            rk.beta->set(8,  5.830861806762871e-01);
            break;

    }
}

void DeleteRKInfo(RKinfo& rk)
{
    delete rk.alpha1;
    delete rk.alpha2;
    delete rk.beta;
    delete rk.gamma;
    delete rk.delta;
}
