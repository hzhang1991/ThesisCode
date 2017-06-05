// --------------------------------------------------------------------------
//  SPECIFICATION FILE (GridPatch.h)
//    A class that describes a 1d grid patch
// --------------------------------------------------------------------------


// 1D GRID PATCH ----------------------------------------
class GridPatch
{
 public:
    GridPatch(int,double,double,int,int,int,int);
    // Constructor
    
    GridPatch(const GridPatch& anotherGridPatch);
    // Copy constructor
    
    ~GridPatch();
    // Destructor

    void getq(dTensorBC3& qout);
    // Put the values of "q" into "qout"

    void getaux(dTensorBC3& auxout);
    // Put the values of "aux" into "auxout"

    void setq(dTensorBC3 qin);
    // Put the values of "q" into "qout"

    void setaux(dTensorBC3 auxin);
    // Put the values of "aux" into "auxout"

    double getxlow();
    // get parameter xlow

    double getxhigh();
    // get parameter xhigh

    double getdx();
    // get parameter dx

    int getmx();
    // get parameter mx

    void setxlow(double);
    // set parameter xlow

    void setxhigh(double);
    // set parameter xhigh

    void resize(int);
    // reset mx and modify q and aux
    // according to this new mx

    int getmeqn();
    // get parameter meqn
    
    int getmaux();
    // get parameter maux

    int getkmax();
    // get parameter kmax

    int getmbc();
    // get parameter mbc from dTensorBC3 q

 private:
    int mx,meqn,maux,kmax,mbc;
    double xlow,xhigh;
    dTensorBC3   q;
    dTensorBC3 aux;
};
// ---------------------------------------------------------
