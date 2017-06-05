#ifndef _VLASOVPARAMS_H_
#define _VLASOVPARAMS_H_

class IniDocument;
class VlasovParams
{   
    public:

        // These two parameters are used to set the background ions
        double rho0   ;  // initial background density for maxwellian
        double area   ;  // area (of configuration space)
//      double temp   ;  // temperature of plasma

//      // Parameters used for the beam problem:
//      double eta    ;  // Ratio of omega / omega0
//      double omega0 ;
//      double n0     ;
//      double vth2   ;  // Square of the thermal velocity
//      double a      ;
//      double R2     ;  // R^2

        //void init(IniDocument& ini_doc);
        void init();

        // Destructor
        ~VlasovParams();

        // Constructor:
        VlasovParams():is_initialized(false){}

    private:
        bool is_initialized;

};

extern VlasovParams vlasovParams;
#endif
