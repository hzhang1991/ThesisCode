#ifndef _ACOUSTICPARAMS_H_
#define _ACOUSTICPARAMS_H_

class IniDocument;
class AcousticParams
{   
    public:

        void init(IniDocument& ini_doc);

        // methods 
        const double& get_c(void);
        void   set_c(double c);

    private:

        double c  ;  // speed of light ( c > 0 )

};
extern AcousticParams acousticParams;
#endif
