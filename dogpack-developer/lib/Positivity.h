#ifndef _POSITIVITY_H__
#define _POSITIVITY_H__

// Class used to describe a "positive" solution.
//
// Currently only used in the 2D-Cartesian plasma parts of the code.  
// (Last checked 2/11/2013 -DS)
class PositivityState
{
    private:
        enum StringencyLevel{
            AVERAGE=0, // checking positivity of cell average
            RIEMANN_POINTS=1, // checking positivity at riemann points
            POSITIVITY_POINTS=2, // enforcing positivity at positivity points
            // using implicit source term to guarantee positivity
            // for sufficiently short time step independent of solution
            IMPLICIT_SOURCE=3,
        };

    private:
        // positivity indicators
        //
        // minimum positivity indicator value
        double minval;
        // used with minval to suggest new positivity CFL
        //double dt;
        //
        // stringency of positivity enforcement
        //
        StringencyLevel stringencyLevel;
        // factor by which desired cfl is being rescaled
        // in order to enforce positivity.
        double cflFactor;
        double dt;
        //
        int stepsWithoutViolation;
        //
        // flags
        //
        bool repeatStageFlag;
        bool repeatStepFlag;
    public:
        PositivityState();
    private:
        void resetStringency()
        {
            stringencyLevel = POSITIVITY_POINTS;
            stepsWithoutViolation=0;
            cflFactor = 1.;
        }
        void increaseStringency();
        void set_dt(double in) { dt = in; }
        void set_minval(double in) { minval = in; }
    public:
        void update(double minval_in, double dt_in);
        void clearFlags()
        {
            repeatStageFlag = false;
            repeatStepFlag = false;
        }
        double get_cflFactor() const {return cflFactor;}
        double get_dt() const {return dt;}
        bool get_repeatStepFlag() const { return repeatStepFlag; }
        bool get_repeatStageFlag() const { return repeatStageFlag; }
        //bool using_split_source() const;
        void write_state(void* file)const;
        void read_state(void* file);
};

#endif
