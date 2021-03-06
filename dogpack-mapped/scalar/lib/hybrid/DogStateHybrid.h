#ifndef _DogStateHybrid_h_
#define _DogStateHybrid_h_

// This class currently only contains three doubles, all related to time
// quantities.

// Why isn't this a struct like DogStateCart2 is?  Or rather, why isn't
// DogStateCart2 a class?
//
// If this is to follow the same class hierarchy as the structured code, this
// should inherit from DogStateTB
class DogStateHybrid
{

    private:
        double time;
        double dt;
        double initial_dt;

    public:

        // The initial values for all private variables are set to zero.
        void init();

        // Constructors:
        DogStateHybrid();
        DogStateHybrid(DogStateHybrid&);

        // Desctructor:
        ~DogStateHybrid();

        // read-access methods
        //
        double get_time() const { return time; }
        double get_dt()   const { return dt; }
        double get_initial_dt() const { return initial_dt; }

        // write-access methods
        //
        void set_time(double in)       { time = in;}
        void set_dt(double in)         { dt = in; }
        void set_initial_dt(double in) { initial_dt = in; }
};

extern DogStateHybrid dogStateHybrid;

#endif
