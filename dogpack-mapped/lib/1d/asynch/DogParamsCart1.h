#ifndef DOGPARAMSCART1_H
#define DOGPARAMSCART1_H

class IniDocument;
struct DogParamsCart1{
private:
  bool is_initialized;
  int mx;
  int mbc;
  double xlow;
  double xhigh;
  double dx; // calculated
  bool read_grid;
  
  //
  // methods
  //
private:
  void checkParameters();
  void setDerivedParameters();
public:
  bool get_is_initialized(){return is_initialized;}
  void reportParameters();
  
  DogParamsCart1(){is_initialized=false;}
  
  void init(IniDocument& ini_doc);
  
  const int   & get_mx()    const{ return  mx;   }
  const int   & get_melems()    const{ return  mx;   }
  const int   & get_mbc()   const{ return  mbc;  }
  const double& get_xlow()  const{ return  xlow; }
  const double& get_xhigh() const{ return  xhigh;}
  const double& get_dx()    const{ return  dx;   }
  const bool& get_read_grid() const{ return read_grid; }
  
  void set_mx(int mx_in){ mx = mx_in; }
  void set_xlims(double,double);
};

extern DogParamsCart1 dogParamsCart1;

#endif
