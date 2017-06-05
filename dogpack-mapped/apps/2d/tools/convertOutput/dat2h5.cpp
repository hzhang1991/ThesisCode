#include<string>
#include"debug.h"
#include<IniDocument.h>
#include<DogParams.h>
#include<DogParamsCart2.h>
#include<tensors.h>

using namespace std;

int input_meqn;
int meqn;
#define size_type int

enum DataType{STATE_DATA=1, SAMPLE_DATA=2};

void usage()
{
  printf("\n"
         "usage: dat2h5 [-meqns <N>] <file1>.dat <file2>.dat ...\n"
         "\n"
         "  <N>: number of state variables in file\n"
         "  <fileX>: an output file, or a state file if it ends\n"
         "           in '_restart.dat' or '_state.dat'.\n"
         "\n"
  );
  exit(1);
}

bool testSuffix(string str,string suffix)
{
   return
      str.length() > suffix.length()
   && str.substr(str.length()-suffix.length(), suffix.length()) == suffix;
}

string stripSuffix(string str,string suffix)
{
   assert(testSuffix(str,suffix));
   return str.substr(0,str.length()-suffix.length());
}

//DataType getDataType(string filename_prefix)
//{
//  DataType dataType;
//
//  if( testSuffix(filename_prefix,"_restart")
//   || testSuffix(filename_prefix,"_state"))
//  {
//    printf("%s is a state file.\n",filename_prefix.c_str());
//    dataType = STATE_DATA;
//  }
//  else
//  {
//    printf("%s is a sample file.\n",filename_prefix.c_str());
//    dataType = SAMPLE_DATA;
//  }
//  return dataType;
//}

bool hasSlash(const string& str)
{
  size_type pos = str.find_first_of('/');
  if(pos==string::npos) return false;
  return true;
  //const char* file_str = str.c_str();
  //while(*file_str)
  //{
  //  if(*file_str++ == '/') return true;
  //}
  //return false;
}

// returns string::npos if nothing is found
size_type findLastSlash(const string& str)
{
   return str.find_last_of('/');
   //int pos = str.length();
   //while(str[pos]!='/' && pos-->0);
   //if(str[pos]!='/') pos=-1;
   //return pos;
}

string removeFromLastSlash(const string& str)
{
   size_type pos = findLastSlash(str);
   if(pos==string::npos) return "";
   return str.substr(0,pos);
}

// strip directory prefix
string stripDirectoryPrefix(const string& str)
{
  // return everything after last slash.
  size_type pos = findLastSlash(str);
  if(pos==string::npos) return str;
  size_type start = pos+1;
  return str.substr(start,str.length()-start);
}

string getDirectoryOfFile(string filename)
{
  // if path contains a slash, delete everything starting with last slash
  if(hasSlash(filename)) return removeFromLastSlash(filename);
  // else take directory as "."
  return string(".");
}

void dat2h5(string filename)
{

  string filename_prefix = stripSuffix(filename,".dat");
  printf("filename_prefix=[%s]\n",filename_prefix.c_str());

  // get first letter of base file name
  //
  size_type pos = filename.find_last_of('/');
  if(pos==string::npos) pos=0;
  else pos++;
  string state_variable_name = filename.substr(pos,1);
  //string basefilename = stripDirectoryPrefix(filename_prefix);
  //char state_variable_name[2];
  //strcpy(state_variable_name,"q");
  //strncpy(state_variable_name,basefilename.c_str(),1);
  //assert(!strcmp(state_variable_name,"q") || !strcmp(state_variable_name,"a"));
  assert(state_variable_name=="q" || state_variable_name=="a");

  const int mx = dogParamsCart2.get_mx();
  const int my = dogParamsCart2.get_my();
  const int mbc = dogParamsCart2.get_mbc();
  const int kmax = dogParams.get_kmax();
  const int space_order = dogParams.get_space_order();
  //DataType dataType = getDataType(filename_prefix);
  //if(dataType==STATE_DATA)
  //{
  dTensorBC4 q(mx,my,meqn,kmax,mbc);
  double QinitRestartASCII(string fname, dTensorBC4& q);
  double t = QinitRestartASCII(filename_prefix,q);
  void WriteStateHDF5(string fname, string varname,
    		   const dTensorBC4& q, double t);
  WriteStateHDF5(filename_prefix, state_variable_name, q, t);

  //}
  //else
  //{
  //  int x_elements=space_order*mx;
  //  int y_elements=space_order*my;
  //  dTensor3 qvals(x_elements,y_elements,meqn);
  //  double ReadOutputASCII(string fname, dTensor3& qvals);
  //  double t = ReadOutputASCII(filename_prefix,qvals);
  //  void WriteOutputHDF5(string fname, const dTensor3& qvals, double t);
  //  WriteOutputHDF5(filename_prefix, qvals, t);
  //}
}

int main(int argc, char** argv)
{
  rlimit rlim;
  getrlimit(RLIMIT_CORE, &rlim);
  //cout << " current core dump size limit: " << rlim.rlim_cur << endl;
  //cout << " maximum core dump size limit: " << rlim.rlim_max << endl;
  rlim.rlim_cur = rlim.rlim_max; //300000000;

  if(!argv[1]) usage();

  char** arg_ptr = argv;

  IniDocument* ini_doc=0;
  string dirname="";
  string olddirname="";
  input_meqn=0;

  //for(int arg = 1; arg < argc; arg++)
  while(*(++arg_ptr))
  {
    //printf("*arg_ptr=%s\n",*arg_ptr);
    // is the argument an option?
    if((*arg_ptr)[0] == '-')
    {
      if(!strcmp(*arg_ptr,"-meqn")
      || !strcmp(*arg_ptr,"-meqns"))
      {
         arg_ptr++;
         sscanf(*arg_ptr, "%d",&input_meqn);
         printf("setting meqn=%d\n",input_meqn);
      }
    }
    else // assume that it is a file
    {
      string dirname = getDirectoryOfFile(*arg_ptr);
      if(olddirname!=dirname)
      {
        delete ini_doc;
        ini_doc = new IniDocument;
        ini_doc->initFromFile(dirname+"/parameters.ini");
        dogParams.init(*ini_doc);
        dogParamsCart2.init(*ini_doc);
        // override with qhelp.dat values
        {
          int qhelp_mx;
          int qhelp_my;

          string qhelp_filename = dirname+"/qhelp.dat";
          FILE* qhelp = fopen(qhelp_filename.c_str(),"r");
          for(int i=0;i<5;i++)
          {
            char char_buffer[101];
            fscanf(qhelp,"%100s",char_buffer);
          }
          fscanf(qhelp,"%d",&qhelp_mx);
          fscanf(qhelp,"%d",&qhelp_my);
          if(qhelp_mx != dogParamsCart2.get_mx()){
            dogParamsCart2.set_mx(qhelp_mx);
            printf("reset mx to value in qhelp.dat: %d\n",qhelp_mx);
          }
          if(qhelp_my != dogParamsCart2.get_my()){
            dogParamsCart2.set_my(qhelp_my);
            printf("reset my to value in qhelp.dat: %d\n",qhelp_my);
          }
        }
        // default value of meqn
        if(input_meqn==0) meqn = dogParams.get_meqn();
        else meqn = input_meqn;
        olddirname=dirname;
      }

      dat2h5(*arg_ptr);
    }
  }
}
