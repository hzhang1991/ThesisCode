#include <string>
#include "dogdefs.h"

using std::string;

#if(COMPILE_WITH_HDF5!=0)
#include "hdf5.h"

// routines to allocate and free C-style arrays.

static void freeArray(int rank, void* array)
{
    if(rank<1)
    {
        eprintf("invalid array rank: %d",rank);
    }

    if(rank>1)
    {
        // first free the array that the first element points to
        freeArray(rank-1,((void**)array)[0]);
    }
    // printf("freeing memory address %x\n",array);
    free(array);
}

static void* allocateArray(
    int rank,
    hsize_t* dims,
    size_t element_size,
    void** elementArr=NULL, // a reference to array of elements
    int* num_elements=NULL  // a reference to the number of elements
    )
{
    void *ret_array=NULL;
    void* parent_array=NULL; // an array of pointers to elements

    if(rank<1)
    {
        printf("invalid rank: %d\n", rank);
        exit(1);
    }
    //
    // in case NULL pointers were passed
    //
    void *elementArr_alternate_referent;
    if(!elementArr) elementArr=&elementArr_alternate_referent;
    int num_elements_alternate_referent;
    if(!num_elements) num_elements=&num_elements_alternate_referent;

    // recursively allocate parent arrays
    // and determine the number of elements in the current array.
    //
    (*num_elements)=dims[rank-1];
    int num_elements_parent=0;
    if(rank>1)
    {
        ret_array=allocateArray(
            rank-1,
            dims,
            sizeof(void*),
            &parent_array,
            &num_elements_parent);
        (*num_elements)*=num_elements_parent;
    }

    // allocate the array of elements
    //
    size_t alloc_size=element_size*(*num_elements);
    (*elementArr)=calloc(1,alloc_size);
    // printf("allocated memory address %x\n",(*elementArr));
    if(rank==1)
    {
        ret_array = (*elementArr);
    }

    // make the parent array elements point to subarrays
    // in the array of elements
    //
    if(rank>1)
    {
        size_t charsize = sizeof(char); // probably 1.
        int char_indices_per_element = element_size/charsize;
        const int& num_elements_per_subarray = dims[rank-1];
        int char_indices_per_subarray
            = char_indices_per_element*num_elements_per_subarray;

        for(int i=0;i<num_elements_parent;i++)
        {
            ((char**) parent_array)[i]
                =&((char*)(*elementArr))[i*char_indices_per_subarray];
        }
    }

    return ret_array;
}

static void WriteSampleValuesHDF5(
    string fname, const dTensor3& q, double t)
{
    fname+=".h5";

    const int RANK=3;
    hsize_t     dims[RANK];            /* dataset dimensions */
    double***   data;                  /* data to write */
    hid_t       file, dataset;         /* file and dataset handles */
    hid_t       datatype, dataspace;   /* handles */
    hid_t       time_attr_id, time_attr;
    herr_t      status;
    int         i, j, m;

    // allocate memory for data
    //
    int num_elements;
    dims[0]=q.getsize(3);
    dims[1]=q.getsize(2);
    dims[2]=q.getsize(1);
    data = (double***)
        allocateArray(RANK,dims,sizeof(double),NULL,&num_elements);
    int num_elements_calc=1;
    for(int i=1;i<=3;i++) num_elements_calc*=q.getsize(i);
    assert(num_elements==num_elements_calc);
    //data  = (double*)calloc(alloc_size);
    if(data==NULL)
    {
        eprintf("HDF5.cpp: WriteOutputHDF5: could not allocate memory");
    }

    // populate data array with q
    //
    // HDF5 should have a way to transfer the data
    // directly from the tensor's vector,
    // but I have not figured it out yet.
    for (i=0; i<q.getsize(1); i++) // x-axis idx
    for (j=0; j<q.getsize(2); j++) // y-axis idx
    for (m=0; m<q.getsize(3); m++) // eqn index
    {
        data[m][j][i]=q.get(i+1,j+1,m+1);
        //printf("data[%d][%d][%d]=%f\n", m,j,i,data[m][j][i]);
    }

    // write data to file
    //
    file = H5Fcreate(fname.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    dataspace = H5Screate_simple(RANK, dims, NULL);
    datatype = H5Tcopy(H5T_NATIVE_DOUBLE);
    // change to put q and aux into one file?
    dataset = H5Dcreate2(file, "qvals", datatype, dataspace,
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE,
        H5S_ALL, H5S_ALL, H5P_DEFAULT, **data);
    // record the time as an attribute (metadata)
    time_attr_id  = H5Screate(H5S_SCALAR);
    time_attr = H5Acreate2(dataset, "t", H5T_NATIVE_DOUBLE,
                     time_attr_id, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Awrite(time_attr, H5T_NATIVE_DOUBLE, &t);

    // free resources
    //
    H5Sclose(time_attr_id);
    H5Sclose(dataspace);
    H5Aclose(time_attr);
    H5Tclose(datatype);
    H5Dclose(dataset);
    H5Fclose(file);
    freeArray(RANK,data);
}

// We use C rather than C++ bindings to avoid requiring the
// user to compile hdf5 with C++ support.
//
// We write the array as a one-dimensional array with
// a time attribute and an array attribute that specifies
// the dimensions.  (Alternatively we could pass a multidimensional
// C array, the construction of which which would require a call to
// malloc for every dimension.)
//
// Probably we should record both q and a in a single file.
//
void WriteStateHDF5(
    string fname,
    string varname,
    const dTensorBC4& q,
    double t)
{
    fname+=".h5";

    const int RANK=4;
    double****  data;
    hsize_t     dims[RANK];            /* dataset dimensions */
    hid_t       file, dataset;         /* file and dataset handles */
    hid_t       datatype, dataspace;   /* handles */
    hid_t       time_attr_id, time_attr;
    herr_t      status;
   
    // allocate memory for data
    //
    // define dimensions; must be consistent with indices below
    //
    int num_elements;
    dims[0]=q.getsize(3);
    dims[1]=q.getsize(4);
    dims[2]=q.getsize(2);
    dims[3]=q.getsize(1);
    data = (double****)
        allocateArray(RANK,dims,sizeof(double),NULL,&num_elements);
    int num_elements_calc=1;
    for(int i=1;i<=4;i++) num_elements_calc*=q.getsize(i);
    assert(num_elements==num_elements_calc);
    //data  = (double*)calloc(alloc_size);
    if(data==NULL)
    {
        eprintf("HDF5.cpp: WriteStateHDF5: could not allocate memory");
    }
    
    // populate data array with q
    //
    // I believe that HDF5 has a way to transfer the data
    // directly from the tensor's vector using the HDF5 notion
    // of a hyperslab.
    // Copying the data allows us to store the data on disk
    // in a different order than in memory.
    // If the eqn idx changes most slowly (i.e. is the first index of
    // data) then if the user only wants to access certain
    // components fewer disk inodes need be accessed.
    // But if it so happens that the user only wants to access
    // part of the domain, then a good strategy might be to make
    // mx and my change most slowly.
    //
    for(int nx=0;  nx<q.getsize(1);  nx++) // x-axis
    for(int ny=0;  ny<q.getsize(2);  ny++) // y-axis
    for(int neq=0; neq<q.getsize(3); neq++) // eqn idx
    for(int nk=0;  nk<q.getsize(4);  nk++) // coef idx
    {
        // must be consistent with indices above
        //int idx=n4+q.getsize(4)*(n3+q.getsize(3)*(n2+q.getsize(2)*n1));
        //data[idx]=q.get(n1+1,n2+1,n3+1,n4+1);
        //data[nx][ny][neq][nk]=q.get(nx+1,ny+1,neq+1,nk+1);
        data[neq][nk][ny][nx]=q.get(nx+1,ny+1,neq+1,nk+1);
    }

    // write data to file
    //
    file = H5Fcreate(fname.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    dataspace = H5Screate_simple(RANK, dims, NULL);
    datatype = H5Tcopy(H5T_NATIVE_DOUBLE);
    // change to put q and aux into one file?
    dataset = H5Dcreate2(file, varname.c_str(), datatype, dataspace,
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE,
        H5S_ALL, H5S_ALL, H5P_DEFAULT, ***data);
    // record the time as an attribute (metadata)
    time_attr_id  = H5Screate(H5S_SCALAR);
    time_attr = H5Acreate2(dataset, "t", H5T_NATIVE_DOUBLE,
                     time_attr_id, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Awrite(time_attr, H5T_NATIVE_DOUBLE, &t);

    // free resources
    //
    H5Sclose(time_attr_id);
    H5Sclose(dataspace);
    H5Aclose(time_attr);
    H5Tclose(datatype);
    H5Dclose(dataset);
    H5Fclose(file);
    freeArray(RANK,data);
}

double ReadStateArrayHDF5(string fname, string varname, dTensorBC4& q)
{
    fname += ".h5";
    hid_t file = H5Fopen(fname.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if(file<0)
    {
        eprintf("HDF5.cpp: unable to open file %s", fname.c_str());
    }
    hid_t dataset = H5Dopen2(file, varname.c_str(), H5P_DEFAULT);
    //
    // get time attribute
    //
   hid_t attr = H5Aopen(dataset, "t", H5P_DEFAULT);
   double time_out;
   H5Aread(attr, H5T_NATIVE_DOUBLE, &time_out);
   H5Aclose(attr);
    //
    // check data type
    //
    hid_t datatypeHandle  = H5Dget_type(dataset);
    H5T_class_t datatypeClass = H5Tget_class(datatypeHandle);
    assert (datatypeClass == H5T_FLOAT);
    H5T_order_t order = H5Tget_order(datatypeHandle);
    assert(order==H5T_ORDER_LE);
    //
    hid_t dataspaceHandle = H5Dget_space(dataset);
    //
    // check dimensions
    //
    int rank = H5Sget_simple_extent_ndims(dataspaceHandle);
    assert(rank==4);
    hsize_t datasetDims[4];
    H5Sget_simple_extent_dims(dataspaceHandle, datasetDims, NULL);
    assert(datasetDims[0]==q.getsize(3));
    assert(datasetDims[1]==q.getsize(4));
    assert(datasetDims[2]==q.getsize(2));
    assert(datasetDims[3]==q.getsize(1));

    // Define hyperslab in the dataset.
    //
    hsize_t nullOffset[4];
    for(int i=0;i<rank;i++) nullOffset[i] = 0;
    H5Sselect_hyperslab(dataspaceHandle, H5S_SELECT_SET, nullOffset, NULL,
        datasetDims, NULL);

    // Define the memory dataspace.
    hid_t memspace = H5Screate_simple(rank,datasetDims,NULL);
    // Define memory hyperslab.
    H5Sselect_hyperslab(memspace, H5S_SELECT_SET, nullOffset, NULL,
       datasetDims, NULL);

    // Read data from hyperslab in the file into the hyperslab in memory.
    double *data_out_vec;
    double ****data_out = (double****)
        allocateArray(rank,datasetDims, sizeof(double),
            (void**)&data_out_vec,NULL);
    H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspaceHandle,
        H5P_DEFAULT, data_out_vec);

    // transfer from data array to data tensor
    for(int nx=0;  nx<q.getsize(1);  nx++) // x-axis
    for(int ny=0;  ny<q.getsize(2);  ny++) // y-axis
    for(int neq=0; neq<q.getsize(3); neq++) // eqn idx
    for(int nk=0;  nk<q.getsize(4);  nk++) // coef idx
    {
        // must be consistent with indices above
        //q.set(nx+1,ny+1,neq+1,nk+1, data_out[nx][ny][neq][nk]);
        q.set(nx+1,ny+1,neq+1,nk+1, data_out[neq][nk][ny][nx]);
    }

    freeArray(rank, data_out);
    H5Sclose(memspace);
    H5Sclose(dataspaceHandle);
    H5Tclose(datatypeHandle);
    H5Dclose(dataset);
    H5Fclose(file);
    return time_out;
}
#else // COMPILE_WITH_HDF5
#include <stdlib.h> // for getenv()
static void hdf5_error()
{
    const char* DOGPACK = getenv("DOGPACK");
    if(!DOGPACK) DOGPACK="${DOGPACK}";
    eprintf(
        "\n  COMPILE_WITH_HDF5 was not set,"
        "\n  yet datafmt=5 (HDF5) in dogpack.data."
        "\n  If you wish to use HDF5 format, you must"
        "\n  (1) make sure that you have HDF5 installed, and"
        "\n  (2) edit %s/config/Makefile.config"
        "\n      (You need to set COMPILE_WITH_HDF5=1), and"
        "\n  (3) rm %s/lib/2d/cart/HDF5.o and recompile."
        , DOGPACK, DOGPACK);
}
void WriteStateHDF5(
    string fname,
    string varname,
    const dTensorBC4& q,
    double t)
{
  hdf5_error();
}
void ReadStateArrayHDF5(string fname, string varname, dTensorBC4& q)
{
  hdf5_error();
}
#endif //COMPILE_WITH_HDF5
