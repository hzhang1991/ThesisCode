This is an example of solving Vlasov on an annulus using the splitting
technique.  To run this example, you must do the following:

1.) Generate the mesh:

cd Unstructured\_Mesh
make; ./mesh.exe; cd ..

In order to accomodate the Poisson solver, you need to set the "sub-factor"
parameter equal to the order of the method.  This solver accomodates orders
1,2, and 3.

2.) Use MATLAB to create the matrix used in the Poisson solve:

Open MATLAB and navigate to /path/to/this/app/matlab

type

CreateMatrix( order ),

where order = same order used in creating the sub-mesh.

3.) Compile and run the code.

After creating the mesh and the matrix that will be used for the Poisson
solves, you may now run the code.  Make sure that the spatial order is the
same as that set in the submesh, which is the same as that set in
CreateMatrix.

