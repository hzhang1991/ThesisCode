function CreateSliceList( )
%CREATESLICELIST  Find and save the list of slices along y=0.
%
% This function finds all the elements that intersect the line y=0 and saves
% them in the outputdirectory.

meshdir   = ['../Unstructured_Mesh/mesh_output'];
[BndyList, Psort, Mid, Left, Right] = FindSliceIndexUnst2( meshdir, meshdir );

end
