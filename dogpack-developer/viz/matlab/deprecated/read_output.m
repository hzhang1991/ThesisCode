% read sampled output
%
% dogpack no longer produces sampled output, so this method
% exists only for the purpose of displaying old data
%
function [qfull,time]=read_output(datafmt, basefilename, ...
  mx, my, meqn, method1, components, READ_ALL_AT_ONCE)

  if(nargin<8)
    READ_ALL_AT_ONCE=1;
  end
  % in matlab it seems less efficient to read in only the needed components
  assert(READ_ALL_AT_ONCE==1);
  
  if(datafmt==1) % ASCII
    [qfull,time]=read_output_ASCII(basefilename, ...
      mx, my, meqn, method1, components, READ_ALL_AT_ONCE);
  elseif(datafmt==5) % HDF5
    [qfull,time]=read_output_HDF5(basefilename, ...
      mx, my, meqn, method1, components, READ_ALL_AT_ONCE);
  else
    error(['unsupported value of datafmt: ' num2str(datafmt)]);
  end
end

function [qfull,time]=read_output_ASCII(basefilename, ...
  mx, my, meqn, method1, components, READ_ALL_AT_ONCE);
  %
  filename = [basefilename '.dat'];
  fids = fopen(filename,'r');
  time = fscanf(fids,'%e',1);
  % in the ascii case, for backwards compatibility,
  % because of the order of storage
  % we have to read all the data
  qsoln = fscanf(fids,'%e',[meqn,inf]);
  fclose(fids);
  % qsoln = transpose(qsoln);
  mx_out = mx*method1;
  my_out = my*method1;
  qfull=zeros(mx_out,my_out,numel(components));
  for me=components %1:meqn
    qfull(:,:,me) = reshape(qsoln(me,:),mx_out,my_out);
  end
end

function [qfull,time]=read_output_HDF5(basefilename, ...
  mx, my, meqn, method1, components, READ_ALL_AT_ONCE)
  %
  filename = [basefilename '.h5'];
  fileinfo=hdf5info(filename);
  %dims=fileinfo.GroupHierarchy.Datasets.Dims;
  time = fileinfo.GroupHierarchy.Datasets.Attributes.Value;
  READ_ALL=1;
  if(isequal(components,1:meqn)||READ_ALL)
    % does this automatically reverse the dimensions?
    qfull=hdf5read(filename,'/qvals');
    %size(qfull)  % check the dimensions
    if(~isequal(components,1:meqn))
      qfull=qfull(:,:,components);
    end
  else
    file = H5F.open(filename, 'H5F_ACC_RDONLY', 'H5P_DEFAULT');
    dataset = H5D.open(file, 'qvals');
    dataspace = H5D.get_space(dataset);
    [rank, dims] = H5S.get_simple_extent_dims(dataspace);
    assert(rank==3);
    mx_out = mx*method1;
    my_out = my*method1;
    expected_dims=[mx_out,my_out,meqn];
    assert(isequal(dims',fliplr(expected_dims)));
    %
    % Check that the data type is double.
    % (this check might fail if reading a file produced
    % on a different architecture)
    %
    datatype=H5D.get_type(dataset);
    datatype_expected=H5T.copy('H5T_NATIVE_DOUBLE');
    assert(logical(H5T.equal(datatype,datatype_expected)));
    assert(logical(H5T.equal(datatype,'H5T_NATIVE_DOUBLE')));
    %
    % just select the desired componenets
    count=[mx_out,my_out,1];
    dims_requested = [mx_out, my_out, numel(components)];
    [sorted_components,perm] = sort(components);
    components_are_sorted=(isequal(components,sorted_components));
    %
    %
    % Why does this succeed versus [read_all_tag] above does not?
    READ_ALL_AT_ONCE=1;
    if (components_are_sorted || READ_ALL_AT_ONCE)
      selection_operator='H5S_SELECT_SET';
      for i=1:numel(components)
        offset=[0,0,components(i)-1];
        H5S.select_hyperslab(dataspace,selection_operator,...
          fliplr(offset),[],fliplr(count),[]);
        selection_operator='H5S_SELECT_OR';
      end
      memspace = H5S.create_simple(rank,fliplr(dims_requested),[]);
      H5S.select_hyperslab(memspace, 'H5S_SELECT_SET', [0,0,0], [],...
        fliplr(dims_requested), []);
      %
      qfull=H5D.read(dataset, datatype, memspace, dataspace,...
       'H5P_DEFAULT');
      if(~isequal(components,sorted_components))
        %iperm(perm(i:numel(perm))) = i:numel(perm);
        %qfull = qfull(:,:,iperm);
        qfull(:,:,perm) = qfull;
      end
      H5S.close(memspace);
    else
      assert(READ_ALL_AT_ONCE==0);
      % read in one component at a time
      qfull=zeros(mx_out,my_out,numel(components));
      for i=1:numel(components)
        offset=[0,0,components(i)-1];
        H5S.select_hyperslab(dataspace,'H5S_SELECT_SET',...
          fliplr(offset),[],fliplr(count),[]);
        memspace = H5S.create_simple(2,fliplr([mx_out,my_out]),[]);
        H5S.select_hyperslab(memspace, 'H5S_SELECT_SET', [0,0], [],...
          fliplr([mx_out,my_out]), []);
        %
        qfull(:,:,i)=H5D.read(dataset, datatype, memspace, dataspace,...
          'H5P_DEFAULT');
        H5S.close(memspace);
      end
    end
    H5D.close(dataset);
    H5S.close(dataspace);
    H5F.close(file);
  end
end
