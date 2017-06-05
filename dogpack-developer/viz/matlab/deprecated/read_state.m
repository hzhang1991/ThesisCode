function [q,time]=read_state(datafmt, basefilename, ...
  mx, my, meqn, method1, components)

  disp('read_state: reading data');
  tic;
  if(datafmt==1) % ASCII
    [q,time]=read_state_ASCII(basefilename, ...
      mx, my, meqn, method1, components);
  elseif(datafmt==5) % HDF5
    [q,time]=read_state_HDF5(basefilename, ...
      mx, my, meqn, method1, components);
  else
    error(['unsupported value of datafmt: ' num2str(datafmt)]);
  end
  toc;
end

function [q,time]=read_state_ASCII(basefilename, ...
    mx, my, meqn, method1, components)
  %
  filename = [basefilename '.dat']; % was '_restart.dat'
  disp(['reading file ' filename]);
  fids = fopen(filename,'r');
  time = fscanf(fids,'%e',1);
  qsoln = fscanf(fids,'%e',[inf]);
  fclose(fids);
  % order of subscripts must match WriteRestartASCII in Output.cpp
  % (but reversed because matlab unlike C uses reverse odometer subscripts)
  dims = [get_kmax(method1),meqn,mx,my];
  q=reshape(qsoln,dims);
  % shift the subscripts to match HDF5 format
  % q=permute(q,[3,4,1,2]); % same as ipermute in this case
  q=q(:,components,:,:);
  q=shiftdim(q,2);
  %q=q(:,:,:,components);
  dims=circshift(dims,[0,2]);
end

function [q,time]=read_state_HDF5(basefilename, ...
  mx, my, meqn, method1, components)
  %
  kmax=get_kmax(method1);
  filename = [basefilename '.h5'];  % was '_restart.h5'
  disp(['reading file ' filename]);
  fileinfo=hdf5info(filename);
  expected_dims=[mx,my,kmax,meqn];
  dims=fileinfo.GroupHierarchy.Datasets.Dims;
  assert(isequal(dims,expected_dims));
  time = fileinfo.GroupHierarchy.Datasets.Attributes.Value;
  READ_ALL=1;
  % it is apparently more efficient to read in the entire data;
  if(isequal(components,1:meqn)||READ_ALL)
    q=hdf5read(filename,'/q');
    [d1,d2,d3,d4]=size(q);
    assert(all(dims==[d1,d2,d3,d4]));
    if(~isequal(components,1:meqn))
      q=q(:,:,:,components);
    end
  else
  % why is this stuff so slow?  is it interpreted?`
    file = H5F.open(filename, 'H5F_ACC_RDONLY', 'H5P_DEFAULT');
    dataset = H5D.open(file, 'q');
    dataspace = H5D.get_space(dataset);
    [rank, dims] = H5S.get_simple_extent_dims(dataspace);
    assert(rank==4);
    expected_dims=[mx,my,kmax,meqn];
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
    count=[mx,my,kmax,1];
    mx_out=mx*method1;
    my_out=my*method1;
    dims_requested = [mx_out, my_out, kmax, numel(components)];
    [sorted_components,perm] = sort(components);
    components_are_sorted=(isequal(components,sorted_components));
    %
    if(0)
    % Why does the following fail while [read_all_tag] below does not?
    %if(components_are_sorted || READ_ALL_AT_ONCE)
      selection_operator='H5S_SELECT_SET';
      for i=1:numel(components)
        offset=[0,0,0,components(i)-1];
        H5S.select_hyperslab(dataspace,selection_operator,...
          fliplr(offset),[],fliplr(count),[]);
        selection_operator='H5S_SELECT_OR';
      end
      memspace = H5S.create_simple(rank,fliplr(dims_requested),[]);
      H5S.select_hyperslab(memspace, 'H5S_SELECT_SET', [0,0,0,0], [],...
        fliplr(dims_requested), []);
      % Fails on this line
      q=    H5D.read(dataset, datatype, memspace, dataspace,...
       'H5P_DEFAULT');
      if(~isequal(components,sorted_components))
        %iperm(perm(i:numel(perm))) = i:numel(perm);
        %qfull = qfull(:,:,:,iperm);
        q(:,:,:,perm) = q;
      end
      H5S.close(memspace);
    else
      %assert(READ_ALL_AT_ONCE==0);
      % read components one at a time
      q=zeros(mx,my,kmax,numel(components));
      for i=1:numel(components)
        offset=[0,0,0,components(i)-1];
        H5S.select_hyperslab(dataspace,'H5S_SELECT_SET',...
          fliplr(offset),[],fliplr(count),[]);
        memspace = H5S.create_simple(3,fliplr([mx,my,kmax]),[]);
        H5S.select_hyperslab(memspace, 'H5S_SELECT_SET', [0,0,0], [],...
          fliplr([mx,my,kmax]), []);
        q(:,:,:,i)=H5D.read(dataset, datatype, memspace, dataspace,...
          'H5P_DEFAULT');
        H5S.close(memspace);
      end
    end
    H5D.close(dataset);
    H5S.close(dataspace);
    H5F.close(file);
  end
end
