% names is a cell array of names or it is the number of components.
% if it is the number of components, then default names are chosen:
% q1, q2, q3, ...
%
function c=ComponentNames(names)

  %
  %if(isnumeric(names))
  %  assert(numel(names)==1);
  %  nameArray = cell(1,names);
  %  for i=1:names
  %    nameArray{i} = ['q' num2str(i)];
  %  end
  %  names = nameArray
  %end

  assert(iscell(names));

  s = struct;
  s.nameArray = names;
  s.name2idx = struct;
  for i=1:numel(names)
    s.name2idx.(names{i}) = i;
  end
  c = struct;
  c.s = s;
  c = class(c,'ComponentNames');
end
