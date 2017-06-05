
% read a number from the next line in the
% file referenced by fid and check that the string
% var_name is to its right

function num=get_num(fid,var_name)
  tline=fgetl(fid);
  if(~ischar(tline))
     error([var_name ': no line to read']);
  end;
  a1=textscan(tline,'%[^%=]',1);
  sz_a1_1=size(a1{1});
  if(sz_a1_1(1)==0)
     error([var_name ': cannot parse ' tline]);
  end;
  a2=char(a1{1}(1));
  % now split a2 using whitespace.
  a3=textscan(a2,'%[^ 	]',20);
  num = str2num(char(a3{1}(1)));
  if(isempty(num))
     error([var_name ': cannot parse num in ' tline]);
  end;
  if(nargin==1)
     return;
  end;
  % do checking if requested
  s3=size(a3{1});
  if(s3(1)~=2)
     error([var_name ': cannot parse ' tline ': expecting two tokens']);
  end;
  var_name_check = char(a3{1}(2));
  var_name_check(var_name_check=='[')='(';
  var_name_check(var_name_check==']')=')';
  %disp(['var_name_check=' var_name_check]);
  if(not(all(var_name_check==var_name)))
     error([var_name ': does not match ' var_name_check]);
  end
end
