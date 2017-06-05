
% read parameters from a file in the format;
%
%    1.2    variable_name_1  % comment 1
%    2      variable_name_2  % comment 2
%
%    str    variable_name_3  = comment 3
%
% For each variable name such as variable_name_1 statements like
%
%   global variable_name_1;
%   variable_name_1 = 1.2;
%
% are evaluated.
%
function read_config_file(readfile)

  %readfile = [outputdir '/param.data'];
  disp(['=== reading config file ' readfile ' ===']);
  fid=fopen(readfile,'r');
  if(fid==-1)
     error(['could not open file ' readfile]);
  end
  while 1
    tline=fgetl(fid);
    if ~ischar(tline), break, end
  
    % delete everything after comment character
    if(length(tline)==0)
       continue;
    end
    a1=textscan(tline,'%[^=%#]',1);
    sz_a1_1=size(a1{1});
    if(sz_a1_1(1)==0)
       continue;
    end
    a2=char(a1{1}(1));
    % change all commas to spaces,
    % all brackets to parentheses.
    s2=size(a2);
    for i=1:s2(1,2)
      if(a2(i)==',') a2(i)=' '; end
      if(a2(i)=='[') a2(i)='('; end
      if(a2(i)==']') a2(i)=')'; end
    end
    % now split a2 using spaces.
    a3=textscan(a2,'%[^ 	]',20);
    s3=size(a3{1});
    if((s3(1)>=2) && (mod(s3(1),2)==0)) 
      nvals=s3(1)/2;
      a4=reshape(a3{1},nvals,2);
      if(~isempty(str2num(char(a4(1,2)))))
      % if we can convert the first token
      % of the second half of the tokens to
      % a number, then do nothing
          ;
      else
      % assume that we have a sequence of names
      % of variables.
         for(i=1:nvals)
           %command0=['global ', char(a4(i,2)),';'];
           %eval(command0)
           % delete semicolon(s) if you want to see the command
           % or its output.
           varname = char(a4(i,2));
           if(isempty(str2num(char(a4(1,1)))))
             % if the first character is not a number,
             % then surround with quotes.
             command=[varname, '=''', char(a4(i,1)),''';'];
           else
             command=[varname, '=', char(a4(i,1)),';'];
           end
           % get a version of the variable name without any indexing
           a5=textscan(varname,'%[^(]',20);
           varbasename=char(a5{1}(1));
           %if(isempty(whos('global',varbasename)))
             global_command = ['global ' varbasename ';'];
             %disp(global_command);
             eval(global_command);
           %end
           disp(command);
           eval(command);
         end
      end
    end
  end
  fclose(fid);
end
