function parameters = read_more_params_from_ini(ini_fileName, parameters, varargin)
%=========================================================================
% Define variables from .ini file ignoring section labels
%
% ini_fileName: ini file (including .ini extension)
% parameters: struct
%==========================================================================

  try
      fid = fopen(ini_fileName);
  catch
      rethrow(lasterr);
  end

  if(fid==-1)
    warning(['no such file: ' ini_fileName]);
    return;
  end

  % hash table of option names
  hash = struct;
  % create a struct with values of NaN for all options listed
  for idx=1:numel(varargin)
      hash.(varargin{idx}) = 1;
  end
  while ~feof(fid)
    tline = fgets(fid);
    % delete any comments and whitespace
    %
    % strip away comment character followed by non-quote characters
    % (we will not allow a quote character in a comment)
    tline = regexprep(tline,'[%#;][^"]*', '');
    % strip away whitespace
    tline = strtrim(tline);
    if(~isempty(tline) ...
      && ~('#'==tline(1)) ...
      && ~(';'==tline(1)) ...
      && ~('%'==tline(1)) ...
      )
      % could also ignore anything that clearly does not match what we are looking for
      section = regexp(tline, {'\[[.\w]+\]'}, 'match');
      if(~isempty(section{1}))
        continue;
      end
      option = parse_ini_option(tline);
      % do not add the parameter unless requested
      %param = regexp(option.param,'(?<name>\w+)(\.\w+)?(?<subscript>\(\d+\)|\{\d+\})?$', 'names');
      param = regexp(option.param,'(?<name>\w+).*$', 'names');
      if(~isfield(hash,param.name))
        continue;
        %break;
      end;
      option_val = extract_option_value(option);
      command = ['parameters.' option.param '=' option_val ';'];
      eval(command);
    end
  end
  fclose(fid);
end

function option = parse_ini_option(tline)
    option = regexp(tline, '(?<param>\S+)\s*=\s*(?<val>.*)', 'names');
    assert(~isempty(option), ['could not parse: ' tline]);
end
