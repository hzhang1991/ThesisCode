function parameter_hierarchy = read_ini(ini_fileName)
  %=========================================================================
  % Reads .ini file and returns a hierarchical structure
  %
  % ini_fileName: ini file (including .ini extension)
  %
  % parameter_hierarchy: structure of section structures.
  %     every member field of a section structure is one
  %     of the following types:
  %     (*) a section structure (thus providing for hierarchy),
  %     (*) a cell array (usually of strings),
  %     (*) an array of numbers, or
  %     (*) an atom, whose value may be a string or number.
  %         
  %   As a usage example, consider the following ini file
  %   and corresponding test code:
  %
  %     === begin file example.ini ===
  %     ; blank lines and comments are allowed
  %
  %     [section1]
  %     number = 5
  %     string1 = unquoted_string
  %     b(1) = 3.12e5 ; the value must be numeric
  %     b(2) = 3      ; the value must be numeric
  %     [section1.subsection]
  %     number = 1.57
  %     [section2]
  %     c{1} = string1 ; a string
  %     c{2} = "string with spaces" ; a quoted string
  %     c{3} = "3.6" ; a string
  %     c{4} = 3.6   ; a number
  %     === end file example.ini ===
  %
  %     === matlab test code illustrating example ===
  %     % apologies for matlab's nonstandard strcmp function
  %     options=read_ini('example.ini')
  %     assert(       options.section1.number == 5)
  %     assert(strcmp(options.section1.string1, 'unquoted_string'))
  %     assert(       options.section1.b(1) == 3.12e5)
  %     assert(       options.section1.b(2) == 3)
  %     assert(       options.section1.subsection.number == 1.57)
  %     assert(strcmp(options.section2.c{1} ,  'string1'))
  %     assert(strcmp(options.section2.c{2} ,  'string with spaces'))
  %     assert(strcmp(options.section2.c{3} ,  '3.6'))
  %     assert(       options.section2.c{4} == 3.6)
  %
  %==========================================================================
  parameter_hierarchy = []; 
  
  try
      fid = fopen(ini_fileName);
  catch
      rethrow(lasterr);
  end
  MATCH_SECTION = false;
  CURR_SECTION = '';
  paramAndValueArray = [];
  j = 1;
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
          section = regexp(tline, {'\[[.\w]+\]'}, 'match');
          if ~isempty(section{1})
              section = regexprep(section{1}{1},{'\[','\]'},'');
              MATCH_SECTION = true;
              if ~isempty(CURR_SECTION)
                  % done wih previous section so save it
                  command=['parameter_hierarchy.' CURR_SECTION '=paramAndValueArray;'];
                  eval(command);
                  %name_hierarchy = regexp(CURR_SECTION, '(\w+)', 'tokens');
                  %parameter_hierarchy.(CURR_SECTION) = paramAndValueArray;
                  j = 1;
                  paramAndValueArray = [];
              end
              CURR_SECTION = section;
          else
              MATCH_SECTION = false;
          end
          if ~MATCH_SECTION
              option = parse_ini_option(tline);
              option_val = extract_option_value(option);
              command = ['paramAndValueArray.' option.param '=' option_val ';'];
              eval(command);
              j=j+1;
          end
      end
  end
  %parameter_hierarchy.(CURR_SECTION) = paramAndValueArray;
  command=['parameter_hierarchy.' CURR_SECTION '=paramAndValueArray;'];
  eval(command);
  fclose(fid);
end % END iniread
