function INI = ConvertIniFile2Struct(IniFileName)
%CONVERTINIFILE2STRUCT.    Parse a single ".ini" file.
%
% INI = CONVERTINIFILE2STRUCT( IniFileName ) converts ".ini" file IniFileName
% to a struct INI that contains a list of key-value pairs for each section
% found in IniFileName.
%
% Lines starting with ';' and '#' are ignored (comments).  Inline comments
% starting at the end of a valid line with ';' are also ignored.
%
% Parameters
% ----------
% 
%       IniFileName - A string identifying the ".ini" file to be read
%
% Returns
% -------
%
%       INI - A structure containing all key-value pairs from the ".ini" file.
%             These can be accessed through INI.[SECTION].[KEY].  All names
%             are converted to lower-case before saving them.
% 
% Author: David Seal, Michigan State University, Fall 2014

% open .ini file for reading.  Return error upon failure
[fid,msg] = fopen(IniFileName,'r');                
if( fid == -1 )
    error(['Could not open IniFileName', IniFileName]);
end

% Struct containing all field variable names
INI   = struct();            

% Read in entire file
while( ~feof(fid) )

    % Scan the next line, and remove any leading/trailing whitespace
    s = strtrim(fgetl(fid));    

    if isempty(s)
        % No information on this line, read next row
        continue;
    end

    % ';' and '#' start comment lines
    if(s(1)==';' || s(1)=='#' )
        continue;
    end

    if ( s(1)=='[' ) && (s(end)==']' )
        % This is a section title.  Create a second struct to hold its
        % information
        SectionName = genvarname( lower( strtrim( s(2:end-1) ) ) );
        INI.(SectionName) = struct();
    else

        % Add in a new field name to the current section
        [key, val] = parse_line(s);

        if( ~isempty(SectionName) )
            % Add this key/val pair to current section []
            INI.(SectionName).(lower(genvarname(key))) = val;
        else
            % Orphan value
            INI.(lower(genvarname(key))) = val;
        end
    end

end
fclose(fid);

end

function [key, val] = parse_line( line_name )

    % Add in a new field name to the current section
    [key, val] = strtok(line_name, '=');    % parse before and after the '=' sign
    key        = strtok(key);               % Remove extra whitespace

    % Clean up any trailing comments in the value
    %
    % TODO - include an option for killing the "#" symbol as well ...
    val        = strtok(val(2:end),';');
    val        = strtrim(val);

end
