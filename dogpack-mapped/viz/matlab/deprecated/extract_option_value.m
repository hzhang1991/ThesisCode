% helper (not called directly by user)

function option_val = extract_option_value(option)
    % strip off leading and trailing space and quotation marks if present
    option.val = strtrim(option.val);
    %
    % determine if the data is a string or a number
    num='';
    if(option.val(1)=='"')
      % should we handle escaped quote marks and escaped escapes?
      % would also want to add this functionality to the C++ parser.
      assert(option.val(end)=='"', ['missing ending quote ' ...
        'in definition of parameter ' option.param]);
      option.val=option.val(2:end-1);
    else
      num=str2num(option.val);
    end
    %
    if(isempty(num))
      option_val = ['''' regexprep(option.val, '''', '''''') ''''];
    else
      option_val = option.val;
    end
end

