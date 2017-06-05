
% cyclically query user for variable and time to display
%
function plotout(c)

  %if(~exist('do_rescale','var')); do_rescale=1; end
  s = c.s;
  plot_params = s.plot_params;
  params = s.params;
  plot_grid = s.plot_grid;
  
  varMap = s.varMap;
  for i=1:varMap.num_variables
    % deprecate? used to interpret user input
    eval([char(varMap.name_arr(i)) '=' num2str(i) ';'])
  end
  %
  % so that user can enter q to quit or qq to completely quit
  q=-1;
  qq=-2;

  mf_display=0;
  while(mf_display~=-1)
    % display choices
    for i=1:varMap.num_variables
      disp(sprintf('%2d = %-15s = %s', i, char(varMap.name_arr(i)), ...
        plot_params.display_names{i}));
      %disp([num2str(i), ' = ', char(varMap.name_arr(i))]);
    end
    
    %default_input_string='[E,Eck,Bxu;Hal,dtJ,pek]''';
    default_input_string='[E,pek;dtui,Bxu]''';
    input_string = input([ 'Which component[s] of q do you want to plot ' ...
            default_input_string ' ) ? '],'s');
    disp(' ')
    if(isempty(input_string)||isempty(mf_display))
      % default components
      %mf_display=[E,Eck,Bxu;Hal,dtJ,pek]';
      mf_display = eval(default_input_string)
    else
      mf_display = eval([input_string ';']);
    end
    %
    if(length(mf_display)==1 && (mf_display==q || mf_display==qq))
      break;
    end
    %
    % restrict to selected subset of arrays
    %
    plot_struct = struct;
    plot_struct.name_arr = varMap.name_arr(mf_display);
    plot_struct.mf_display = mf_display;
    
    % create/raise the plot windows to be used
    plot_struct.use_subplots=0;
    pane_num=0;
    if(size(mf_display,1)>1)
      plot_struct.use_subplots=1;
    end
    if(plot_struct.use_subplots)
    % make array of subplots
      pane_num = varMap.num_variables+1;
      %disp(['figure(' num2str(pane_num) ')']);
      figure(pane_num);
      clf;
    else
      for i=1:numel(mf_display)
        %disp(['figure(mf_display(' num2str(i) ')']);
        figure(mf_display(i));
        % leave the contents of the figure untouched
        % in case the user wants to use this feature
        % to make a comparison.
        %clf;
      end
    end

    nf = 0;
    n1 = -1;
    while(nf~=-1)
      nf  = input([ ' Plot which frame ( 0 - ',num2str(params.nout),...
            ' ) [type -1 or q to quit] ? ']);
      if (nf==q)
        break
      end
      if(nf==qq)
        return
      end
      if isempty(nf)
        nf = n1 + 1;
      end
      if max(nf)> params.nout
        disp(' ');
        disp(' End of plots ');
        disp(' ');
        nf = params.nout;
      end

      % in case the user entered something like 0:5:40 (every fifth frame)
      next_frame=0;
      for midx=1:length(nf)
        n1=nf(midx);
        if(next_frame==0)
          next_frame=1;
        else
          go_next_frame  = input([  ' Plot frame ' num2str(n1) ' (or q)?']);
          if (go_next_frame==qq)
            return;
          elseif (go_next_frame==q)
            break
          end
        end
        [state,time] = read_state(c,n1);
        for i=1:numel(mf_display)
           if(mf_display(i)==0)
              %subplot(size(mf_display,1),size(mf_display,2),i);
              %clf; % need something that will only clear the subplot
           else

              disp(['calculating term ' plot_struct.name_arr{i}]);
              tic;
              var_name = s.varMap.name_arr{mf_display(i)};
              var_vals  = sample_var_from_state(c, var_name, state);
              display_panel_frame(c, i, time, var_vals, plot_struct);

              if(0)
                filename = ['/Users/evanjohnson/talk/figures/B' ...
                  num2str(n1) '_' num2str(params.mx) 'x' num2str(params.my) ...
                    '_ion_iso_period=' num2str(params.ion_iso_period)]
                filename = regexprep(filename, '\.','_');
                disp(['print to file: ' filename]);
                print('-dpng',filename);
              end
           end
        end
      end
    end
  end
  disp(' ')
end

function display_panel_frame(c, i, time, var_vals, plot_struct)
  plot_grid = c.s.plot_grid;
  params = c.s.params;
  % display the data
  disp(['displaying data']);
  if(plot_struct.use_subplots)
    disp(['selecting subplot ' num2str(i)]);
    subplot(size(plot_struct.mf_display,1),size(plot_struct.mf_display,2),i);
    display_fig=0;
  else
    display_fig=plot_struct.mf_display(i);
  end
  %
  var_name = char(plot_struct.name_arr(i));
  display_frame(c, time, var_vals, var_name, 'fig', display_fig);
end

