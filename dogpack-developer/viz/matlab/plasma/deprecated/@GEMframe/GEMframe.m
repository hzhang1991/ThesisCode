
% constructor for GEMframe (GEM Vframe class).
% It extends the Vframe class.
function frame=GEMframe(varargin)

  frame=struct;
  % pass all arguments through to the parent constructor
  vframe = Vframe(varargin{:});
  % inherit from Vframe
  frame = class(frame, 'GEMframe', vframe);
end
