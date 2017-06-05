CTAGS_ARGS = \
  ctags --extra=+f \
        --langdef=file \
        --langmap='file:.ini.m.defs' \
        -R 

# the regular expressions below do not effectively match my functions
MATLAB_CTAGS_ARGS = \
  --file-tags \
  --langdef=matlab \
  --langmap=matlab:.m \
  '--regex-matlab=/^function[ \t]*([a-zA-Z0-9_]+)[ \t]*=[ \t]*([a-zA-Z0-9_]+)/\2/f,functions/' \
  '--regex-matlab=/^function[ \t]*([a-zA-Z0-9_]+)[^=]*$$/\1/f,functions/' \

retags:
	find `pwd` -name '*.m' | xargs ctags $(MATLAB_CTAGS_ARGS)

tags: retags
