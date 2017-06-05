def parse_grid_section( config, params ):
    """Parse the grid section of the parameters file."""

    params['xlow']    = config.getfloat ('grid', 'xlow'   )
    params['xhigh']   = config.getfloat ('grid', 'xhigh'  )
    params['mx']      = config.getint   ('grid', 'mx'     )
    params['mbc']     = config.getint   ('grid', 'mx'     )


#   def parameters_get_int( config, config_default, sec_name, key ):

#       from ConfigParser import NoOptionError

#       try:
#           val = config.getint(sec_name, key )    
#       except(NoOptionError):
#           val = config_default.getint( sec_name, key )
#       return val

def parse_ini_parameters(parameters_file, params={} ):
    """Parse a parameters.ini file.

    This function parses an input file parameters.ini and returns a dictionary
    containing key-value pairs for each required object in a DoGPack
    parameters.ini file.

    TODO - test for errors on checking these files - read in from the default
    parameters file if that object is missing.

    Input
    -----

        parameters_file - a string pointing to a valid parameters.ini file.

    Returns
    -------

        params - a dictionary containing key-value pairs.  If params already
        exists as a dictionary, this routine will add new keys.
    """

    import ConfigParser

    config = ConfigParser.RawConfigParser()
    config.read( parameters_file )

    #print(parameters_file)
    #print( config.sections() )

    # This file should be used to pull information from a section
#   defaults_file = config.get('dogParams', 'defaults_file' )
#   config_default = ConfigParser.RawConfigParser()
#   config_default.read( defaults_file )

#   print( defaults_file )
#   print(config_default.sections() )

#   params['ndims']       = parameters_get_int( config, config_default, 'dogParams', 'ndims'       )

    params['ndims']       = config.getint  ('dogParams', 'ndims'       )
    if params['ndims']==1:
        params['mesh_type'] = 'Cartesian'
    else:
        params['mesh_type']   = config.get     ('dogParams', 'mesh_type'   )
    params['nout']        = config.getint  ('dogParams', 'nout'        )
    params['tfinal']      = config.getfloat('dogParams', 'tfinal'      )
    params['dt_init']     = config.getfloat('dogParams', 'dtv(1)'      )
    params['dt_max']      = config.getfloat('dogParams', 'dtv(2)'      )
    params['cfl_max']     = config.getfloat('dogParams', 'cflv(1)'     )
    params['cfl_des']     = config.getfloat('dogParams', 'cflv(2)'     )
    params['nv']          = config.getint  ('dogParams', 'nv'          )

    params['time_stepping_method'] = config.get('dogParams', 'time_stepping_method' )
    params['limiter_method']       = config.get('dogParams', 'limiter_method'       )
    params['space_order'] = config.getint  ('dogParams', 'space_order' )
    params['time_order']  = config.getint  ('dogParams', 'time_order'  )
    params['use_limiter'] = config.getint  ('dogParams', 'use_limiter' ) 
    params['mcapa']       = config.getint  ('dogParams', 'mcapa'       )
    params['maux']        = config.getint  ('dogParams', 'maux'        )
    params['source_term'] = config.getint  ('dogParams', 'source_term' )
    params['meqn']        = config.getint  ('dogParams', 'meqn'        )
    #mrestart    = 0  ; restart from old data (1-yes, 0-no)
    #nstart      = 0  ; if mrestart==1: from file q(nstart)_restart.dat
    params['datafmt']     = config.getint ('dogParams', 'datafmt'      )

    if( config.has_section('grid') ):
        parse_grid_section( config, params )
    return params

