# ============================================================================
#  Distributed under the terms of the Berkeley Software Distribution (BSD) 
#  license
#                     http://www.opensource.org/licenses/
# ============================================================================

import sys
import os
import getopt

def write_environment_variable(csh_handle,bash_handle,var,value):
    csh_handle.write('setenv %s "%s"\n' % (var.upper(),value))
    bash_handle.write('export %s="%s"\n' % (var.upper(),value))

def write_shortcut_alias(csh_handle,bash_handle,var,value):
    csh_handle.write( 'alias %s "%s"\n' % (var,value))
    bash_handle.write('alias %s="%s"\n' % (var,value))


# Paths
dogpack_path = os.path.abspath(os.curdir)
outfile_base="setenv"
    
# Open output files
csh_file = open(os.path.join(dogpack_path,".".join((outfile_base,"csh"))),'w')
bash_file = open(os.path.join(dogpack_path,".".join((outfile_base,"bash"))),'w')

# Write out boiler plate
boiler_plate = ("# DoGPack environment settings\n")
csh_file.write(boiler_plate)
bash_file.write(boiler_plate)

# Write out variables
python_path = "${PYTHONPATH}"
matlab_path = "${MATLABPATH}"

dogpackP = "".join((dogpack_path,"/viz/python"))
dogpackM = "".join((dogpack_path,"/viz/matlab"))

python_path = ":".join((dogpackP,python_path))
matlab_path = ":".join((dogpackM,matlab_path))

print ""
print "The following variables will be set:"
print "  DOGPACK = %s" % dogpack_path
write_environment_variable(csh_file,bash_file,"DOGPACK",dogpack_path)

print "  PYTHONPATH = %s" % python_path
write_environment_variable(csh_file,bash_file,"PYTHONPATH",python_path)

print "  MATLABPATH = %s" % matlab_path
write_environment_variable(csh_file,bash_file,"MATLABPATH",matlab_path)

pdog1_command = "python $DOGPACK/viz/python/plotdog1.py"
print "  plotdog1 = %s" %pdog1_command
write_shortcut_alias(csh_file,bash_file,"plotdog1",pdog1_command)

pdog2_command = "python $DOGPACK/viz/python/plotdog2.py"
print "  plotdog2 = %s" %pdog2_command
write_shortcut_alias(csh_file,bash_file,"plotdog2",pdog2_command)
print ""


# Close output files
csh_file.close()
bash_file.close()
