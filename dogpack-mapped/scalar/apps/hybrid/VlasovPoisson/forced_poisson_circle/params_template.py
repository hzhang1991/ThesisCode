# parameters template. 
#
# This should agree with the current params.ini with the exception of the
# fields we wish to modify.
#
# This template is called by run_sample.py for producing many outputs that are
# needed to run a convergence study.
#

# template for running dog.exe
dogpack_data_template = '''
; Parameters common to dogpack applications
[dogParams]
defaults_file = "$DOGPACK/config/dogParams_defaults.ini"
ndims       = 4         ; 1 or 2
mesh_type   = Hybrid    ; (either Hybrid or Unstructured) 
nout        = 1          ; number of output times to print results
tfinal      = 1.0       ; final time
;dtv(1)      = %(dt)e      ; initial dt
;dtv(2)      = %(dt)e      ; max allowable dt 
dtv(1)      = 0.0001    ; initial dt              (First step will always be thrown away)
dtv(2)      = 1.0e5     ; max allowable dt        
cflv(1)     = 0.21      ; max allowable Courant number 1/(2k-1) (for local problem)
cflv(2)     = 0.20      ; desired Courant number                (for local problem)
nv          = 500000    ; max number of time steps per call to DogSolve
time_stepping_method = Runge-Kutta ; (e.g., Runge-Kutta, SDC, Lax-Wendroff, User-Defined)
limiter_method = moment ; (e.g., moment, viscosity)
space_order = %(s_order)i   ; =method(1)= order of accuracy in space
time_order  = %(t_order)i   ; =method(2)= order of accuracy in time
use_limiter = 0  ; =method(3)= use limiter (1-yes, 0-no)
verbosity   = 1  ; =method(4)= verbosity of output
mcapa       = 0  ; =method(5)= mcapa (capacity function index in aux arrays)
maux        = 0  ; =method(6)= maux (number of aux arrays, maux >= mcapa)
source_term = 1  ; =method(7)= source term (1-yes, 0-no)
meqn        = 1  ; number of equations
mrestart    = 0  ; restart from old data (1-yes, 0-no)
nstart      = 0  ; if mrestart==1: from file q(nstart)_restart.dat
datafmt     = 1  ; 1 for ascii, 5 for hdf5.
ic_quad_order = 20  ; order of quadrature for L2-projection of initial conditions
; -------------------------------------------------
;   Cartesian grid data (ignored if Unstructured)
; -------------------------------------------------
[grid]
mx    =  %(mx)i  ; number of grid elements in x-direction
my    =  %(my)i  ; number of grid elements in y-direction
mbc   = 0     ; number of ghost cells on each boundary  (periodicity hard coded)
xlow  = -6.0e0 ; left end point
xhigh =  6.0e0 ; right end point
ylow  =  0.0e0 ; lower end point
yhigh =  1.0e0 ; upper end point
'''

# template for running the unstructured code
input2d_template = '''
Unstructured          : Grid type (Structured vs. Unstructured)
%(h0)f                : Initial grid spacing (H0)
50000                 : Maximum allowed iterations
2                     : Create a submesh with given integer refinement factor
----------------------------------------------------------------
Bounding box:
 -1.0e0                : xmin
  1.0e0                : xmax
 -1.0e0                : ymin
  1.0e0                : ymax
----------------------------------------------------------------
4                     : Number of fixed points
----------------------------------------------------------------
Fixed points (x y):
 0.0e0   -1.0e0
 0.0e0    1.0e0
-1.0e0    0.0e0
 1.0e0    0.0e0
'''
