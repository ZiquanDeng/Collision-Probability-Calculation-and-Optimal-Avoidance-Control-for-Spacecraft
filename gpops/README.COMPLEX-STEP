----------------------------------------------
Use of Complex-Step Differentiation with GPOPS
----------------------------------------------

One of the choices in GPOPS for computation of the objective function
gradient and constraint Jacobian is complex-step differentiation.
Complex-step differentiation is perfectly well-defined for virtually
all functions, but some functions create problems.  In particular, the
built-in MATLAB functions ABS, MIN, MAX, and DOT do not work correctly
with complex-step differentiation.  Also, the standard transpose does
not work.  In the case of the transpose, it is necessary to use a
dot-transpose to ensure that a real transpose is taken (and not a
complex conjugate transpose).  Future releases of GPOPS will address
the issues with these limitations.


