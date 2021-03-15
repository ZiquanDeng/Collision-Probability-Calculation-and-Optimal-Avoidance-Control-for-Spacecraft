-------------------------------------------
Use of INTLAB for Automatic Differentiation
-------------------------------------------

One of the choices in GPOPS for computation of the objective function
gradient and constraint Jacobian is INTLAB (which stands for Interval
Laboratory).  INTLAB was created by Professor Siegfried M. Rump of
Hamburg University of Technology.  INTLAB is available at no charge
for internal or non-commercial use by visiting the following website:

    http://www.ti3.tu-harburg.de/rump/intlab/

All other uses of INTLAB require a license which must be obtained from
Dr. Siegfried M. Rump (e-mail:  rump@tu-harburg.de).

-----------------------------------------
Installation of INTLAB for Use with GPOPS
-----------------------------------------

While there is a setup routine in INTLAB, all of the directories in
INTLAB are NOT required when using GPOPS.  In fact, if you use the
default setup routine in INTLAB, you will end up with a filename
conflict due to the built-in MATLAB function 'legendre'.  As a result,
the following procedure should be following when installing INTLAB:

(1)  Download INTLAB from the above URL
(2)  Uncompress the archive (Windows or Unix, depending upon your OS)
(3)  run the file 'gpopsIntlabSetup.m' from within MATLAB.  Make sure
     when you run 'gpopsIntlabSetup.m' that you provide a string with
     the INTLAB home directory, e.g.,

              gpopsIntlabSetup('/home/user/Intlab/')
     where '/home/user/Intlab/' is the home directory of the INTLAB
     installation.  If successful, the following directories will be
     added to your MATLAB path:

     	/home/user/Intlab/intval
	/home/user/Intlab/gradient
	/home/user/Intlab/hessian
	/home/user/Intlab/slope
	/home/user/Intlab/utility
	/home/user/Intlab/long
	/home/user/Intlab/AccSum



