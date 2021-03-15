% -------------------- %
% Setup File for GPOPS %
% -------------------- %
% Notes: 
%   (1) This file assumes that you have write permissions
%       to change the MATLAB path.  If you do cannot change the
%       MATLAB path, contact your system administrator.
%   (2) SNOPT is NOT included in the GPOPS distribution and must be
%       obtained and installed separately (see the GPOPS manual for
%       more information).
%   (2) Built-in automatic differentiation is included in the GPOPS
%       distribution and is installed by default using this setup
%       file.
%   (3) MAD and INTLAB are NOT included in the GPOPS distribution
%       and must be obtained separately (see the GPOPS manual for
%       more information).
currdir = pwd;
addir = strcat(currdir,'/ad/');
addpath(addir,path,'-begin');
libdir = strcat(currdir,'/lib/');
addpath(libdir,path,'-begin');
nlpdir = strcat(currdir,'/nlp/');
addpath(nlpdir,path,'-begin');
savepath;
