function gpopsIntlabSetup(intlab_home_directory)

% -----------------------------------------------
% This function sets up INTLAB for Use with GPOPS
% -----------------------------------------------
% 
% Input: intlab_home_directory=Base Directory of INTLAB
% (e.g. Windows: c:\Intlab;  Unix: /home/user/Intlab)
dirs = {'intval','gradient','hessian','slope','utility','long'};
for i=1:length(dirs)
  currdir = fullfile(intlab_home_directory,dirs{i});
  addpath(currdir);
end;
savepath;
