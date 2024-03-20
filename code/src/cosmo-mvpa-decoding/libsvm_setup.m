this_dir = pwd;
%make sure libsvm repo is cloned 
%https://github.com/cjlin1/libsvm
main_github_path = '/Users/jacekmatuszewski/Documents/GitHub';

cd main_github_path
cd libsvm; % change this to the directory where you put LIBSVM
cd matlab  % go to matlab sub-directory
make       % compile libsvm mex functions; requires a working compiler
rmpath(pwd)   % } ensure directory is on top
addpath(pwd)  % } of the search path

% verify it worked.
cosmo_check_external('libsvm'); % should not give an error

cd this_dir