%HYBRIDMT_INSTALL Add current directory to MATLAB path and build searcheable database.
rmpath(pwd);
addpath(pwd);
builddocsearchdb([pwd '\html']);
