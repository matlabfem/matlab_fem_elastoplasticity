
% This script adds the absolute location of
% the shared functions to the path

path1 = cd;
%cd ..           
path2 = cd;
if ( isunix )
    addpath(genpath([path2,'.']));
else
    addpath(genpath([path2,'.']));
end
cd(path1)
%clear all