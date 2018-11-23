%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Energy and Reserve Dispatch with Distributionally Robust Joint Chance Constraints
%   Christos Ordoudis, Viet Anh Nguyen, Daniel Kuhn, Pierre Pinson
% 
%   This script adds the folder src into the path of matlab
%   Notice that the src folder has to be in the same level as the Experiment1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

path1 = pwd;
path2 = strsplit(path1,filesep);
path3 = strjoin(path2(1:end-1),filesep);
path_src = strcat(path3,filesep,'src');
path(genpath(path_src),path); 
clear path1 path2 path3 path_src