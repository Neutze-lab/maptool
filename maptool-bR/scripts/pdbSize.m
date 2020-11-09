clear
close all

pdb = pdbread('/home/adams/data_proc/mapTool/maptool-bR-refined/input/5b6v-H_protein_ret_h2o.pdb');
margin = 5; 
gridSpacing = 0.25;

X = [pdb.Model.Atom(:).X]'; 
Y = [pdb.Model.Atom(:).Y]'; 
Z = [pdb.Model.Atom(:).Z]'; 

pdbLimits = [min(X) max(X) min(Y) max(Y) min(Z) max(Z)];

% Limits with border
limits = ceil([pdbLimits(1)-margin pdbLimits(2)+margin pdbLimits(3)-margin pdbLimits(4)+margin pdbLimits(5)-margin pdbLimits(6)+margin]);

CELLWORK = [limits(2)-limits(1) limits(4)-limits(3) limits(6)-limits(5) 90 90 90];

GRIDWORK = CELLWORK(1:3)*1/gridSpacing;

XYZLIM = limits/0.25;
%%% 
% bR
% CELL WORK 42 54 77 90 90 90
% GRID WORK 168 216 308
% XYZLIM -144 24 -176 40 -136 172
%
% Gives output margin 7
% for original pdb
% size =   -38.2480    8.8620  -46.5140   12.4340  -36.3970   44.5210

% for protein+retinal+waters+ 4 scaling waters margin 7
% size =   -38.2480    7.5060  -46.5140   12.0000  -36.3970   44.5210

% for protein+retinal+waters+ 4 scaling waters margin 5
% size =    -36.2480    5.5060  -44.5140   10.0000  -34.3970   42.5210

% Thus minimal size should be 2*38 2*46 2*44 -> 
% use at least CELL WORK 80 100 100
% 