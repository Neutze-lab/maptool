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
%%
% protomer A
% CELL WORK 44 58 82 90 90 90
% GRID WORK 176 232 328
% XYZLIM -176 0 16 248 8 336
%
% protomer B
% CELL WORK 43 57 82 90 90 90
% GRID WORK 172 228 328
% XYZLIM -48 124 -20 208 -28 300