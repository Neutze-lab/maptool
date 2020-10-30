clear
close all

% Find limits in pdb x y z 
%pdb = pdbread('/home/cecilia/Desktop/Documents/StrRef/BR/SACLA2016/map_0510/bR1_refine_106_without_H.pdb');
% X = [pdb.Model.Atom(:).X]'; [pdb.Model.HeterogenAtom(:).X]';
% Y = [pdb.Model.Atom(:).Y]'; [pdb.Model.HeterogenAtom(:).Y]';
% Z = [pdb.Model.Atom(:).Z]'; [pdb.Model.HeterogenAtom(:).Z]';

pdb = pdbread('/home/adams/data_proc/mapTool/maptool-rhodopsin/input/RHO_chainB.pdb')
border = 5; 
dist = 0.25;
%%
X = [pdb.Model.Atom(:).X]'; 
Y = [pdb.Model.Atom(:).Y]'; 
Z = [pdb.Model.Atom(:).Z]'; 
%%
pdblim = [min(X) max(X) min(Y) max(Y) min(Z) max(Z)]
%%
% limits with border
limits = ceil([pdblim(1)-border pdblim(2)+border pdblim(3)-border pdblim(4)+border pdblim(5)-border pdblim(6)+border])
%%
cellwork = [limits(2)-limits(1) limits(4)-limits(3) limits(6)-limits(5) 90 90 90]
%%
gridwork = cellwork(1:3)*1/dist
%%
xyzlim = limits/0.25
%%
% rhodopsin A
% CELL WORK 44 58 82 90 90 90
% GRID WORK 176 232 328
% XYZLIM -176 0 16 248 8 336
%
% rhodopsin B
% CELL WORK 43 57 82 90 90 90
% GRID WORK 172 228 328
% XYZLIM -48 124 -20 208 -28 300
%

