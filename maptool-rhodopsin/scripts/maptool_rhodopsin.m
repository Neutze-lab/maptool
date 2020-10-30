% maptool_rhodopsin
% -------------------------------------------------------------------------
% written by Cecilia Wickstrand
% Modified by Adams Vallejos
%
% publication: 
% "A tool for visualizing protein motions in time-resolved crystallography"
% published online 01-04-2020, in Structural Dynamics (Vol.7, Issue 2).
% https://doi.org/10.1063/1.5126921 
%
% For analysis of a set of difference Fourier electron density maps from 
% time-resolved serial femtosecond crystallography (TR-SFX) experiments on
% bacteriorhodopsin. 
% - Difference electron density within a sphere around every atom in 
%   the resting state structure is extracted 
% - The average positive and negative density over the sphere is
%   calculated for each map, with or without a sigma cut-off
% - The resulting dual amplitude function may be plotted along the trace of
%   selected atoms or used for further analysis
%
% -------------------------------------------------------------------------


% INPUT
% -------------------------------------------------------------------------

% Sphere and sigma settings
radius = 2; %Å
distance = 0.5; %Å, how dense grid within sphere
sigmacutoff = 3; % exclude data below this sigma level


% Files
here = '/PATH/TO//maptool-rhodopsin/';
pdbpath = [here '/input/PDBFILE.pdb']; % resting state pdb
indir = [here 'output/']; % where to find .h5 maps

mapnames = {'FOO_1'; 'FOO_2'};

timeticks = {'TIMEPOINT_1', 'TIMEPOINT_2'}';



% START CALCULATIONS
% ------------------------------------------------------------------------
time = tic;


% LOAD ATOMS
pdb = pdbread(pdbpath);

% In this case all HETATM are written as ATOM in pdb
atomcoord = [[pdb.Model.Atom.X]' [pdb.Model.Atom.Y]' [pdb.Model.Atom.Z]'];...
             %[pdb.Model.HeterogenAtom.X]' [pdb.Model.HeterogenAtom.Y]' [pdb.Model.HeterogenAtom.Z]'];
nratoms =  size(atomcoord,1);


% CALCULATE SPHERE COORDINATES
calcdist = @(x,y,z) sqrt(x^2 + y^2 + z^2);

% radius 2 distance 0.5 -> spots = [-2 -1.5 -1 -0.5 0 0.5 1 1.5 2]
spots = distance:distance:ceil(radius/distance)*distance;
spots = [-fliplr(spots) 0 spots];
nrspots = length(spots);

% spherelist column 1 2 3 = coordinates
spherelist= zeros(nrspots^3,7);
count = 0;
for i = 1:nrspots
   for j = 1:nrspots
        for k = 1:nrspots
            dist = calcdist(spots(i),spots(j),spots(k));
            if dist <= radius
               count = count + 1;
               spherelist(count,:) = [spots(i) spots(j) spots(k) i j k dist];
            end
        end
   end
end
spherelist = spherelist(1:count,:); 
nrpoints = size(spherelist,1);


% PRECALCULATE ALL COORDINATES IN ALL SPHERES
% (example: 1834 rows = atoms, 2109 columns = points) 
X = repmat(atomcoord(:,1),1,nrpoints)+repmat(spherelist(:,1)',nratoms,1);
Y = repmat(atomcoord(:,2),1,nrpoints)+repmat(spherelist(:,2)',nratoms,1);
Z = repmat(atomcoord(:,3),1,nrpoints)+repmat(spherelist(:,3)',nratoms,1);

% Reshape to a single column starting with first point for each atom, then second etc. 
X=reshape(X, nratoms*nrpoints,1);
Y=reshape(Y, nratoms*nrpoints,1);
Z=reshape(Z, nratoms*nrpoints,1);


% LOAD GRID
% (using first experimental XYZinfo file - here all files have the same grid)
XYZinfo = dlmread([indir 'log/' mapnames{1} '_XYZinfo.dat']);

celldimensions = XYZinfo(1,1:3);
gridpoints = XYZinfo(2,1:3); 
axislimits = XYZinfo(3,:);
axisorder = XYZinfo(4,1:3);

% grid distance
dX = celldimensions(1)/gridpoints(1);
dY = celldimensions(2)/gridpoints(2);
dZ = celldimensions(3)/gridpoints(3);

tmp = [axisorder(1) axislimits(1:2); axisorder(2) axislimits(3:4); axisorder(3) axislimits(5:6)];
axes = sortrows(tmp);

% points at each side
sX = [axes(1,2):axes(1,3)]*dX;
sY = [axes(2,2):axes(2,3)]*dY;
sZ = [axes(3,2):axes(3,3)]*dZ;

[gY,gZ,gX] = meshgrid(sY,sZ,sX);


% LOAD MAPS AND EXTRACT DENSITIES
% ------------------------------------------------------------------------
nrmaps = length(mapnames);
sigma = zeros(nrmaps, 2);
mapd0 = zeros(nrmaps, nratoms, nrpoints);

for m = 1:nrmaps
    fprintf(['Currently at map ' num2str(m) ' time ' num2str(toc(time)) ' s\n'])

    % Load sigma info
    XYZinfo = dlmread([indir 'log/' mapnames{m} '_XYZinfo.dat']);
    sigma(m,1:2) = XYZinfo(5:6,1);

    
    % LOAD MAP AND CALCULATE DENSITY AT POINTS BY INTERPOLATION
    % reshape back to rows = atoms, columns = points
    map = hdf5read([indir mapnames{m} '_cartesian.h5'],'map');
    tmp = interp3(gY,gZ,gX,map,Y,Z,X);

    mapdensities = reshape(tmp, nratoms, nrpoints);
    % divide with sigma of original map
    mapd0(m,:,:) = mapdensities/sigma(m,2); 
end


% CALCULATE AVERAGE DENSITIES AND CORRELATIONS
%------------------------------------
% Calculate mean positive and negative densities
mapd = mapd0;
mapd(abs(mapd0) < sigmacutoff) = 0;

mapd_pos = mapd;
mapd_pos(mapd < 0) = 0;
meanposden = mean(mapd_pos,3);

mapd_neg = mapd;
mapd_neg(mapd > 0) = 0;
meannegden = mean(mapd_neg,3);

% Calculate Pscore = pearson correlation (<A+> <A->, <B+> <B->)
Pscore = zeros(nrmaps,nrmaps);

for m = 1:nrmaps
    for n = 1:nrmaps
        Pscore(m,n) = corr2([meanposden(m,:) meannegden(m,:)],[meanposden(n,:) meannegden(n,:)]);
    end
end


% PREPARE EXAMPLE PLOTS
%------------------------------------
fprintf('Preparing plots.\n')

% Plot settings
golden = [0.83, 0.65, 0.13]*0.8;
slate = [0.5 0.5 1]*0.8;
fontsize = 14;
fontname = 'helvetica narrow';
set(0,'DefaultAxesFontName',fontname,'DefaultTextFontName',fontname);

% avoid overlapping in plots        
scale_E = max(max(meanposden-meannegden));

% For plotting all atoms, find where helices start and stop
helix_limits_res = [35 64; 71 100; 108 138; 150 174; 201 228; 246 276; 286 308; 310 321];
res_all = [pdb.Model.Atom.resSeq]'; 
limits_all = zeros(size(helix_limits_res));
for h = 1:size(helix_limits_res,1)
    limits_all(h,1) = find(res_all==helix_limits_res(h,1), 1, 'first');
    limits_all(h,2) = find(res_all==helix_limits_res(h,2), 1, 'last');
end

% For plotting C alphas, find where helices start and stop
selected_Ca = find(strcmp({pdb.Model.Atom.AtomName}','CA'));
res_Ca = [pdb.Model.Atom(selected_Ca).resSeq]'; 
helix_limits_res = [35 64; 71 100; 108 138; 150 174; 201 228; 246 276; 286 308; 310 321];
limits_Ca = zeros(size(helix_limits_res));
for h = 1:size(helix_limits_res,1)
    limits_Ca(h,1) = find(res_Ca==helix_limits_res(h,1), 1, 'first');
    limits_Ca(h,2) = find(res_Ca==helix_limits_res(h,2), 1, 'last');
end

% For plotting selected residues, find where residues start and stop
res = [pdb.Model.Atom.resSeq]';
resofinterest = [82 85 89 182 212 216 300];
selected_res = find(ismember(res, resofinterest));
limits_res = zeros(length(resofinterest),2);
for i = 1:length(limits_res)
    limits_res(i,1) = find(res(selected_res)==resofinterest(i), 1, 'first');
    limits_res(i,2) = find(res(selected_res)==resofinterest(i), 1, 'last');
    ticks_res{i} = num2str(resofinterest(i));
end
ticks_res{end} = 'RET';


% PLOT
%------------------------------------
figure('units','normalized','outerposition',[0 0 1 1],'name',['radius ' num2str(radius) ' Å, ' num2str(sigmacutoff) ' sigma,'])

% 1. mean pos/neg density, selected residues around site
subplot(2,3,1)
scale_E_set = 4;
    hold all
    for n = 1:nrmaps
        line([1 length(selected_res)], [1/scale_E_set+(n) 1/scale_E_set+(n)],'color', [0.8 0.8 0.8],'linestyle','--')
        line([1 length(selected_res)], [-1/scale_E_set+(n) -1/scale_E_set+(n)],'color', [0.8 0.8 0.8],'linestyle','--')
        plot(1:length(selected_res), -meanposden(n,selected_res)/scale_E_set+(n),'color', slate) 
        plot(1:length(selected_res), -meannegden(n,selected_res)/scale_E_set+(n),'color', golden)
    end
    for i = 1:nrmaps
        line([1 length(selected_res)], [(i) (i)],'color', [0.8 0.8 0.8])
    end
    for h = 2:length(limits_res)
        line([limits_res(h,1)-0.5 limits_res(h,1)-0.5],[0 nrmaps+1],'color', 'k')
    end 
    ylim([0 nrmaps+1])
    xlim([1 length(selected_res)])
    title('Selected residues')   
    set(gca,'XTickLabel',ticks_res,'XTick', mean(limits_res,2)','Ytick',1:nrmaps,'Yticklabel', timeticks)
    set(gca,'Ydir','reverse', 'XAxisLocation', 'top','TickDir','out', 'box','on','FontSize',fontsize)

 % 2. mean pos/neg density, all atoms    
 subplot(2,3,[2 3])
    hold all
    for n = 1:nrmaps
        plot(1:nratoms, -meanposden(n,:)/scale_E+(n),'color', slate,'LineWidth',1.5) 
        plot(1:nratoms, -meannegden(n,:)/scale_E+(n),'color', golden,'LineWidth',1.5)
    end
    for i = 1:nrmaps
        line([1 nratoms], [(i) (i)],'color', [0.8 0.8 0.8])
        hold all
    end 
    for h = 1:length(limits_all)
        line([limits_all(h,1)-0.5 limits_all(h,1)-0.5],[0 nrmaps+1],'color', 'k')
        line([limits_all(h,2)+0.5 limits_all(h,2)+0.5],[0 nrmaps+1],'color', 'k')
    end 
    set(gca,'TickDir','out','Ytick',1:nrmaps,'Yticklabel', timeticks)
    set(gca,'XTickLabel',{'1', '2', '3', '4', '5', '6','7','8'},'XTick', mean(limits_all,2))
    ylim([0 nrmaps+1])
    xlim([1 2527])% A 2565 %  B 2527
    title('All atoms (protein only)')
    set(gca,'Ydir','reverse', 'XAxisLocation', 'top', 'box','on','FontSize',fontsize)
    
% 3. mean pos/neg density, C alpha atoms
subplot(2,3,[5 6])
    hold all
    for n = 1:nrmaps
        plot(1:length(selected_Ca), -meanposden(n,selected_Ca)/scale_E+(n),'color', slate,'LineWidth',1.5) 
        plot(1:length(selected_Ca), -meannegden(n,selected_Ca)/scale_E+(n),'color', golden,'LineWidth',1.5)
    end
    for i = 1:nrmaps
        line([1 length(selected_Ca)], [(i) (i)],'color', [0.8 0.8 0.8])
    end   
    for h = 1:length(limits_Ca)
        line([limits_Ca(h,1)-0.5 limits_Ca(h,1)-0.5],[0 nrmaps+1],'color', 'k')
        line([limits_Ca(h,2)+0.5 limits_Ca(h,2)+0.5],[0 nrmaps+1],'color', 'k')
    end 
    set(gca,'XTickLabel',{'1', '2', '3', '4', '5', '6','7','8'},'XTick', mean(limits_Ca,2)')
    set(gca,'TickDir','out','Ytick',1:nrmaps,'Yticklabel', timeticks)
    ylim([0 nrmaps+1])
    xlim([1 length(selected_Ca)])
    title('C alphas')
    set(gca,'Ydir','reverse', 'XAxisLocation', 'top', 'box','on','FontSize',fontsize)

% 4. Pearson correlation
subplot(2,3,4)
    imagesc(Pscore)
    colorbar
    axis('square')
    caxis([0 1])
    set(gca,'TickDir','out','Ydir','reverse', 'XAxisLocation', 'top', 'box','on','FontSize', fontsize)
    set(gca,'Ytick',1:nrmaps,'Yticklabel', timeticks,'Xtick',1:nrmaps,'Xticklabel', timeticks)
    title('Correlation')
    

    
% OPTIONAL CONTROL SECTION
%**************************************************************************
% - CHECK GRID
%   load a map
%       map3 = hdf5read([indir mapnames{3} '.h5'],'map');
%   check that the size of the map is the same as for gX / gY / gZ 
%   if not, change the order of Y Z X in the meshgrid command 
%
% - CHECK PDB IS COVERED BY MAP
%   pdb = pdbread(pdbpath)
%   lim_map = [sX(1) sX(end) sY(1) sY(end) sZ(1) sZ(end)]
%   lim_pdb = [min([pdb.Model.Atom(:).X]) max([pdb.Model.Atom(:).X]) min([pdb.Model.Atom(:).Y]) max([pdb.Model.Atom(:).Y]) min([pdb.Model.Atom(:).Z]) max([pdb.Model.Atom(:).Z])]
%
% - CHECK DENSITY CALCULATIONS
%   Select a map and a pick few random atoms
%   Calculate the interpolated density value at each of the atoms
%   Open coot, go to the atoms and set the contour level to see that the
%   calculated densities are correct.
% 
%     testset = [157;307;457;1657]
%     testpdb.Model.Atom = pdb.Model.Atom(testset)
%     aX = [testpdb.Model.Atom.X]';
%     aY = [testpdb.Model.Atom.Y]';
%     aZ = [testpdb.Model.Atom.Z]';
%     map3 = hdf5read([indir mapnames{3} '_cartesian.h5'],'map');
%     testdensities = interp3(gY,gZ,gX,map3, aY, aZ , aX)
% 
%     For this test set (THR 24 C; TYR 43 OH; GLY 63 C; VAL 217 O)
%     the calculated densites (-0.0269    0.0157    0.0657   -0.0954)
%     are identical to those seen in coot.
%**************************************************************************