% maptool_bR_ref
%% -------------------------------------------------------------------------
% written by Cecilia Wickstrand, 05-04-2020
% modified by Adams Vallejos, 09-11-2020
%
% publication: 
% "A tool for visualizing protein motions in time-resolved crystallography"
% published online 01-04-2020, in Structural Dynamics (Vol.7, Issue 2).
% https://doi.org/10.1063/1.5126921 
%
% For analysis of a set of difference Fourier electron density maps from 
% time-resolved serial femtosecond crystallography (TR-SFX) experiments on
% rhodopsin. 
% - Difference electron density within a sphere around every atom in 
%   the resting state structure is extracted 
% - The average positive and negative density over the sphere is
%   calculated for each map, with or without a sigma cut-off
% - The resulting dual amplitude function may be plotted along the trace of
%   selected atoms or used for further analysis
%
% This script is an example for analysis of three maps (16 ns, 760 ns and 1.725 ms).
% Before analysis, execute "process_maps.sh" to convert the maps to 
% cartesian coordinates and .h5 format.
% -------------------------------------------------------------------------


%% INPUT
% -------------------------------------------------------------------------
clear
clc
% SPHERE AND SIGMA SETTINGS
radius = 2; % Å
gridSpacing = 0.5; %Å, how dense grid within sphere
rmsCutoff = 3; % exclude data below this sigma level


% Files
% Main folder of MAPTOOL
mainPath = '/GIVE_PATH_TO_MAPTOOL/maptool/';% 

% Resting state pdb file
pdbFile = [mainPath '/input/foo.pdb'];

% Path to .h5 maps obtained from 'processMaps.sh'
inputDirectory = [mainPath 'output/'];

% Names of map files without extension
mapFileName = {'nameOfMapFile_1'; 'nameOfMapFile_2'};

% Corresponding time-points of maps above
timePointLabels = {'timePointOfMapFile_1', 'timePointOfMapFile_1'}';

%%  START CALCULATIONS
% ------------------------------------------------------------------------
% Clock starts ticking
tic;

% Number of atoms in pdb file
restingState = pdbread(pdbFile);

% -------------------------------------------------------------------------
% Only ATOM label in pdb file, here HETATM is substituted by ATOM
atomicCoordinates = [...
    [restingState.Model.Atom.X]'...
    [restingState.Model.Atom.Y]'...
	[restingState.Model.Atom.Z]'...
];

% With HETATM in pdb file
% atomicCoordinates = [...
%     [restingState.Model.Atom.X]'...
%     [restingState.Model.Atom.Y]'...
%     [restingState.Model.Atom.Z]'];...
%     [restingState.Model.HeterogenAtom.X]'...
%     [restingState.Model.HeterogenAtom.Y]'...
%     [restingState.Model.HeterogenAtom.Z]'];
% -------------------------------------------------------------------------

% Get number of atoms in the restingstate
numberOfAtoms = size(atomicCoordinates, 1);

% Calculate distance from sphere center
distanceWithinSphere = @(x,y,z) sqrt(x^2 + y^2 + z^2);

% -------------------------------------------------------------------------
% Generate spots of cubic grid with sphere inscribed
spots = linspace(-radius, radius, 2*(radius/gridSpacing) + 1);
numberOfSpots = length(spots);

% List grid spots within the sphere
sphereList= zeros(numberOfSpots^3, 7);
count = 0;
for i = 1 : numberOfSpots
   for j = 1 : numberOfSpots
        for k = 1 : numberOfSpots
            distance = distanceWithinSphere( spots(i) , spots(j) , spots(k) );
            if distance <= radius
               count = count + 1;
               sphereList(count,:) = [spots(i) spots(j) spots(k) i j k distance];
            end
        end
   end
end
sphereList = sphereList(1:count, :); 
numberOfPoints = size(sphereList, 1);
% -------------------------------------------------------------------------

% PRECALCULATE ALL COORDINATES IN ALL SPHERES
% (example: 1834 rows = atoms, 2109 columns = points) 
X = repmat(atomicCoordinates(:,1), 1, numberOfPoints)+repmat(sphereList(:,1)', numberOfAtoms, 1);
Y = repmat(atomicCoordinates(:,2), 1, numberOfPoints)+repmat(sphereList(:,2)', numberOfAtoms, 1);
Z = repmat(atomicCoordinates(:,3), 1, numberOfPoints)+repmat(sphereList(:,3)', numberOfAtoms, 1);

% Reshape to a single column starting with first point for each atom, then second etc. 
X=reshape(X, numberOfAtoms*numberOfPoints, 1);
Y=reshape(Y, numberOfAtoms*numberOfPoints, 1);
Z=reshape(Z, numberOfAtoms*numberOfPoints, 1);

% LOAD GRID FROM MAP
% (using first experimental XYZinfo file - here all files have the same grid)
XYZinfo = dlmread([inputDirectory 'log/' mapFileName{1} '_XYZinfo.dat']);

% Dimensions of cartesian cell
cellDimensions = XYZinfo(1,1:3);
% Grid dimensions of cartesian map
gridPoints = XYZinfo(2,1:3); 
% Boundaries of cell axes
axesBoundaries = XYZinfo(3,:);
% Order of axes
axesOrder = XYZinfo(4,1:3);

% Grid spacing
dX = cellDimensions(1)/gridPoints(1);
dY = cellDimensions(2)/gridPoints(2);
dZ = cellDimensions(3)/gridPoints(3);

% Re-arrange axes according to ordering of coordinates order from maps
axes = sortrows([axesOrder(1) axesBoundaries(1:2);
                 axesOrder(2) axesBoundaries(3:4);
                 axesOrder(3) axesBoundaries(5:6)]);

% Number of points on each side of the grid
sX = (axes(1,2) : axes(1,3)) * dX;
sY = (axes(2,2) : axes(2,3)) * dY;
sZ = (axes(3,2) : axes(3,3)) * dZ;

% 3D grid coordinates with dimensions from above
[gY, gZ, gX] = meshgrid(sY, sZ, sX);


% LOAD MAPS AND EXTRACT DENSITIES
% ------------------------------------------------------------------------
numberOfMaps = length(mapFileName);

% Read rms(sigma) from maps
rmsFromMaps = zeros(numberOfMaps, 2);

% Populate with mal densities
normalizedDensities = zeros(numberOfMaps, numberOfAtoms, numberOfPoints);
for m = 1:numberOfMaps
    formatSpec= 'Currently at map %s time %g s\n';
    fprintf(formatSpec,mapFileName{m},toc(tic))

%   Load sigma info
    XYZinfo = dlmread([inputDirectory 'log/' mapFileName{m} '_XYZinfo.dat']);
    rmsFromMaps(m,1:2) = XYZinfo(5:6,1);

%   LOAD MAP AND CALCULATE DENSITY AT POINTS BY INTERPOLATION
%   reshape back to: rows = atoms, columns = points
    cartesianMap = h5read([inputDirectory mapFileName{m} '_cartesian.h5'],'/map');
    interpolatedDensities = interp3(gY,gZ,gX,cartesianMap,Y,Z,X);

    mapDensities = reshape(interpolatedDensities, numberOfAtoms, numberOfPoints);
%   Divide with sigma of original map
    normalizedDensities(m,:,:) = mapDensities/rmsFromMaps(m,2); 
end


% CALCULATE AVERAGE DENSITIES AND CORRELATIONS
%------------------------------------
% Calculate mean positive and negative densities
mapDensitiesAbsolute = normalizedDensities;
mapDensitiesAbsolute(abs(normalizedDensities) < rmsCutoff) = 0;

mapDensities_positive = mapDensitiesAbsolute;
mapDensities_positive(mapDensitiesAbsolute < 0) = 0;
meanpositiveDensities = mean(mapDensities_positive,3);

mapDensities_negative = mapDensitiesAbsolute;
mapDensities_negative(mapDensitiesAbsolute > 0) = 0;
meanNegativeDensitues = mean(mapDensities_negative,3);

% Calculate pearson correlation (<A+> <A->, <B+> <B->)
pearsonCorrelation = zeros(numberOfMaps,numberOfMaps);

for m = 1:numberOfMaps
    for n = 1:numberOfMaps
        pearsonCorrelation(m,n) = corr2(...
            [meanpositiveDensities(m,:) meanNegativeDensitues(m,:)],...
            [meanpositiveDensities(n,:) meanNegativeDensitues(n,:)]);
    end
end

% PREPARE PLOTS
%------------------------------------
fprintf('Preparing plots.\n')

% Plot settings
% RGB Colors normalized to 256
golden = [0.7, 0.5, 0];
purple = [0.5, 0.4, 1];
silver = [0.8, 0.8, 0.8];
fontSize = 12;
fontName = 'helvetica narrow';
lineWidth = 1; % Default value 0.5
labelAll= {'1', '2', '3', '4', '5', '6','7','8'};
labelCA = {'1', '2', '3', '4', '5', '6','7','8'};
% Set fonts in all figures
set(0,'DefaultAxesFontName',fontName,'DefaultTextFontName',fontName);

% Avoid overlapping in plots        
scale_E = max(max(meanpositiveDensities-meanNegativeDensitues));

% For plotting all atoms, find where helices start and stop
bR_helix = [35 64;
            71 100;
            108 138;
            150 174;
            201 228;
            246 276;
            286 308;
            310 321];

% Get residues index from resting state
allResidues = [restingState.Model.Atom.resSeq]'; 

% Select residues from limits above
limitsHelix = zeros(size(bR_helix));
for h = 1:size(bR_helix,1)
    limitsHelix(h,1) = find(allResidues==bR_helix(h,1), 1, 'first');
    limitsHelix(h,2) = find(allResidues==bR_helix(h,2), 1, 'last');
end

% For plotting C alphas, find where helices start and stop
selectedCA = find(strcmp({restingState.Model.Atom.AtomName}','CA'));
residuesCA = [restingState.Model.Atom(selectedCA).resSeq]'; 

bR_helix = [35 64;
            71 100;
            108 138;
            150 174;
            201 228;
            246 276;
            286 308;
            310 321];

limitsCA = zeros(size(bR_helix));
for h = 1:size(bR_helix,1)
    limitsCA(h,1) = find(residuesCA==bR_helix(h,1), 1, 'first');
    limitsCA(h,2) = find(residuesCA==bR_helix(h,2), 1, 'last');
end

% Select functionally relevant residues at the active site
residuesOfInterest = [82 85 89 182 212 216 300];

selectedResidues = find(ismember(allResidues, residuesOfInterest));
limitsResidues = zeros(length(residuesOfInterest),2);
ticks = cell(1,length(residuesOfInterest));
for i = 1:length(limitsResidues)
    limitsResidues(i,1) = find(allResidues(selectedResidues)==residuesOfInterest(i), 1, 'first');
    limitsResidues(i,2) = find(allResidues(selectedResidues)==residuesOfInterest(i), 1, 'last');
    ticks{i} = num2str(residuesOfInterest(i));
end
ticks{end} = 'RET';


% PLOT
%------------------------------------
% figure('units','normalized','outerposition',[0 0 1 1],'name',['radius '...
%num2str(radius) ' Å, ' num2str(rmsCutoff) ' sigma,'])
figure('name',['radius ' num2str(radius) ' Å, ' num2str(rmsCutoff) ' sigma,'])


% 1. mean pos/neg density, selected residues around site
subplot(2,3,1)
set_scale_E = 4;
    hold all
    for n = 1:numberOfMaps
        line([1 length(selectedResidues)],...
             [1/set_scale_E+(n) 1/set_scale_E+(n)],...
             'color', silver,...
             'linestyle','--')

        line([1 length(selectedResidues)],...
             [-1/set_scale_E+(n) -1/set_scale_E+(n)],...
             'color', silver,...
             'linestyle','--')

        plot(1:length(selectedResidues),...
            -meanpositiveDensities(n,selectedResidues)/set_scale_E+(n),...
            'color', purple,...
            'LineWidth',lineWidth) 
        
        plot(1:length(selectedResidues),...
             -meanNegativeDensitues(n,selectedResidues)/set_scale_E+(n),...
             'color', golden,...
             'LineWidth',lineWidth)
    end

    for i = 1:numberOfMaps
        line([1 length(selectedResidues)], [(i) (i)],'color', silver)
    end

    for h = 2:length(limitsResidues)
        line([limitsResidues(h,1)-0.5 limitsResidues(h,1)-0.5],...
             [0 numberOfMaps+1],'color', 'k')
    end 
    ylim([0 numberOfMaps+1])
    xlim([1 length(selectedResidues)])
    title('Selected residues')   
    set(gca,...
        'XTickLabel',ticks,...
        'XTick', mean(limitsResidues,2)',...
        'Ytick',1:numberOfMaps,...
        'Yticklabel', timePointLabels)

    set(gca,...
        'Ydir','reverse',...
        'XAxisLocation', 'top',...
        'TickDir','out',...
        'box','on',...
        'FontSize',fontSize)

 % 2. Mean pos/neg density, all atoms    
 subplot(2,3,[2 3])
    hold all
    for n = 1:numberOfMaps
        plot(1:numberOfAtoms,...
             -meanpositiveDensities(n,:)/scale_E+(n),...
             'color', purple,...
             'LineWidth',lineWidth)
         
        plot(1:numberOfAtoms,...
             -meanNegativeDensitues(n,:)/scale_E+(n),...
             'color', golden,...
             'LineWidth',lineWidth)
    end

    for i = 1:numberOfMaps
        line([1 numberOfAtoms], [(i) (i)],'color', silver)
        hold all
    end 
    
    for h = 1:length(limitsHelix)
        line([limitsHelix(h,1)-0.5 limitsHelix(h,1)-0.5],[0 numberOfMaps+1],'color', 'k')
        line([limitsHelix(h,2)+0.5 limitsHelix(h,2)+0.5],[0 numberOfMaps+1],'color', 'k')
    end 
    set(gca,'TickDir','out','Ytick',1:numberOfMaps,'Yticklabel', timePointLabels)
    set(gca,'XTickLabel',labelAll,'XTick', mean(limitsHelix,2))
    ylim([0 numberOfMaps+1])
    xlim([1 numberOfAtoms])
    title('All atoms (protein only)')
    set(gca,'Ydir','reverse', 'XAxisLocation', 'top', 'box','on','FontSize',fontSize)

% 3. Mean pos/neg density, C-alpha atoms
subplot(2,3,[5 6])
    hold all
    for n = 1:numberOfMaps
        plot(1:length(selectedCA),...
             -meanpositiveDensities(n,selectedCA)/scale_E+(n),...
             'color', purple,...
             'LineWidth',lineWidth)
         
        plot(1:length(selectedCA),...
             -meanNegativeDensitues(n,selectedCA)/scale_E+(n),...
             'color', golden,...
             'LineWidth',lineWidth)
    end

    for i = 1:numberOfMaps
        line([1 length(selectedCA)], [(i) (i)],'color', silver)
    end

    for h = 1:length(limitsCA)
        line([limitsCA(h,1)-0.5 limitsCA(h,1)-0.5],[0 numberOfMaps+1],'color', 'k')
        line([limitsCA(h,2)+0.5 limitsCA(h,2)+0.5],[0 numberOfMaps+1],'color', 'k')
    end 
    set(gca,'XTickLabel',labelCA,'XTick', mean(limitsCA,2)')
    set(gca,'TickDir','out','Ytick',1:numberOfMaps,'Yticklabel', timePointLabels)
    ylim([0 numberOfMaps+1])
    xlim([1 length(selectedCA)])
    title('C alphas')
    set(gca,'Ydir','reverse', 'XAxisLocation', 'top', 'box','on','FontSize',fontSize)

% 4. Pearson correlation
subplot(2,3,4)
    imagesc(pearsonCorrelation)
    colorbar
    axis('square')
    caxis([0 1])
    set(gca,...
        'TickDir','out',...
        'Ydir','reverse',...
        'XAxisLocation', 'top',...
        'box','on',...
        'FontSize', fontSize)

    set(gca,...
        'Ytick',1:numberOfMaps,...
        'Yticklabel',timePointLabels,...
        'Xtick',1:numberOfMaps,...
        'Xticklabel', timePointLabels)
    title('Correlation')
    
% OPTIONAL CONTROL SECTION
%**************************************************************************
% - CHECK GRID
%   load a map
%       map3 = h5read([indir mapnames{3} '.h5'],'/map');
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
%     map3 = h5read([indir mapnames{3} '_cartesian.h5'],'/map');
%     testdensities = interp3(gY,gZ,gX,map3, aY, aZ , aX)
% 
%     For this test set (THR 24 C; TYR 43 OH; GLY 63 C; VAL 217 O)
%     the calculated densites (-0.0269    0.0157    0.0657   -0.0954)
%     are identical to those seen in coot.
%**************************************************************************

   
