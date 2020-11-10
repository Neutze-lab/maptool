#!/bin/sh
# PROCESS MAPS OF BACTERIORHOPDOPSIN FOR ANALYSIS WITH MAPTOOL
#----------------------------------------------------------------------------------
# - Converts map to cartesian coordinates (with approx. 5 Å border around protein)
# - Converts map to .h5 format
# - Reads map dimensions using mapdump
#----------------------------------------------------------------------------------


# INPUT
# Note! Also set cell work, grid work, xyzlim in maprot command.
#----------------------------------------------------------------------------------
# Declare MAPTOOL's main folder
mainFolder=/home/adams/data_proc/mapTool/maptool-rhodopsin
# Declare map extension: '.map' or '.ccp4'
fileExtension='ccp4'

# Declare input, output, and scripts directories
inputDirectory=$mainFolder/input
outputDirectory=$mainFolder/output
scriptsDirectory=$mainFolder/scripts
#----------------------------------------------------------------------------------



# PROCESS FILES
#----------------------------------------------------------------------------------
# Generate 'output' folder plus 'log' and 'tmp' subfolders
mkdir -p $outputDirectory
mkdir -p $outputDirectory/log
mkdir -p $outputDirectory/tmp

# Navigate into input directory
cd $inputDirectory


# Find all input maps and process them
for map in *.$fileExtension 
do

# Read name of map files
mapName=${map%.*}
echo "Preparing file $map"


# 1. EXTEND INITIAL MAP TO COVER FULL CELL (necessary for mapprot to run)
#----------------------------------------
# Run MAPMASK from CCP4 suite
mapmask MAPIN $map \
        MAPOUT $outputDirectory/tmp/$mapName'__fullCell.ccp4' << EOF > $outputDirectory/log/$mapName'__fullCell.log'
XYZLIM CELL
EOF


# 2. TRANSLATE MAP TO CARTESIAN COORDINATES
#    Set CELL WORK / GRID WORK / XYZLIM to cover the volume of interest
#-----------------------------------------
# CELL WORK: A B C 90 90 90
#   ->  A B C sides of cell in Ångström (= size of pdb + 5Å border)
#   angles 90 degree 
# GRID WORK: 4*A 4*B 4*C 
#   -> 0.25 Å distance between grid points
# XYZLIM: 4*xmin 4*xmax 4*ymin 4*ymax 4*zmin 4*zmax
#-----------------------------------------
# The numbers below are obtained by running pdbSize.m on the resting state (.pdb)
maprot MAPIN $outputDirectory/tmp/$mapName'__fullCell.ccp4' \
       WRKOUT $outputDirectory/$mapName'_cartesian.map' << EOF >> $outputDirectory/log/$mapName'_cartesian.log'
MODE FROM
CELL WORK 43 57 82 90 90 90
GRID WORK 172 228 328
XYZLIM -48 124 -20 208 -28 300
SYMM WORK 1
AVER
ROTA POLAR 0 0 0
TRANS 0 0 0
EOF


# 3. CONVERT MAPS TO HDF5 (for matlab)
#----------------------------------------------------------------------------
$scriptsDirectory/map_to_h5.py $outputDirectory/$mapName'_cartesian.map'  $outputDirectory/$mapName'_cartesian.h5'


# 4. EXTRACT INFORMATION ON MAP DIMESIONS
#    cell dimensions, grid sampling, start and stop, axis order
#    and sigma (=RMSD) from new and original map
#    (Sigma of the original map is used in the matlab calculations since the 
#    sigma value for the cartesian map is slightly changed in the conversion.)
#----------------------------------------------------------------------------
# Example of data in mapdump output:
# Cell dimensions .................................   42.0000    54.0000    77.0000    90.0000    90.0000    90.0000
# Grid sampling on x, y, z ........................  168  216  308
# Start and stop points on columns, rows, sections  -136  172 -144   24 -176 
# Fast, medium, slow axes  3  1  2
# Rms deviation from mean density .................     0.01195
#----------------------------------------------------------------------------

mapdump MAPIN $outputDirectory/$mapName'_cartesian.map' << EOF > $outputDirectory/log/$mapName'_cartesian_header.txt'
EOF

mapdump MAPIN $map << EOF > $outputDirectory/log/$mapName'_original_header.txt'
EOF

grep 'Cell dimensions' $outputDirectory/log/$mapName'_cartesian_header.txt' | awk '{ print $4,$5,$6,$7,$8,$9}' > $outputDirectory/log/$mapName'_XYZinfo.dat'
grep 'Grid sampling on x, y, z' $outputDirectory/log/$mapName'_cartesian_header.txt' | awk '{ print $8,$9,$10}' >> $outputDirectory/log/$mapName'_XYZinfo.dat'
grep 'Start and stop points on columns, rows, sections' $outputDirectory/log/$mapName'_cartesian_header.txt' | awk '{ print $9,$10,$11,$12,$13,$14}' >> $outputDirectory/log/$mapName'_XYZinfo.dat'
grep '^ Fast, medium, slow axes' $outputDirectory/log/$mapName'_cartesian_header.txt' | awk '{ print $5,$6,$7}' >> $outputDirectory/log/$mapName'_XYZinfo.dat'

grep 'Rms deviation from mean density' $outputDirectory/log/$mapName'_cartesian_header.txt' | awk '{ print $7}'  >> $outputDirectory/log/$mapName'_XYZinfo.dat'
grep 'Rms deviation from mean density' $outputDirectory/log/$mapName'_original_header.txt' | awk '{ print $7}'  >> $outputDirectory/log/$mapName'_XYZinfo.dat'

done