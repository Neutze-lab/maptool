# maptool-bR

Maptool was developed for analysis of a set of difference Fourier electron density maps from time-resolved serial femtosecond crystallography (TR-SFX) experiments on bacteriorhodopsin. 
The 3D maps are represented in 1D, to aid interpretation on where and when structural changes occur, by the following:
  - Difference electron density within a sphere around every atom in the resting state structure is extracted 
  - The average positive and negative densities over the sphere are calculated for each map, with or without a sigma cut-off
The resulting dual amplitude function may be plotted along the trace of selected atoms or used for further analysis.

### If you use the scripts, please cite:
> [Wickstrand et al. (2020), "A tool for visualizing protein motions in time-resolved crystallography". Published online 01-04-2020 in Structural Dynamics (Vol.7, Issue 2).](https://doi.org/10.1063/1.5126921)
---

The uploaded set of scripts give an example for analysis of three maps (16 ns, 760 ns and 1.725 ms). 

Scripts:
- pdbSize.m
	Obtains grid dimensions for MAPROT
- processMaps.sh	
	Translates maps to cartesian coordinates
	Converts cartesian map files to h5 format
	Saves information on map dimensions etc.

- map_to_h5.py (used in processMaps.sh)
	Used for conversion from map to h5 format

- maptool_bR.m
	The main script that performs all calculations and plots the dual amplitude functions 

To run the analysis:
- Download the files and change the path to the "maptool" directory within the scripts
- Execute processMaps.sh 
- Open matlab and run maptool_bR.m

Comments on using the tool on other datasets:
- processMaps.sh
	If all angels are already 90 degrees, remove the conversion to cartesian coordinates.
	If a map is not covering the atoms of interest including a border to cover the spheres it may be extended using mapmask BORDER.

- maptool_bR.m
 	Maps from different experiments may have different grid dimensions and different resting state pdbs. The script may be rewritten to load a new grid and/or pdb for each map.
