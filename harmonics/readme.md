# Harmonic Constituents
Included: Code and related files for extracting harmonic constituents from adcirc tidal runs.  

#### Files
1. Source code: ADCIRC_db_extract_2012.F90
2. Executables:
	1. dbex (linux)
	2. dbex_w7 (windows 7)

3. .mod files: required for harmonic extraction
4. Required files: poi.in, fort.14, fort.53

#### Usage Notes
	1. Add required files to working directory
	2. Select option 1 & follow the prompts.
	3. In Windows OS, executable created using cygwin; file endings can be buggy. In Linux OS, recommend to re-complie (.mod files generated on compilations).
	4. Move output ('elev_hc.out') to outputs directory 
