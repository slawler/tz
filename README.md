# tz

Development repo of scripts to analyze & predict surface water in transition zones. 

### Directories
1. harmonics: 
	data &  binaries for extracting harmonic constituents from ADCIRC runs
	[source code](http://adcirc.org/products/adcirc-tidal-databases/)

2. predictions:
    data & binaries for producing tidal predictions using NOAA's ntp4 program
	[source code](https://tidesandcurrents.noaa.gov/faq2.html#65)

3. utilities: python scripts to autmoate processing 

4. validation: python scripts to retrieve NOAA coops predictions for validation.

### Usage Notes
##### 1. Extract Harmonic Constituents for points of interest (poi):
a. cd to harmonics, see usage notes

##### 2. Create Control File(s) 
a. cd to utilities, run MakeCTL.py (follow instructions in script header)

##### 3. Run NTP4 to create tide tables
a. cd to predictions, see usage notes

##### 4. Format NTP4 output, plot results
a. cd to utilities, run PlotTides.py

##### 5. Validate results against gaged locations
a. Under Constructions==> add script to write new file for:
	1. Official Predictions
	2. Predcitions created using ADCIRC results & programs included in this repo	
