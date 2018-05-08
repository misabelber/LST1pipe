This repository stores codes for the implementation of the analysis of LST1 data in ctapipe.

====================DEPENDENCIES===================================================================
This repository works in the framework of ctapipe, activate environment:
     - source activate cta-dev

====================CONFIGURATION==================================================================

Change "DATA_PATH" in the script "stractdata.sh"

====================INSTRUCTIONS==================================================================

"LST_Hillas.py"
	Usage: python LST_Hillas [Particle] [Filename]
	
	Particle: Gamma, Proton, Electron.
	
	Filename: Simtelarray file containing raw event data.

	This program reduces the data until obtaining the image Hillas Parameters and then stores 
	them in a fits and ascii file, losing pixel information. 

"extractdata.sh"
	Allows running "LST_Hillas.py" over all files contained in a certain folder. 
	
