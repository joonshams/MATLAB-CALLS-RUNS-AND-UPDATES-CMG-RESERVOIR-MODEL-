This folder contains early EnKF data assimilation trials that I rand using a lab-scale
reservoir model 'Sandbox'. More specfically, this section of the code calls,runs,and updates 
CMG files (.out,.dat) to generate forecast models. 

main script (EnKFfprvpqandpa.m) 

NOTE: The main script uses RPHTOOLS that are available open-source:
 (https://pangea.stanford.edu/departments/geophysics/dropbox/SRB/public/data/RPHtools.htm)
 
 You will also need a CMG licence (i use 2019/2017) and you will have to make those changes
 in SWINT() and functions that SWINT() calls. 
 CMG FILES NEEDED: gm_version.exe,SR3SimInterface.dll,libiomp5md.dll,and binarrayfile.dll 

The main file uses a series of observations: seismic attributes
(p-wave seismic velocity and attenuation) and pressure (kPa) to assimilate and update
pressure and saturation fields.

CMG FILES: sandtank.(dat/irf/mrf/out/rst/sr3)

Only .dat,.out are readable in notepad and are used to change the gas saturation and
pressure properties to generate forecasts.

------------------------------------------------------------
Shams Joon
Ph.D. Candidate| Pronouns: he/him/his
Department of Energy and Mineral Engineering
svj5235@psu.edu
