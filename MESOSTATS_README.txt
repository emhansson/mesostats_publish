This MESOSTATS_README.txt file was generated on 2022-06-27 by EM HANSSON

GENERAL INFORMATION

1. Title of Dataset: Mesostats

2. Contact Author Information
		Name: Erika M Hansson
		Institution: The University of Sheffield
		Email: emhansson1@sheffield.ac.uk

3. Date of data collection: 2018-2022 

4. Geographic location of data collection: Sheffield, South Yorkshire, United Kingdom 

5. Information about funding sources that supported the collection of the data: University of Sheffield PhD scholarship


SHARING/ACCESS INFORMATION

1. Licenses/restrictions placed on the data: CC-BY

2. Pre-print server link for paper: https://doi.org/10.1101/2022.04.20.488910


DATA & FILE OVERVIEW

1. File List:

mesostats.Rproj -- Rproject file

- Data
  \df.csv -- control data
      cell_count is cell density in cells/ml
      exp_day is days since experiment start
      chamber is biological replicate (population) within experiment
      test is experiment identifier
      rpm is pump speed, equivalent to dilution rate (2.5RPM = 0.3/day, 1.25 RPM = 0.15/day)
      chamberID is a unique identifier for chamber population among experiments
      
  \df2.csv -- clumping data
      cell_count is cell density in cells/ml
      exp_day is days since experiment start
      chamber is biological replicate (population) within experiment
      clumping identifies whether clumping was observed, N for no, Y for yes
      
  \df3.csv -- glyphosate pilot data
      cell_count is cell density in cells/ml
      exp_day is days since introduction of glyphosate treatment
      chamber is biological replicate (population) within experiment
      treatment_group is concentration of glyphosate

- Scripts
  \code.R -- R script document containing all analyses as well as plots included in manuscript

METHODOLOGICAL INFORMATION

1. Description of methods used for collection/generation of data: 
Flow cytometry (Beckman Coulter CytoFLEX)

2. Methods for processing the data: 
CytExpert (Beckman Coulter) was used to gate and count events detected in the PerCP-A channel (Excitation: 488nm, Emission: 690/50 BP). 