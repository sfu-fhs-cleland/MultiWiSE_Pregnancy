# MultiWiSE_Pregnancy
Code for calculating the MultiWiSE metrics during pregnancy

## Code
This code provides an example for how to use the MultiWiSE approach to estimate WFS exposure during pregnancy. The script follows the approach outlined in the manuscript 'Multiyear Wildfire Smoke Exposure (MultiWiSE) metrics:  A data-driven approach to characterizing episodic PM2.5 exposures for epidemiologic research on chronic health effects', and is adapted from the code provided in the 'MultiWiSE.R' script (https://github.com/sfu-fhs-cleland/MultiWiSE_Metrics_Manuscript). The script was modified to use information on key maternal dates (date of conception and date of birth) and individual-level residential history to assemble a individual-level PM2.5 exposure profile and calculate the MultiWiSE metrics during pregnancy and each trimester. This approach uses data on daily total PM2.5 in the 3 years prior to the date of birth to calculate the counterfactual value and modified z-scores, which are then used to estimate WFS and non-WFS PM2.5 and calculate the MultiWiSE metrics during pregnancy. 

## Input Data
Example input data is provided to run the script. The example maternal dates are provided in 'maternal_dates.csv' and the example individual-level residential history are provided in 'residential_history.csv'. The script relies on the publicly available CanOSSEM estimates of daily total PM2.5 (https://www.sfu.ca/fhs/canossem/about-canossem.html), provided in the 'Daily PM25 Data' folder. 
