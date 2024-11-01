# CCDC
Algorithm developed for Continuous Change Detection and Classification (CCDC) of land cover using all available Landsat data. Please contact Zhe Zhu (zhe@uconn.edu) at Department of Natural Resources and the Environment, University of Connecticut if you have any questions.

CCDC has been implemented on [Google Earth Engine (GEE)](https://developers.google.com/earth-engine/apidocs/ee-algorithms-temporalsegmentation-ccdc) and the related tools are at [this link](https://gee-ccdc-tools.readthedocs.io/en/latest).

The Classification software is not provided at the moment, as it requires training data to run the software.

CCDC Assistor 1.02 is also available at [here](https://drive.google.com/drive/folders/1iZmKlSNjJtb6DkinyOiPJtfT74YCE_eF), which is a user interface tool for assisting in data preparation and map extraction for CCDC (more functions are on the way).

Note that the output from CCDC will be thousands of Matlab files that contains all sorts of information for each time series models as follows: 

1. "t_start": when the time series model gets started

2. "t_end": when the time series model gets ended
 
3. "t_break": when the first break (change) is observed

4. "coefs": the coefficients for each time series model for each spectral band. Rows are a0 c1 a1 b1 a2 b2 a3 b3 (Zhu and Woodcock, 2014, 2015). Columns refer to Blue, Green, Red, NIR, SWIR1, SWIR2, and BT (Brightness Temperature) bands, respectively.

5. "rmse": the RMSE for each time series model for each spectral band

6. "pos": the position of each time series model (location)
 
7. "change_prob": the probability of a pixel that have undergone change (between 0 and 100)
 
8. "num_obs": the number of "good" observations used for model estimation

9. "category": the quality of the model estimation (what model is used, what process is used)
 
10. "magnitude": the magnitude of change (difference between model prediction and observation for each spectral band)

You need to extract those information from thousand of Matlab file to generate change maps or used as input for change detection. 

If you want to run on N cores, you will need to write script to submit job to each individual core by CCDC_ChangeARD13_01 i n (i=1,2,3...n; where n is the total number of cores, and i is which core to run the current job). CCDC is extremly computational expensive. Please use as many cores as you can on your Linux clusters. The computing time for one line of ARD data (5,000 pixels) takes around 1 hours for 1 core (CCDC process line-by-line). CCDC default parameters are 0.99 change probability, 6 consecutive observations, and a maximum of 8 coefficients for time series models. If you want to specify your parameters, you just need to create a .txt file named 'CCDC_Parameters.txt' within the images folder, in which the first variable specify change probability, the second specify number of consecutive days, and the last variable is the maximum number of coefficients (can be 4, 6, or 8), such as 0.95 5 6. 

Please cite the following papers

paper 1: Zhu, Z. and Woodcock, C. E., Continuous monitoring of forest disturbance using all available Landsat imagery, Remote Sensing of Environment (2012), doi:10.1016/j.rse.2011.10.030.(paper for CMFDA algorithm = CCDC version 1.0)

paper 2: Zhu, Z. and Woodcock, C. E., Continuous change detection and classification of land cover using all available Landsat data, Remote Sensing of Environment (2014), doi.org/10.1016/j.rse.2014.01.011.(paper for CCDC version 7.3)

paper 3: Zhu, Z., Woodcock, C. E., Holden, C., and Yang, Z., Generating synthetic Landsat images based on all available Landsat data: Predicting Landsat surface reflectance at any given time, Remote Sensing of Environment (2015), doi.org/10.1016/j.rse.2015.02.009.(paper for CCDC version 11.4)

This algorithm has been applied to many parts of the world and you can see all located where it has been applied [here](https://github.com/bullocke/Landsat-Database/blob/master/PRmap.geojson)  

You can download the PPT with GIF images that explain the CCDC algorithm at this [link](https://www.dropbox.com/s/1jzfte8mjy4qzzr/CCDC_algorithm_intro.pptx?dl=0)
