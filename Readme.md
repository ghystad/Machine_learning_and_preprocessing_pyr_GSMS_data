
# The R-code is created for the paper

## A Robust, Agnostic Molecular Biosignature Based on Machine Learning. 
H. James Cleaves II^1,2,3, Grethe Hystad^4, Anirudh Prabhu^1, Michael L. Wong ^1,5, George D. Cody^1, Sophia Economon^6, and Robert M. Hazen^1

*1. Earth and Planets Laboratory, Carnegie Science, Washington, DC 20015, USA*
*2. Earth Life Science Institute, Tokyo Institute of Technology, Tokyo, Japan*
*3. Blue Marble Space Institute for Science, Seattle, WA 98104, USA*
*4. Mathematics and Statistics, Purdue University Northwest, Hammond, IN 46323, USA*
*5. NHFP Sagan Fellow, NASA Hubble Fellowship Program, Space Telescope Science Institute, Baltimore, MD 21218, USA*
*6. Department of Earth and Planetary Sciences, Johns Hopkins University, Baltimore, MD 21218, USA*

Three-dimensional (scan number /mass-to-charge ratio/intensity) data from biotic and abiotic samples are obtained by pyrolysis-gas chromatography mass spectrometry. The R-code created is for preprocessing these data and to use machine learning to predict whether a sample is biotic or abiotic (Preprocessing_training_and_test_data_random_forest).
Nested resampling is used to obtained an estimate for the prediction performance of the model (Nested_resampling_random_forest).

The resulting peaks from the preprocessing steps are aligned based on the paper by 
Tibshirani, R., et al. (2004) 
Sample classification from proteing mass spectrometry by peak probability contrasts.
Bioinformatis, 20(17):3034-44
https://academic.oup.com/bioinformatics/article/20/17/3034/186323

The pyr-GC-MS data, named "Cleavesetal.pyrGCMSData.zip", can be found at https://accounts.osf.io/login?service=https://osf.io/embh8/

The application is released under GNU GPL version 3 license, noting that the following R-libraries are licensed as follows:

dplyr: MIT
MALDIquant: GPL ≥ 3
caret: GPL-2, GPL-3 (expanded from GPL-2)
mlr3: LGPL-3
mlr3verse: LGPL-3
mlr3learners: LGPL-3
mlr3tuning: LGPL-3
data.table: MPL-2.0
ggplot2: MIT
