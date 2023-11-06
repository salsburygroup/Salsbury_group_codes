Introduction
This repository contains the code and data for a molecular dynamics simulation study focusing on the influence of Thrombomodulin (TM) binding on thrombin. TM is a critical modulator of thrombin's function, switching its role from procoagulant to anticoagulant. We employed eight 1 μs all-atom simulations for the free thrombin, thrombin-TM56 complex, and thrombin-TM456 complex, respectively, using the CHARMM36 force field, significantly extending the simulation duration compared to previous works. The study reveals important insights into conformational changes and correlated movements within thrombin upon TM binding, enriching our understanding of blood coagulation and facilitating future research in anticoagulant therapies and oncology.

Data Availability
The trajectory data is too large to be uploaded to GitHub but is available upon request.

Repository Structure
Prepare/: Contains all input files for MD simulations
python/: Contains all Python scripts.
TM_thrombin/, TM56_Na0.125/, TM456_Na0.125/: Directories containing analyses for specific complexes.
pca/: PCA analyses files.
corr/: Correlation matrix files.
hbond/ and hbond_separate/: Directories for hydrogen bond data.
logistic_regression_Hbond/: Contains R scripts for logistic regression analysis on hydrogen bonds.


How to Run the Analyses

Initial Setup

Add the Analysis/ sub-directory, located within the python/ directory, to your Python path. To do this, execute the following command: 'export PYTHONPATH=$PYTHONPATH:/your_path/Analysis'

Edit the .sh files to correctly point to the Python scripts in the python/ directory.


Basic Analyses

Navigate to the specific analysis directory (e.g., TM_thrombin/, TM56_Na0.125/, or TM456_Na0.125/). Run Process1.sh to Process4.sh to conduct basic analyses such as RMSF and correlation matrix calculations.


Plot RMSF

To plot RMSF, run: 'python Plot_deviation.py'


Clustering Analysis

Execute ./clustering.sh. Visualize the clusters by running begin.sh in the clustering/ directory.


PCA Analysis

Use the Jupyter notebooks pca_calculate.ipynb, pca_plot.ipynb, and pca_plot_examine_wells.ipynb located in the pca/ directory. For visualization, execute begin.sh in the pca/visualization/ directory.


Correlation Matrix

Navigate to the corr/ directory and run Visualize.sh.


Hydrogen Bond Analysis

Execute the Python script hbond_analysis_resid.py to extract hydrogen bonds from the trajectory. Use select_Hbonds.ipynb to filter hydrogen bonds. For logistic regression analyses, navigate to the logistic_regression_Hbond/ directory and use the R scripts logistic_regression_Hbond.Rmd and logistic_regression_Hbond3.Rmd.
