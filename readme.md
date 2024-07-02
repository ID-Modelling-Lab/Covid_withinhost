# Covid_withinhost
This repository contains the code for the within-host model of SARS-CoV-2 described in ['Differences in virus and immune dynamics for SARS-CoV-2 Delta and Omicron infections by age and vaccination histories'](https://doi.org/10.1186/s12879-024-09572-x).

We developed a mathematical model describing SARS-CoV-2 infection dynamics, with a simple representation of immunity. The model is fitted to URT viral load measurements from Delta and Omicron cases in Singapore, of whom the majority only had one nasopharyngeal swab measurement. By grouping individuals by age and vaccination history (vaccine brand, number of doses, time since last vaccination), we could recreate similar trends in URT virus dynamics observed in past within-host modelling studies despite the lack of longitudinal patient data.

A basic summary of the files in this repository follows.

## Within-host models
There are 4 model types described in the manuscript, of which fitting results of model type 2 were analysed and presented. The codes required to run the models are available in the 'main models' file.
Each model comprises two parts:
1. Stan script
2. R-script for model fitting

The data used in this analysis is not provided due to confidentiality issues. ['VL_standata.xlsx'](https://github.com/ID-Modelling-Lab/Covid_withinhost/blob/main/main_models/VL_standata.xlsx) is a mock dataset, available in the 'main models' file, with format similar to our processed data that is ready for Stan.

Results of model fitting were saved as RDS files, but could not be uploaded due to large file size. These are available upon email request (maxine@u.nus.edu).

## Figures in the manuscript
- R-scripts to generate the manuscript figures: ['Fig_2.R'](https://github.com/ID-Modelling-Lab/Covid_withinhost/blob/main/Fig_2.R), ['Fig_3-6.R'](https://github.com/ID-Modelling-Lab/Covid_withinhost/blob/main/Fig_3-6.R).
- R-scripts to generate supplementary figures: ['SupFig_3-5.R'](https://github.com/ID-Modelling-Lab/Covid_withinhost/blob/main/SupFig_3-5.R), ['SupFig_7.R'](https://github.com/ID-Modelling-Lab/Covid_withinhost/blob/main/SupFig_7.R), ['SupFig_10.R'](https://github.com/ID-Modelling-Lab/Covid_withinhost/blob/main/SupFig_10.R), ['SupFig_13-14.R'](https://github.com/ID-Modelling-Lab/Covid_withinhost/blob/main/SupFig_13-14.R).
- RData files necessary for figures in main text ['labels.RData'](https://github.com/ID-Modelling-Lab/Covid_withinhost/blob/main/labels.RData) and supplement ['labels_1sampleonly.RData'](https://github.com/ID-Modelling-Lab/Covid_withinhost/blob/main/labels_1sampleonly.RData), ['labels_1sampperpat.RData'](https://github.com/ID-Modelling-Lab/Covid_withinhost/blob/main/labels_1sampperpat.RData), ['labels_swapIP.RData'](https://github.com/ID-Modelling-Lab/Covid_withinhost/blob/main/labels_swapIP.RData): These files are necessary to provide labelling 

