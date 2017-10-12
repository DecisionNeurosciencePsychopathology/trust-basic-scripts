# trust_basic_scripts: 
1. trustbehavior.m: load in the data file from e-prime to create an individual .mat file.
2. trust_group.m: grouping the data from multiple subjects to save in spss format file.
3. trust_makeregressor_group.m: making regresssor files from .mat subject files.
4. asterisk.m: inserting asterisks into the existing regressor files (necessary to indicate multiple runs).
5. plot_trust_data.m: plotting behavioral data from trust.
6. subsample_demo.m: combines excel demographics with the matlab file of group data.
7. trustsurvey.m: load the survey data for a single subject and save as .mat file
8. trust_group_survey.m: collects the survey data (.mat) from all subjects and saves in a table w/ means
9. trust_combsurveyregs.m: after reading in survey data for new participants, combines their survey data w/ .mat subject reg files.
10. trust_survey_visualize.m: helps to visualize survey data.


helper functions:
1. read_in_trust.m: reads in e-prime file, used by trustbehavior.m
