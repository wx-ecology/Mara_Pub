[![DOI](https://zenodo.org/badge/740960158.svg)](https://zenodo.org/doi/10.5281/zenodo.10688128)

# Mara_Pub
Code and data to replicate the Mara analysis (DOI to be updated)

## Rethinking livestock encroachment at a protected area boundary

## **Data**

- 1 `mara_animal_compiled.rds` _dataset with animal occurrence by species and covariates used for Hmsc and EcoCopula models to test occurrence correlation and spatial displacement_
- 2 `mara_resource_cond.rds` _dataset with animal occurrence and resource conditions used for spaMM models to test resource compression_
- 3 `mara_gee_rain_20yr.csv` _TerraClimate precipitation history of the study area_
- 4 `Hmsc_VP_R1.RDS` _Hmsc variance partitioning output_
- 5 `Hmsc_prediction_R1.RDS` _Hmcs model prediction output_
- 6 `Hmsc_network_95_R1.RDS` _Hmsc species association network output_
- 7 `ecoCopula_partcor_network.RDS` _ecoCopula species association network output_


## **Code**

- 1 `1_HMSC` _code for running HMSC model, model validation, and model prediction_
- 2 `2_ecoCopula` _code for running ecocopula model_
- 3 `3_spaMM_table1` _code for spaMM models for testing resource compression_
- 4 `99_fig1_animal_stack_bars` _creating animal occurrence stack bar figure (figure 1)_
- 5 `99_fig3_figS3_spp_association_networks_and_matrix` _creating species association network and matrix figures (figure 3 and figure S3)_
- 6 `99_figs_HMSC_prediction_occurence_dist_to_bd` _creating HMSC prediction figures (figure 1 and figure S6)_
- 7 `99_figS5_HMSC_Variance_Partitioning` _creating HMSC variance partitioning figure (figure S5)_
- 8 `99_tableS1_tableS2_figS1_descriptive_stats` _creating table S1, table S2, and figure S1 regarding descriptive statistics on animal co-occurrence_


