# DielLeafModeling
A repository housing scripts that can be used to generate C3, CAM and C4 diel leaf models using the PlantCoreMetabolism core model

### Links to google colab to test diel FBA leaf models
[Adding Arabidopsis Gene-Protein-Reaction Associations to the PlantCoreMetabolism model](https://colab.research.google.com/drive/1ziN3L51226-YH9JL-WoOVuerPDr-vasN?usp=sharing)  
[Generate C3 and C4 model from PlantCoreMetabolism](https://colab.research.google.com/drive/1ZCU4f5AkWnWa7IaWoS0czV7hsgNpbfSr?usp=sharing)  
[Model C3 leaf metabolism](https://colab.research.google.com/drive/1cWz1xdf4a2G19CyOA2HV03jqyw5MaWR7?usp=drive_link)  
[Tailoring diel leaf C3 model for any particular plant species](https://colab.research.google.com/drive/1x9n5cJV7wvIc1gdQHRV0iylhtquGDQKa?usp=sharing)  
[Model CAM leaf metabolism](https://colab.research.google.com/drive/1YLDYLRVifLH0N1vt8qFrL7EgpFInGhwg?usp=sharing)  
[Model C4 leaf metabolism](https://colab.research.google.com/drive/17ztwZvQ7urK_zM2SZ7pHNRdCxeMfocjv?usp=sharing)  
[Case study and flux visualization](https://colab.research.google.com/drive/1vXCq1sA72_pQM_V1tugicJUqPi-nGn2j?usp=sharing)  
[Application: studying how daylength affects day-time starch accumulation rate in C3 leaves](https://colab.research.google.com/drive/1F1v2fpFF_15e_zip4YzfziIQkYY58mV6?usp=sharing)  
[Application: studying how choice of biomass composition affects fluxes through central metabolism](https://colab.research.google.com/drive/14IfVlYSanbiivlJ5lc_iJ2VMV_tqHTYv?usp=sharing)  
[Application: studying the effect of photorespiratory bypasses in C3 leaves](https://drive.google.com/file/d/1FPcIwfuXtOt90k0WqVeBNmePMGkK2xt8/view?usp=sharing)  

## Repository contents
- Functions.py - a file containing all user-defined functions required in the study.The file contains functions used in the study to describe different case uses. Diel leaf model     
  construction required a function to make the PlantCoreMetabolism model to a day-night C3plant condition. This function can be used in any other leaf metabolic mode to make it diel C3. We 
  created another function for creating a diel C4 model, as mentioned in the use cases part of the manuscript. A general function was created to visualise model results and represented as 
  generateFluxMap. We have updated the pH of the model compartments and included fractionally charged forms of metabolite using two different functions, convertToClassicalModel and 
  convertToFractionalCharges.
### Applications
- Studying_how_starch_accumulation_changes_with_daylength_in_C3_leaves.ipynb - an IPython notebook that demonstrates how daylength can be varied in a 2-phase diel leaf model, and how such 
  models can be used to study the effect of daylength on metabolism
- Studying_the_effect_of_Photorespiratory_bypass_C3.ipynb - an IPython notebook which uses different alternate photorespiratory pathways for analysing the photosynthetic efficiency
- Studying_the_effect_of_biomass_compositions_on_leaf.ipynb - an IPython notebook to explore different biomass compositions in C3 model.
- Tailloring_diel_leaf_model_for__specific_plant_species.ipynb- an Ipython notebook to make species specific plant model
### Case Studies
- Scrip_for_generating_all_models_in_sbml.ipynb - an IPython notebook used to generate C3 and CAM models
- dielC3_script.ipynb - an IPython notebook that uses the C3 sbml file to model C3 metabolism
- dielCAM_script.ipynb - an IPython notebook that uses the C3 sbml file to model CAM metabolism
- dielC4_script.ipynb - an IPython notebook that uses the C4 sbml file to model C4 metabolism
- NocturnalCitrateInC3Leaves.ipynb - an IPython notebook that uses the C3 sbml to describe the role of nocturnal citrate and its visualization using fluxmap
### Data
- PlantCoreMetabolism2.cys - manually organized Cytoscape network of the core model with styles configured to visualize fluxes from 2-phase diel FBA model
- Adding_Arabidopsis_gene_associations_to_PlanCoreMetabolism_models.ipynb - an IPython notebook that demonstrates how gene associations can be added to the PlantCoreMetabolism models
- Uniprotkb_AND_model_organism_3702_2025_05_13.tsv - a tsv file with compiled information of EC numbers and subcellular localizations of genes.
- MetaboliteChargedStates.xls - excel containing fractional charged form of important metabolites
### Models
- PlantCoreMetabolism_v2_0_0.xml - The PlantCoreMetabolism version 2.0.0 model in SBML format
- C3_model.sbml - diel leaf model constrained for C3 photosynthesis
- C4_model.sbml - diel leaf meodel constrained for C4 photosynthesis

  



