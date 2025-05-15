# DielLeafModeling
A repository housing scripts that can be used to generate C3, CAM and C4 diel leaf models using the PlantCoreMetabolism core model

### Links to google colab to test diel FBA leaf models
[Adding Arabidopsis Gene-Protein-Reaction Associations to the PlantCoreMetabolism model](https://colab.research.google.com/drive/1ziN3L51226-YH9JL-WoOVuerPDr-vasN?usp=sharing)  
[Generate C3 and C4 model from PlantCoreMetabolism](https://colab.research.google.com/drive/1ZCU4f5AkWnWa7IaWoS0czV7hsgNpbfSr?usp=sharing)  
[Model C3 leaf metabolism](https://colab.research.google.com/drive/1cWz1xdf4a2G19CyOA2HV03jqyw5MaWR7?usp=drive_link)  
[Model CAM leaf metabolism](https://colab.research.google.com/drive/1YLDYLRVifLH0N1vt8qFrL7EgpFInGhwg?usp=sharing)  
[Model C4 leaf metabolism](https://colab.research.google.com/drive/17ztwZvQ7urK_zM2SZ7pHNRdCxeMfocjv?usp=sharing)  
[Case study and flux visualization](https://colab.research.google.com/drive/1vXCq1sA72_pQM_V1tugicJUqPi-nGn2j?usp=sharing)  
[Case study: studying how daylength affects day-time starch accumulation rate in C3 leaves](https://colab.research.google.com/drive/1F1v2fpFF_15e_zip4YzfziIQkYY58mV6?usp=sharing)

### Repository contents
- Functions.py - a file containing all user-defined functions required in the study
- PlantCoreMetabolism_v2_0_0.xml - The PlantCoreMetabolism version 2.0.0 model in SBML format
- C3_model.sbml - diel leaf model constrained for C3 photosynthesis
- C4_model.sbml - diel leaf meodel constrained for C4 photosynthesis
- PlantCoreMetabolism2.cys - manually organized Cytoscape network of the core model with styles configured to visualize fluxes from 2-phase diel FBA model
- Generating_all_models.ipynb - an IPython notebook used to generate C3 and CAM models
- dielC3_script.ipynb - an IPython notebook that uses the C3 sbml file to model C3 metabolism
- dielCAM_script.ipynb - an IPython notebook that uses the C3 sbml file to model CAM metabolism
- dielC4_script.ipynb - an IPython notebook that uses the C4 sbml file to model C4 metabolism
- Studying_how_stach_accumulation_changes_with_daylength_in_C3_leaves.ipynb - an IPython notebook that demonstrates how daylength can be varied in a 2-phase diel leaf model, and how such models can be used to study the effect of daylength on metabolism
- Adding_Arabidopsis_gene_associations_to_PlanCoreMetabolism_models.ipynb - an IPython notebook that demonstrates how gene associations can be added to the PlantCoreMetabolism models
  



