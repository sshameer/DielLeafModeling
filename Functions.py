#################################################################################
#The file contains functions used in the study to describe different case uses. 
#Diel leaf model construction required a function to make the PlantCoreMetabolism-
#model to a day-night C3plant condition. This function can be used in any other -
#leaf metabolic model to make it diel. We created another function for creating a-
#diel C4 model, as mentioned in the use cases part of the manuscript. A general -
#function was created to visualise model results and represented as generateFluxMap.
#We have updated the pH of the model compartments and included fractionally charged- 
#forms of metabolites using two different functions, convertToClassicalModel and  -
#convertToFractionalCharges, at the bottom of this documentation file.
##################################################################################

#1. This function contains the basic constraints for making a diel C3 leaf model
    #It duplicates the model for the day-night phase and adds the tags as 1 and 2 on 
    #compartments, metabolites and reactions. After that, merge them. Then customizable
    #constraints were included for its simplicity and flexibility. 
########################################################
#This function was used to set up a C3 leaf diel model #
########################################################
def setupC3DielModel(core_model,transferMets="",starch_sucrose_ratio=None):
    '''
    This function can be used to generate a fully constrained diel C3 leaf model
    from a core model.
    Inputs: 1) a cobra model 2) a list of metabolites allowed to accumulate 3)
    starch to sucrose accumulation rate ratio
    Outputs: a fully constrained diel C3 leaf model
    '''
    from cobra.core import Metabolite, Reaction
    import re

    #create two copies of model elements for day and night
    tempCompDict = dict()
    for comp in core_model.compartments:
        tempCompDict[comp+"1"] = core_model.compartments[comp]+" day"
        tempCompDict[comp+"2"] = core_model.compartments[comp]+" night"

    cobra_model2 = core_model.copy()
    for met in cobra_model2.metabolites:
        met.id = met.id+"1"
        met.name = met.name+" Day"
        met.compartment = met.compartment+"1"
    for rxn in cobra_model2.reactions:
        rxn.id = rxn.id+"1"
        rxn.name = rxn.name+" Day"

    cobra_model3 = core_model.copy()
    for met in cobra_model3.metabolites:
        met.id = met.id+"2"
        met.name = met.name+" Night"
        met.compartment = met.compartment+"2"
    for rxn in cobra_model3.reactions:
        rxn.id = rxn.id+"2"
        rxn.name = rxn.name+" Night"

    #merge the day and night model
    cobra_model = cobra_model2.merge(cobra_model3)
    for met in cobra_model3.metabolites:
        if not cobra_model.metabolites.__contains__(met.id):
            cobra_model.add_metabolites(met.copy())
    for comp in cobra_model.compartments:
        cobra_model.compartments[comp] = tempCompDict[comp]

    met1 = Metabolite("X_Phloem_contribution_t1",name="Phloem output during the day",compartment="b1")
    cobra_model.reactions.get_by_id("Phloem_output_tx1").add_metabolites({met1:1})
    met2 = Metabolite("X_Phloem_contribution_t2",name="Phloem output during at night",compartment="b1")
    cobra_model.reactions.get_by_id("Phloem_output_tx2").add_metabolites({met2:1})

    rxn = Reaction("diel_biomass")
    rxn.add_metabolites({met1:-3,met2:-1})
    rxn.lower_bound = 0
    rxn.upper_bound = 1000
    cobra_model.add_reactions([rxn,])

    #Adding reactions to allow for day-night metabolite accumulations
    if transferMets!="":
        tmfile = open(transferMets,"r")
        tmset=set()
        for line in tmfile:
            tmset.add(line.replace("\n",""))
    else:
        tmset=set(["STARCH_p","SUCROSE_v","MAL_v","aMAL_v","NITRATE_v","CIT_v",
        "aCIT_v","GLN_v","ASN_v","SER_v","GLN_v","GLY_v","THR_v","L_ALPHA_ALANINE_v",
        "4_AMINO_BUTYRATE_v","VAL_v","ILE_v","PHE_v","LEU_v","LYS_v","ARG_v",
        "L_ASPARTATE_v","GLT_v","HIS_v","bHIS_v","MET_v","PRO_v","TRP_v","TYR_v",
        "CYS_v","FRUCTAN_v","AMMONIUM_v","PROTON_v"])

    for met in tmset:
        if met == "AMMONIUM_v" or met=="FRUCTAN_v":
            continue
        tempRxn = Reaction(met+"_dielTransfer")
        tempRxn.add_metabolites({cobra_model.metabolites.get_by_id(met+"1"):-1,cobra_model.metabolites.get_by_id(met+"2"):1})
        tempRxn.lower_bound=-1000
        if not ((met == "STARCH_p") or (met == "SUCROSE_v") or (met == "MAL_v") or (met == "aMAL_v") or (met == "NITRATE_v") or (met == "CIT_v") or (met == "aCIT_v") or (met == "PROTON_v")):
            tempRxn.lower_bound=0
        tempRxn.upper_bound=1000
        cobra_model.add_reactions([tempRxn,])

    fractionMets=dict()
    for rxn in cobra_model.reactions:
        for met in rxn.metabolites.keys():
            prefix=""
            a=re.search("^a{1,3}",met.id)
            anion=""
            if a:
                anion=a.group(0)
                prefix=anion
            b=re.search("^b{1,3}",met.id)
            basic=""
            if b:
                basic=b.group(0)
                prefix=basic
            if ((not prefix == "") and met.compartment == "v1"):
                fractionMets[met]=prefix

    temp=cobra_model.copy()
    for met in fractionMets.keys():
        for rxn in met.reactions:
            if rxn.id.__contains__("_dielTransfer"):
                continue
            else:
                mainMet = met.id[len(fractionMets[met]):]
                coeff1 = temp.reactions.get_by_id(rxn.id).metabolites.get(temp.metabolites.get_by_id(mainMet))
                coeff2 = temp.reactions.get_by_id(rxn.id).metabolites.get(temp.metabolites.get_by_id(met.id))
                if not coeff1:
                    coeff1=0
                if not coeff2:
                    coeff2=0
                total = coeff1 + coeff2
                coeff1 = float(coeff1)/total
                coeff2 = float(coeff2)/total
                if cobra_model.reactions.has_id(met.id[0:len(met.id)-1]+"_dielTransfer"):
                    ub = temp.reactions.get_by_id(met.id[0:len(met.id)-1]+"_dielTransfer").upper_bound
                    lb = temp.reactions.get_by_id(met.id[0:len(met.id)-1]+"_dielTransfer").lower_bound
                    temp.reactions.get_by_id(met.id[0:len(met.id)-1]+"_dielTransfer").remove_from_model()
                    temp.reactions.get_by_id(mainMet[0:len(mainMet)-1]+"_dielTransfer").remove_from_model()
                    Reac = Reaction(mainMet[0:len(mainMet)-1]+"_dielTransfer",name=mainMet+"_dielTransfer")
                    Reac.add_metabolites({temp.metabolites.get_by_id(met.id[0:len(met.id)-1]+"1"):-coeff2,temp.metabolites.get_by_id(met.id[0:len(met.id)-1]+"2"):coeff2,temp.metabolites.get_by_id(mainMet[0:len(mainMet)-1]+"1"):-coeff1,temp.metabolites.get_by_id(mainMet[0:len(mainMet)-1]+"2"):coeff1})
                    Reac.lower_bound=lb
                    Reac.upper_bound=ub
                    temp.add_reactions([Reac,])
                break
    ####ADD CONSTRAINTS TO MODEL####
    cobra_model = temp.copy()

    #objective function
    cobra_model.reactions.get_by_id("diel_biomass").objective_coefficient=1
    #Leaves - light
    cobra_model.reactions.get_by_id("Sucrose_tx1").lower_bound=0
    cobra_model.reactions.get_by_id("Sucrose_tx1").upper_bound=0
    cobra_model.reactions.get_by_id("GLC_tx1").lower_bound=0
    cobra_model.reactions.get_by_id("GLC_tx1").upper_bound=0
    cobra_model.reactions.get_by_id("NH4_tx1").lower_bound=0
    cobra_model.reactions.get_by_id("NH4_tx1").upper_bound=0
    #Leaves - dark
    cobra_model.reactions.get_by_id("Sucrose_tx2").lower_bound=0
    cobra_model.reactions.get_by_id("Sucrose_tx2").upper_bound=0
    cobra_model.reactions.get_by_id("GLC_tx2").lower_bound=0
    cobra_model.reactions.get_by_id("GLC_tx2").upper_bound=0
    cobra_model.reactions.get_by_id("Photon_tx2").lower_bound=0
    cobra_model.reactions.get_by_id("Photon_tx2").upper_bound=0
    cobra_model.reactions.get_by_id("NH4_tx2").lower_bound=0
    cobra_model.reactions.get_by_id("NH4_tx2").upper_bound=0

    #Turn off PTOX
    cobra_model.reactions.get_by_id("Plastoquinol_Oxidase_p1").lower_bound=0
    cobra_model.reactions.get_by_id("Plastoquinol_Oxidase_p1").upper_bound=0

    #nitrate uptake constrain
    Nitrate_balance = Metabolite("Nitrate_bal_c", name = "Weights to balance nitrate uptake", compartment = "c1")
    cobra_model.reactions.get_by_id("Nitrate_ec1").add_metabolites({Nitrate_balance:-2})
    cobra_model.reactions.get_by_id("Nitrate_ec2").add_metabolites({Nitrate_balance:3})

    #Rubisco balance
    Rubisco_balance = Metabolite("rubisco_bal_p1", name = "Weights to balance RuBP carboxygenase oxygenase balance", compartment = "p1")
    cobra_model.reactions.get_by_id("RXN_961_p1").add_metabolites({Rubisco_balance:3})
    cobra_model.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p1").add_metabolites({Rubisco_balance:-1})

    #generic ATPase and NADPH oxidase
    Maintenance_constraint = Metabolite("ATPase_NADPHoxidase_constraint_c1",name =  "ATPase_NADPHoxidase_constraint_c1", compartment = "c1")
    Maintenance_constraint2 = Metabolite("ATPase_NADPHoxidase_constraint_c2",name =  "ATPase_NADPHoxidase_constraint_c2", compartment = "c2")
    Maintenance_constraint3 = Metabolite("Light_dark_maintainence_constraint",name =  "Light_dark_maintainence_constraint", compartment = "c1")
    cobra_model.reactions.get_by_id("ATPase_tx1").add_metabolites({Maintenance_constraint:1,Maintenance_constraint3:1})
    cobra_model.reactions.get_by_id("ATPase_tx2").add_metabolites({Maintenance_constraint2:1,Maintenance_constraint3:-1})
    cobra_model.reactions.get_by_id("NADPHoxc_tx1").add_metabolites({Maintenance_constraint:-3})
    cobra_model.reactions.get_by_id("NADPHoxc_tx2").add_metabolites({Maintenance_constraint2:-3})
    cobra_model.reactions.get_by_id("NADPHoxm_tx1").add_metabolites({Maintenance_constraint:-3})
    cobra_model.reactions.get_by_id("NADPHoxm_tx2").add_metabolites({Maintenance_constraint2:-3})
    cobra_model.reactions.get_by_id("NADPHoxp_tx1").add_metabolites({Maintenance_constraint:-3})
    cobra_model.reactions.get_by_id("NADPHoxp_tx2").add_metabolites({Maintenance_constraint2:-3})

    ##constrain sucrose and starch storage
    if starch_sucrose_ratio is not None:
        Sucorse_starch_balance = Metabolite("sucrose_starch_bal_c", name = "Weights to balance sucrose-starch uptake", compartment = "c1")
        cobra_model.reactions.get_by_id("SUCROSE_v_dielTransfer").add_metabolites({Sucorse_starch_balance:-1*starch_sucrose_ratio})
        cobra_model.reactions.get_by_id("STARCH_p_dielTransfer").add_metabolites({Sucorse_starch_balance:1})

    #Plastid enolase was not detected in Arabidopsis mesophyll tissue
    cobra_model.reactions.get_by_id("2PGADEHYDRAT_RXN_p1").lower_bound=0
    cobra_model.reactions.get_by_id("2PGADEHYDRAT_RXN_p1").upper_bound=0
    cobra_model.reactions.get_by_id("2PGADEHYDRAT_RXN_p2").lower_bound=0
    cobra_model.reactions.get_by_id("2PGADEHYDRAT_RXN_p2").upper_bound=0

    #Setting chloroplastic NADPH dehydrogenase to 0  ((Yamamoto et al., 2011)
    cobra_model.reactions.get_by_id("NADPH_Dehydrogenase_p1").lower_bound=0
    cobra_model.reactions.get_by_id("NADPH_Dehydrogenase_p1").upper_bound=0
    cobra_model.reactions.get_by_id("NADPH_Dehydrogenase_p2").lower_bound=0
    cobra_model.reactions.get_by_id("NADPH_Dehydrogenase_p2").upper_bound=0

    #ATP_ADP_Pi constrained to 0 because while there is evidence for its existance, it does not carry high flux
    cobra_model.reactions.get_by_id("ATP_ADP_Pi_pc1").lower_bound = 0
    cobra_model.reactions.get_by_id("ATP_ADP_Pi_pc1").upper_bound = 0
    cobra_model.reactions.get_by_id("ATP_ADP_Pi_pc2").lower_bound = 0
    cobra_model.reactions.get_by_id("ATP_ADP_Pi_pc2").upper_bound = 0

    return cobra_model


    #####################################################################################
#2. This function contains the basic constraints for making a diel C4 leaf model
    #It duplicates the model for the day-night phase and adds the tags as 1 and 2 on 
    #compartments, metabolites and reactions. After that, merge them. Then customizable
    #constraints were included for its simplicity and flexibility.

########################################################
#This function was used to set up a C4 leaf diel model #
########################################################
def setupC4DielModel(core_model,transferMets="",M_BS_transferMets="",starch_sucrose_ratio=None):
    '''
    This function can be used to generate a fully constrained diel C3 leaf model
    from a core model.
    Inputs: 1) a cobra model 2) a list of metabolites allowed to accumulate 3)
    starch to sucrose accumulation rate ratio
    Outputs: a fully constrained diel C3 leaf model
    '''
    from cobra.core import Metabolite, Reaction
    import re

    M_model = setupC3DielModel(core_model,transferMets,starch_sucrose_ratio)
    for met in M_model.metabolites:
      if met.id[-1]=="1":
        met.name = met.name.replace("Day","M Day")
      elif met.id[-1]=="2":
        met.name = met.name.replace("Night","M Night")
    for rxn in M_model.reactions:
      if rxn.id[-1]=="1":
        rxn.name = rxn.name.replace("Day","M Day")
      elif rxn.id[-1]=="2":
        rxn.name = rxn.name.replace("Day","M Night")
      elif "dielTransfer" in rxn.id:
        rxn.id = rxn.id.replace("dielTransfer","dielTransfer_M")
    M_model.reactions.diel_biomass.remove_from_model()
    M_model.reactions.Phloem_output_tx1.remove_from_model()
    M_model.reactions.Phloem_output_tx2.remove_from_model()

    BS_model = setupC3DielModel(core_model,transferMets,starch_sucrose_ratio)
    for met in BS_model.metabolites:
      if met.id[-1]=="1":
        met.id = met.id[0:-1]+"3"
        # print(met.id)
        met.name = met.name.replace("Day","BS Day")
        met.compartment=met.compartment[0:-1]+"3"
        # print(met.id)
      elif met.id[-1]=="2":
        # print(met.id)
        met.id = met.id[0:-1]+"4"
        met.name = met.name.replace("Night","BS Night")
        met.compartment=met.compartment[0:-1]+"4"
        # print(met.id)
    for rxn in BS_model.reactions:
      if rxn.id[-1]=="1":
        # print(rxn.id)
        rxn.id = rxn.id[0:-1]+"3"
        rxn.name = rxn.name.replace("Day","BS Day")
        # print(rxn.id)
      elif rxn.id[-1]=="2":
        rxn.id = rxn.id[0:-1]+"4"
        rxn.name = rxn.name.replace("Day","BS Night")
      elif "dielTransfer" in rxn.id:
        rxn.id = rxn.id.replace("dielTransfer","dielTransfer_BS")

    temp = M_model.merge(BS_model)

    if M_BS_transferMets!="":
        tmfile = open(M_BS_transferMets,"r")
        tmset=set()
        for line in tmfile:
            tmset.add(line.replace("\n",""))
    else:
        tmset=set(["MAL_c","PYRUVATE_c","SUCROSE_c","MALTOSE_p","ASN_c","GLT_c",
                   "L_ALPHA_ALANINE_c","LEU_c","MET_c","LYS_c","HIS_c","THR_c",
                   "VAL_c","PHE_c","4_AMINO_BUTYRATE_c","SER_c","ARG_c","GLN_c",
                   "GLY_c","L_ASPARTATE_c","ILE_c","PRO_c","CYS_c","TRP_c",
                   "TYR_c","G3P_c","GAP_c","Pi_c","aPi_c","PHOSPHO_ENOL_PYRUVATE_c"])

    for met in tmset:
        if met == "AMMONIUM_v" or met=="FRUCTAN_v":
            continue
        tempRxn = Reaction(met+"_MBStransfer_Day")
        tempRxn.add_metabolites({temp.metabolites.get_by_id(met+"1"):-1,temp.metabolites.get_by_id(met+"3"):1})
        tempRxn.lower_bound=-1000
        tempRxn.upper_bound=1000
        temp.add_reactions([tempRxn,])

    for met in tmset:
        if met == "AMMONIUM_v" or met=="FRUCTAN_v":
            continue
        tempRxn = Reaction(met+"_MBStransfer_Night")
        tempRxn.add_metabolites({temp.metabolites.get_by_id(met+"2"):-1,temp.metabolites.get_by_id(met+"4"):1})
        tempRxn.lower_bound=-1000
        tempRxn.upper_bound=1000
        temp.add_reactions([tempRxn,])

    ####ADD CONSTRAINTS TO MODEL####
    C4_model = temp.copy()
    C4_model.reactions.diel_biomass.remove_from_model()

    #setting the total output with 3:1 ratio
    reac = Reaction("Total_Phloem_output_tx")
    reac.name = "Total Phloem Output"
    C4_model.add_reactions([reac,])
    C4_model.reactions.get_by_id("Total_Phloem_output_tx").add_metabolites({"X_Phloem_contribution_t3" : -0.75, "X_Phloem_contribution_t4" : -0.25})
    reac.objective_coefficient = 1

    # Bundle sheath CO2 uptake constraint to 0
    C4_model.reactions.CO2_tx3.bounds=(0,0)
    C4_model.reactions.CO2_tx4.bounds=(0,0)

    #Rubisco balance for bundle sheath
    C4_model.metabolites.get_by_id("rubisco_bal_p3").remove_from_model()
    Rubisco_balance = Metabolite("rubisco_bal_p3", name = "Weights to balance RuBP carboxygenase oxygenase balance", compartment = "p3")
    C4_model.reactions.get_by_id("RXN_961_p3").add_metabolites({Rubisco_balance:29.2})
    C4_model.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p3").add_metabolites({Rubisco_balance:-1})

    return C4_model
    
#################################################################################
#3. This function helps to visualize the model results by generating Fluxmaps
#################################################################################

def generateFluxMap(cobra_model,solution, outfile,phases = 2):
    import cobra
    #solution = cobra.flux_analysis.parsimonious.optimize_minimal_flux(cobra_model)
    #solution = cobra.flux_analysis.parsimonious.pfba(cobra_model)          #If the previous line returns error comment it out and uncomment this line instead

    #open output file for writing
    f = open(outfile,"w");

    #use rxnSet to identify reaction that have already been processed
    rxnSet = set()

    mult=set()
    #Looping through all reactions in the model
    for rxn in cobra_model.reactions:
        #Get the ID
        RXN=rxn.id
        #declare a boolean variable multFlag to keep track of whether this reaction is present in multiple models
        multFlag=False

        #check if the reaction has already been processed before and if yes skip this run in the loop
        if(rxnSet.__contains__(RXN)):
            continue
        if rxn.id.__contains__("EX") or rxn.id.__contains__("Transfer"):
            multFlag=False
        #check if the reaction ends with one or two i.e it is present more than once in the model
        elif(["1","2","3","4","5","6","7","8","9"].__contains__(rxn.id[len(rxn.id)-1])):
            #change the id to without the suffix 1-9 and declare it as a reaction which has multiple instances
            RXN = rxn.id[0:len(rxn.id)-1]
            multFlag=True
        elif rxn.id[len(rxn.id)-2:] == "10":
            #change the id to without the suffix 10 and declare it as a reaction which has multiple instances
            RXN = rxn.id[0:len(rxn.id)-2]
            multFlag=True

        #if metabolite has multiple instances
        values = dict()
        status1 = dict()
        status2 = dict()
        if(multFlag):
            tempvalue = list()
            temp1 = list()
            temp2 = list()
            mult.add(RXN)
            #add the reaction we are about to process to the reactions processed list
            for i in range(1,phases+1):
                rxnSet.add(RXN+str(i))
                if(round(float(solution.fluxes.get(RXN+str(i)))*10000000) == 0):
                    tempvalue.append(0)
                    temp1.append("none")
                    temp2.append("none")
                elif(float(solution.fluxes.get(RXN+str(i)))*10000 > 0):
                    tempvalue.append(solution.fluxes.get(RXN+str(i)))
                    temp1.append("produced")
                    temp2.append("consumed")
                elif(float(solution.fluxes.get(RXN+str(i)))*10000 < 0):
                    tempvalue.append(solution.fluxes.get(RXN+str(i)))
                    temp1.append("consumed")
                    temp2.append("produced")
            values[RXN] = tempvalue
            status1[RXN] = temp1
            status2[RXN] = temp2

            #select 1 reaction so that we can identify the reactants and products which can be then used to generate the edge shared_name
            rxn=cobra_model.reactions.get_by_id(RXN+"1")

            for reac in rxn.reactants:
                REAC=reac.id
                if(REAC.__contains__("1")):
                    if(REAC.rindex("1")==len(REAC)-1):
                        REAC=REAC[0:len(REAC)-1]
                    f.write("R_"+RXN+" (reaction-reactant) M_"+REAC)
                    for i in range(1,phases+1):
                        f.write("\t"+str(values[RXN][i-1])+"\t"+str(status2[RXN][i-1]))
                    f.write("\n")
                if(RXN.__contains__("biomass")):
                    f.write("R_"+RXN+" (reaction-product)) M_"+REAC)
                    for i in range(1,phases+1):
                        f.write("\t"+str(values[RXN][i-1])+"\t"+str(status1[RXN][i-1]))
                    f.write("\n")
            for prod in rxn.products:
                PROD=prod.id
                if(PROD.__contains__("1")):
                    if(PROD.rindex("1")==len(PROD)-1):
                        PROD=PROD[0:len(PROD)-1]
                f.write("R_"+RXN+" (reaction-product) M_"+PROD)
                for i in range(1,phases+1):
                    f.write("\t"+str(values[RXN][i-1])+"\t"+str(status1[RXN][i-1]))
                f.write("\n")
            if(RXN.__contains__("biomass")):
                f.write("R_"+RXN+" (reaction-reactant) M_"+REAC)
                for i in range(1,phases+1):
                    f.write("\t"+str(values[RXN][i-1])+"\t"+str(status2[RXN][i-1]))
                f.write("\n")
        else:
            #add the reaction we are about to process to the reactions processed list
            rxnSet.add(RXN)
            if(round(float(solution.fluxes.get(rxn.id))*10000000) == 0):
                value = 0;
                status1= "none";
                status0= "none";
            elif(solution.fluxes.get(rxn.id)*10000 > 0):
                value = solution.fluxes.get(rxn.id)*1000;
                status1= "produced";
                status0= "consumed";
            elif(solution.fluxes.get(rxn.id)*10000 < 0):
                value = solution.fluxes.get(rxn.id)*1000;
                status1= "consumed";
                status0= "produced";

            for reac in rxn.reactants:
                REAC=reac.id
                if(REAC.__contains__("1")):
                    if(REAC.rindex("1")==len(REAC)-1):# or (met.id.rindex("2")==len(rxn.id)-1):
                        REAC=REAC[0:len(REAC)-1]
                f.write("R_%s (reaction-reactant) M_%s\t%s\t%s\t0\tnone\n" % (RXN,REAC,value,status0));
            for prod in rxn.products:
                PROD=prod.id
                if(PROD.__contains__("1")):
                    if(PROD.rindex("1")==len(PROD)-1):# or (met.id.rindex("2")==len(rxn.id)-1):
                        PROD=PROD[0:len(PROD)-1]
                f.write("R_%s (reaction-product) M_%s\t%s\t%s\t0\tnone\n" % (RXN,PROD,value,status1));

    f.close();
######################################################################################
#4.  This function converts a model which fractionally charged metabolites to a classical 
     #model by setting the pH to default in any particular compartment.
######################################################################################

def convertToClassicalModel(cobra_model2,comp="*",updateCharges=""):
    if not updateCharges == "":
        #override charges from file
        fin = open(updateCharges,"r")
        ChargeDict=dict()
        for line in fin:
            met=line.replace("\n","").split("\t")[0]
            met = met.replace("-","_")
            charge = line.replace("\n","").split("\t")[1]
            ChargeDict[met]=charge
        fin.close()

        temp=cobra_model2.copy()
        for met in temp.metabolites:
            if(met.compartment=="b"):
                cobra_model2.metabolites.get_by_id(met.id).remove_from_model()

        for met in cobra_model2.metabolites:
            tempMet=met.id
            if(["1","2","3","4","5","6","7","8","9"].__contains__(met.id[len(met.id)-1])):
                tempMet = met.id[0:len(met.id)-1]
            elif(met.id[len(met.id)-2:]=="10"):
                tempMet = met.id[0:len(met.id)-2]
            if(ChargeDict.keys().__contains__(tempMet)):
                met.charge = ChargeDict.get(tempMet)
            if met.charge is None:
                met.charge=0

    import re
    #List metabolites that were added
    fractionMets=set()
    for rxn in cobra_model2.reactions:
        for met in rxn.metabolites.keys():
            a=re.search("^a{1,3}",met.id)
            anion=""
            if a:
                anion=a.group(0)
            b=re.search("^b{1,3}",met.id)
            basic=""
            if b:
                basic=b.group(0)
            if (abs(rxn.metabolites.get(met)) % 1 > 0 and (not anion == "" or not basic== "")):
                if not comp == "*":
                    if met.compartment == comp:
                        fractionMets.add(met.id)
                else:
                    fractionMets.add(met.id)
    cobra_model = cobra_model2.copy()

    for met in fractionMets:
        for rxn in cobra_model.metabolites.get_by_id(met).reactions:
            a=re.search("^a{1,3}",met)
            anion=""
            if a:
                anion=a.group(0)
            b=re.search("^b{1,3}",met)
            basic=""
            if b:
                basic=b.group(0)
            prefix = anion
            if prefix == "":
                prefix = basic
                
            coeff1 = rxn.metabolites.get(cobra_model.metabolites.get_by_id(met[len(prefix):]))
            if coeff1 is None:
                  raise ValueError(f"Metabolite {met[len(prefix):]} not found in reaction {rxn.id}")
            rxn.add_metabolites({cobra_model.metabolites.get_by_id(met[len(prefix):]):-1*coeff1})
            coeff2 = rxn.metabolites.get(cobra_model.metabolites.get_by_id(met))
            rxn.add_metabolites({cobra_model.metabolites.get_by_id(met):-1*coeff2})
            Charge=(coeff1*float(cobra_model.metabolites.get_by_id(met[len(prefix):]).charge))+(coeff2*float(cobra_model.metabolites.get_by_id(met).charge))-((coeff1+coeff2)*float(cobra_model.metabolites.get_by_id(met[len(prefix):]).charge))
            if rxn.id.__contains__("_Transfer"):
                compSet = set()
                for m in rxn.metabolites:
                    compSet.add(m.compartment)
                this = cobra_model.metabolites.get_by_id(met).compartment
                rxn.add_metabolites({cobra_model.metabolites.get_by_id(met[len(prefix):]):coeff1+coeff2,cobra_model.metabolites.get_by_id("PROTON_"+cobra_model.metabolites.get_by_id(met).compartment):round(Charge,4)})
                other = compSet.difference(this)
                for c in other:
                    if cobra_model.metabolites.get_by_id("PROTON_"+c).reactions.__contains__(rxn):
                        prot_bal = rxn.metabolites.get(cobra_model.metabolites.get_by_id("PROTON_"+c)) + rxn.metabolites.get(cobra_model.metabolites.get_by_id("PROTON_"+this))
                        if prot_bal > 0:
                            tempCoeff = rxn.metabolites.get(cobra_model.metabolites.get_by_id("PROTON_"+c))
                            rxn.add_metabolites({cobra_model.metabolites.get_by_id("PROTON_"+c):-1*tempCoeff})

                            tempCoeff = rxn.metabolites.get(cobra_model.metabolites.get_by_id("PROTON_"+this))
                            rxn.add_metabolites({cobra_model.metabolites.get_by_id("PROTON_"+this):-1*tempCoeff})

                            rxn.add_metabolites({cobra_model.metabolites.get_by_id("PROTON_"+c):prot_bal})
                        elif prot_bal < 0:
                            tempCoeff = rxn.metabolites.get(cobra_model.metabolites.get_by_id("PROTON_"+c))
                            rxn.add_metabolites({cobra_model.metabolites.get_by_id("PROTON_"+c):-1*tempCoeff})

                            tempCoeff = rxn.metabolites.get(cobra_model.metabolites.get_by_id("PROTON_"+this))
                            rxn.add_metabolites({cobra_model.metabolites.get_by_id("PROTON_"+this):-1*tempCoeff})

                            rxn.add_metabolites({cobra_model.metabolites.get_by_id("PROTON_"+c):prot_bal})
                        else:
                            tempCoeff = rxn.metabolites.get(cobra_model.metabolites.get_by_id("PROTON_"+c))
                            rxn.add_metabolites({cobra_model.metabolites.get_by_id("PROTON_"+c):-1*tempCoeff})

                            tempCoeff = rxn.metabolites.get(cobra_model.metabolites.get_by_id("PROTON_"+this))
                            rxn.add_metabolites({cobra_model.metabolites.get_by_id("PROTON_"+this):-1*tempCoeff})
            else:
                rxn.add_metabolites({cobra_model.metabolites.get_by_id(met[len(prefix):]):coeff1+coeff2,cobra_model.metabolites.get_by_id("PROTON_"+cobra_model.metabolites.get_by_id(met).compartment):round(Charge,4)})
    uncharged = cobra_model.copy()
    return uncharged

############################################################################################
#5.    As a result of the above function the model will be in its default state, this function-
       #helps to update the pH to the compartment. Under different pHs the metabolites will have 
       #different fractional charge forms that also can be update using this function by uploading 
       #an excel containing the different fractional charge forms at different pHs.
############################################################################################

#This function apply the new pH to the compartment
def convertToFractionalCharges(uncharged,infile="MetaboliteChargedStates.xlsx",compH={"v2":5.5},TransferTag=""):
    model = uncharged.copy()

    from xlrd import open_workbook
    from cobra.core import Metabolite

    wb = open_workbook(infile)

    FractionDict = dict()
    FractionCharge = dict()
    for comp in compH.keys():
        pH = compH.get(comp)
        tempDict1 = dict()
        for sheet in wb.sheets():
            met = str(sheet.name)
            metID = met+"_"+comp
            tempDict2 = dict()
            for i in range(1,sheet.ncols):
                charge = sheet.cell(0,i).value
                if charge == "" or charge == " ":
                   continue
                row = int((pH*10)-9)
                pH = sheet.cell(row,0).value
                n = int(model.metabolites.get_by_id(metID).charge) - charge
                prefix = ""
                if n < 0:
                    prefix = "a"*int(abs(n))
                elif n > 0:
                    prefix = "b"*int(abs(n))
                met_spcl = prefix+met
                FractionCharge[met_spcl]=int(charge)
                fra = sheet.cell(row,i).value
                if fra > 0:
                    if tempDict2.keys().__contains__(met_spcl):
                        tempDict2[met_spcl] = tempDict2[met_spcl]+fra
                    else:
                        tempDict2[met_spcl] = fra
            tempDict1[met]=tempDict2
        FractionDict[pH] = tempDict1

    rxnUpdated = set()

    for comp in compH.keys():
        pH = compH.get(comp)
        for met in FractionDict[pH].keys():
            for rxn in model.metabolites.get_by_id(met+"_"+comp).reactions:
                if rxn.id.__contains__("Transfer"):
                    rxnUpdated.add(rxn)
                coeff = rxn.metabolites.get(model.metabolites.get_by_id(met+"_"+comp))
                rxn.add_metabolites({model.metabolites.get_by_id(met+"_"+comp):-1*coeff})
                for met_spcl in FractionDict[pH][met].keys():
                    coeff_new = coeff*FractionDict[pH][met][met_spcl]*0.01
                    f = 0
                    for m in model.metabolites.query(met_spcl+"_"+comp):
                        if m.id==met_spcl+"_"+comp:
                            f=1
                            break
                    if f == 0:
                        m=Metabolite(met_spcl+"_"+comp)
                        m.name = met_spcl+"_"+comp
                        m.formula = model.metabolites.get_by_id(met+"_"+comp).formula
                        if "a" in met_spcl:
                            if m.formula.__contains__("H"):
                                parts = m.formula.split("H")
                                m.formula = parts[0]+"H"+str(int(parts[1][0])+(met_spcl.count("a")))+parts[1][1:]
                            else:
                                m.formula = m.formula+"H"+(len(met_spcl)*1)
                        else:
                            if m.formula.__contains__("H"):
                                parts = m.formula.split("H")
                                m.formula = parts[0]+"H"+str(int(parts[1][0])-(met_spcl.count("a")))+parts[1][1:]
                        m.notes = model.metabolites.get_by_id(met+"_"+comp).notes
                        m.charge = FractionCharge[met_spcl]
                        m.compartment = comp
                    rxn.add_metabolites({m:round(coeff_new,3)})
                Charge = 0
                for tmet in rxn.metabolites.keys():
                    if tmet.charge is not None:
                        Charge = Charge + (int(tmet.charge) * rxn.metabolites.get(tmet))
                    else:
                        print(f"Warning: Metabolite {tmet} has no charge information.")
                if Charge != 0:
                    rxn.add_metabolites({model.metabolites.get_by_id("PROTON_"+model.metabolites.get_by_id(met+"_"+comp).compartment):-1*Charge})
    for rxn in rxnUpdated:
        compSet = set()
        n = 0
        for met in rxn.metabolites.keys():
            if met.id.__contains__("PROTON"):
                n=n+1
            compSet.add(met.compartment)
        c=list(compSet)
        a=c[0][len(c[0])-1]
        b=c[1][len(c[1])-1]
        if a > b:
            thisC = c[1]
            nextC = c[0]
        else:
            thisC = c[0]
            nextC = c[1]
        if rxn.id.__contains__("MAL_v_Transfer"):
           print("This ="+thisC)
           print("Next ="+nextC)
           print(rxn.reaction)
           print(f"thisC: {thisC}, nextC: {nextC}")
        if nextC in ["v1", "v2"]:
            proton_id = "PROTON_" + nextC
        else:
            raise ValueError(f"Unexpected compartment: {nextC}")
        if rxn.metabolites.keys().__contains__(model.metabolites.get_by_id("PROTON_"+thisC)):
            this_prot = rxn.metabolites.get(model.metabolites.get_by_id("PROTON_"+thisC))
            tempCoeff = rxn.metabolites.get(model.metabolites.get_by_id("PROTON_"+thisC))
            rxn.add_metabolites({model.metabolites.get_by_id("PROTON_"+thisC):-1*tempCoeff})
        else:
            this_prot = 0
        if rxn.metabolites.keys().__contains__(model.metabolites.get_by_id("PROTON_"+nextC)):
            next_prot = rxn.metabolites.get(model.metabolites.get_by_id("PROTON_"+nextC))
            tempCoeff = rxn.metabolites.get(model.metabolites.get_by_id("PROTON_"+nextC))
            rxn.add_metabolites({model.metabolites.get_by_id("PROTON_"+nextC):-1*tempCoeff})
        else:
            next_prot = 0
        net = this_prot+next_prot
        if net > 0 and "_v_dielTransfer" not in rxn.id:
            rxn.add_metabolites({model.metabolites.get_by_id("PROTON_"+nextC):net})
        elif net < 0 and "_v_dielTransfer" not in rxn.id:
            rxn.add_metabolites({model.metabolites.get_by_id("PROTON_"+thisC):net})
        elif net > 0 and "_v_dielTransfer" in rxn.id:
            rxn.add_metabolites({model.metabolites.get_by_id("PROTON_"+thisC):net})
        elif net < 0 and "_v_dielTransfer" in rxn.id:
            rxn.add_metabolites({model.metabolites.get_by_id("PROTON_"+nextC):net})
    return model
    
