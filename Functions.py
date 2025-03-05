
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
