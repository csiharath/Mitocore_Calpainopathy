# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from cobra.flux_analysis import flux_variability_analysis
from cobra.flux_analysis.parsimonious import pfba
import cobra


def get_sum_fluxes(metabolite):
	summarylines = str(metabolite.summary()).split("\n")
	stop = 0
	if summarylines[6] == "Empty DataFrame":
		print("No producing reactions")
		return 0
	i = 0
	for line in summarylines[7:]:
		i += 1
		if line == "Consuming Reactions":
			stop = i-2
	sommefluxes = 0.0
	for line in summarylines[7:7+stop]:
		flux = float(line.split()[1])
		sommefluxes += flux
	return sommefluxes


def get_fluxes_from_mitocore_metabolite(metabolite, model):
	"""
	fonction qui a partir d'un métabolite, récupère tout les flux producteurs du
	métabolites, et fais la somme des flux par sous système
	"""
	summarylines = str(metabolite.summary()).split("\n")
	# print(metabolite.summary())
	stop = 0
	i = 0
	dictflux = dict()
	if summarylines[6] == "Empty DataFrame":
		print("No producing reactions")
		return 0
	for line in summarylines[7:]:
		i += 1
		if line == "Consuming Reactions":
			stop = i-2
	sommefluxes = 0.0
	for line in summarylines[7:7+stop]:
		flux = float(line.split()[0].strip("%"))
		sommefluxes += flux
		r_id = line.split()[2]
		subsystem = model.reactions.get_by_id(r_id).notes['SUBSYSTEM']
		if subsystem not in dictflux.keys():
			dictflux[subsystem] = 0
		dictflux[subsystem] += flux
	# print(sommefluxes)
	return dictflux


def get_mitocore_respiratory_exchange_ratio(model):
	flux_model_o2 = model.reactions.O2t.flux
	flux_model_Co2 = abs(model.reactions.CO2t.flux)
	print("o2 : "+str(flux_model_o2))
	print("Co2 : "+str(flux_model_Co2))
	return flux_model_Co2/flux_model_o2


def define_boundary_and_run_model(model, glucose_upper_bound, FA_upper_bound,
	O2_upper_bound, Acetoacetate_upper_bound=0, hydroxybutyrate_upper_bound=0,
	lactate_upper_bound=0, atp_value=1000,
	glucose_id="GLCt1r", FA_id="HDCAtr", O2_id="O2t", Acetoacetate_id="ACACt2",hydroxybutyrate_id="BHBt", lactate_id="L_LACt2r", OF="OF_ATP_MitoCore",
	glycolyse_target_fonction="PGM", beta_oxydation_target_fonction="r0732", ketone_bodies_target_fonction="BDHm", leucine_target_fonction="LEUTAm", isoleucine_target_fonction="ILETAm",
	acetyl_CoA="accoa_m", FVA=False):
	"""
	fonction qui fait tourner un modèle avec des valeurs d'entrées de métabolites
	principaux choisies. Entrée d'oxygene,Glucose et acide gras requis.
	Le modèle est considéré comme basé sur MitoCore. si ce n'est pas le cas,
	alors il faut modifier les id des réactions associées.
	"""

	# model.reactions.r0732.bounds = (5.257609e-01,1000)

	# Récupération des réactions
	
	glucose_reaction = model.reactions.get_by_id(glucose_id)
	FA_reaction = model.reactions.get_by_id(FA_id)
	O2_reaction = model.reactions.get_by_id(O2_id)
	Acetoacetate_reaction = model.reactions.get_by_id(Acetoacetate_id)
	hydroxybutyrate_reaction = model.reactions.get_by_id(hydroxybutyrate_id)
	lactate_reaction = model.reactions.get_by_id(lactate_id)
	OF_ATP_reaction = model.reactions.get_by_id(OF)

	# définition de la fonction objectif
	model.objective = model.reactions.get_by_id(OF)
	# Changement des upper_bound:
	glucose_reaction.upper_bound = glucose_upper_bound
	FA_reaction.upper_bound = FA_upper_bound
	O2_reaction.upper_bound = O2_upper_bound
	Acetoacetate_reaction.upper_bound = Acetoacetate_upper_bound
	hydroxybutyrate_reaction.upper_bound = hydroxybutyrate_upper_bound
	lactate_reaction.upper_bound = lactate_upper_bound

	if atp_value != 1000:
		OF_ATP_reaction.upper_bound = atp_value
		# OF_ATP_reaction.bounds=(atp_value, atp_value)

	# fva
	if FVA :
		fva = flux_variability_analysis(model)
	# if calpain:
	# 	print(fva.loc['LDH_L']['minimum'])
	# 	model.reactions.LDH_L.upper_bound = fva.loc['LDH_L']['minimum']
	
	# calcul de la solution
	solution = model.optimize()
	# solution = pfba(model)

	# récupération des valeurs de flux de chaque voie métabolique
	dict_voie_metabolique = dict()
	dict_voie_metabolique["glycolysis"]=abs(model.reactions.get_by_id(glycolyse_target_fonction).flux)
	dict_voie_metabolique["beta oxydation"]=abs(model.reactions.get_by_id(beta_oxydation_target_fonction).flux)
	dict_voie_metabolique["ketone bodies"]=abs(model.reactions.get_by_id(ketone_bodies_target_fonction).flux)
	dict_voie_metabolique["leucine degradation"]=abs(model.reactions.get_by_id(leucine_target_fonction).flux)
	dict_voie_metabolique["isoleucine degradation"]=abs(model.reactions.get_by_id(isoleucine_target_fonction).flux)
	# récupération des pourcentagest de l'acétyl-CoA produit par différentes voies métaboliques.
	# dict_accoa = 1
	# dict_atp = 1
	dict_accoa = get_fluxes_from_mitocore_metabolite(model.metabolites.get_by_id(acetyl_CoA),model)
	dict_atp = get_fluxes_from_mitocore_metabolite(model.metabolites.get_by_id("atp_c"),model)


	if FVA:
		return (dict_accoa, dict_voie_metabolique, model, solution, dict_atp, fva)
	else:
		return (dict_accoa, dict_voie_metabolique, model, solution, dict_atp)


def plot_voie_metabolique_ou_accoa(dict_results, nom_graph, accoa=False, prop=False):
	"""
	fonction qui plot les résultats produits par la fonction define_boundary_and_run_model
	il faut lui donner un dictionnaire de dictionnaire résultats
	accoa = True si le dictionnaire donné est un dictionnaire de valeur d'accoa
	"""
	list_voie = list()
	list_isoleucine_degration = []
	list_ketone_bodies = []
	list_ketogenesis_leucine_degradation = []
	list_FA_metabolism = []
	list_Glycolysis = []
	list_intensite = []
	# changement des clé du dictionnaire en fonction de si on analyse l'accoa ou les voies métaboliques
	if not accoa:
		cle_isoleucine_degradation = "isoleucine degradation"
		cle_ketone_bodies = 'ketone bodies'
		cle_leucine_degradation = 'leucine degradation'
		cle_beta_oxydation = 'beta oxydation'
		cle_glycolyse = "glycolysis"
	else:
		cle_isoleucine_degradation = "Isoleucine degradation"
		cle_ketone_bodies = 'FA and ketone body metabolism / Ketogenesis'
		cle_leucine_degradation = 'Ketogenesis / Leucine degradation'
		cle_beta_oxydation = 'FA metabolism'
		cle_glycolyse = "TCA cycle"

	for intensité in dict_results:
		# print(intensité)
		for voie in dict_results[intensité]:
			if voie not in list_voie:
				list_voie.append(voie)
	for intensité in dict_results:
		for voie in list_voie:
			if voie not in dict_results[intensité].keys():
				dict_results[intensité][voie] = 0

	for intensité in dict_results:

		print(dict_results[intensité])
		try:
			list_isoleucine_degration.append(dict_results[intensité][cle_isoleucine_degradation])
		except KeyError:
			list_isoleucine_degration = [0, 0, 0]
			pass

		try:
			list_ketone_bodies.append(dict_results[intensité][cle_ketone_bodies])
		except KeyError:
			list_ketone_bodies = [0, 0, 0]
			pass

		try:
			list_ketogenesis_leucine_degradation.append(dict_results[intensité][cle_leucine_degradation])
		except KeyError:
			list_ketogenesis_leucine_degradation = [0, 0, 0]
			pass

		try:
			list_FA_metabolism.append(dict_results[intensité][cle_beta_oxydation])
		except KeyError:
			list_FA_metabolism = [0, 0, 0]
			pass

		try:
			list_Glycolysis.append(dict_results[intensité][cle_glycolyse])
		except KeyError:
			list_Glycolysis = [0, 0, 0]
			pass

		list_intensite.append(intensité)

	list_isoleucine_degration = np.array(list_isoleucine_degration)
	list_ketone_bodies = np.array(list_ketone_bodies)
	list_ketogenesis_leucine_degradation = np.array(list_ketogenesis_leucine_degradation)
	list_FA_metabolism = np.array(list_FA_metabolism)
	list_Glycolysis = np.array(list_Glycolysis)
	list_intensite = np.array(list_intensite)
	print(list_intensite)

	# width = 0.5
	plt.figure(figsize=(3,5))
	try:
		plt.bar(list_intensite, list_Glycolysis, color='#024547', label='glycolysis') #DCCCA3 
	except ValueError:
		pass
	try:
		plt.bar(list_intensite, list_FA_metabolism, bottom=list_Glycolysis, color='#037E81', label='FA metabolism') #824C71 
	except ValueError:
		pass
	try:
		plt.bar(list_intensite, list_ketone_bodies, bottom=list_FA_metabolism+list_Glycolysis, color='#4CA051', label='ketone bodies') #90AA86 
	except ValueError:
		pass
	try:
		plt.bar(list_intensite, list_isoleucine_degration, bottom=list_FA_metabolism+list_Glycolysis+list_ketone_bodies, color='#95C121', label='isoleucine degradation') #darkslategrey 
	except ValueError:
		pass
	try:
		plt.bar(list_intensite, list_ketogenesis_leucine_degradation, bottom=list_FA_metabolism+list_Glycolysis+list_ketone_bodies+list_isoleucine_degration, color='#D9E6B1', label='leucine degradation') #CEABB1 
	except ValueError:
		pass
	plt.xlabel("Exercise intensity (VO2max percentage)", fontsize=13)
	plt.xticks(fontsize=12)
	plt.yticks(fontsize=12)
	if not accoa and not prop:
		plt.ylabel("Flux", fontsize=13)
	elif prop:
		plt.ylabel("ATP production pathways proportions", fontsize=13)
		plt.legend(loc=(1.04, 0))
		# plt.ylim((0, 1))
	else:
		plt.ylabel("Percentage", fontsize=13)
		plt.legend(loc=(1.04, 0))
		plt.ylim = 1
	
	plt.savefig('Plots/'+nom_graph, bbox_inches='tight')
	plt.show()
	return dict_results

def run_model_calpainopathy(model,glucose_upper_bound,FA_upper_bound,O2_upper_bound,modification_file,Acetoacetate_upper_bound=0,hydroxybutyrate_upper_bound=0,
	lactate_upper_bound=0,
	glucose_id="GLCt1r",FA_id="HDCAtr",O2_id="O2t",Acetoacetate_id="ACACt2",hydroxybutyrate_id="BHBt",lactate_id="L_LACt2r",OF="OF_ATP_MitoCore",
	glycolyse_target_fonction="PGM",beta_oxydation_target_fonction="r0732",ketone_bodies_target_fonction="BDHm",leucine_target_fonction="LEUTAm",isoleucine_target_fonction="ILETAm",
	acetyl_CoA="accoa_m"):
	"""
	fonction qui récupère les valeurs de flux d'un modèle sain
	puis modifie les Upperbounds des réactions dans le fichier model.convertion.csv avec le facteur associé

	WARNING : Il faut avoir fait tourner le modèle sain d'abord (model.optimize())
	"""
	##############################################################################

	##############################################################################
	glucose_reaction = model.reactions.get_by_id(glucose_id)
	FA_reaction = model.reactions.get_by_id(FA_id)
	O2_reaction = model.reactions.get_by_id(O2_id)
	Acetoacetate_reaction = model.reactions.get_by_id(Acetoacetate_id)
	hydroxybutyrate_reaction = model.reactions.get_by_id(hydroxybutyrate_id)
	lactate_reaction = model.reactions.get_by_id(lactate_id)
	# définition de la fonction objectif
	model.objective = model.reactions.get_by_id(OF)
	# Changement des upper_bound:
	glucose_reaction.upper_bound = glucose_upper_bound
	FA_reaction.upper_bound = FA_upper_bound
	O2_reaction.upper_bound = O2_upper_bound
	Acetoacetate_reaction.upper_bound = Acetoacetate_upper_bound
	hydroxybutyrate_reaction.upper_bound = hydroxybutyrate_upper_bound
	lactate_reaction.upper_bound = lactate_upper_bound
	# calcul de la solution
	solution = model.optimize()
	# modification des réactions par rapport à la littérature
	list_modified_reactions = list()
	count = 0
	with open(modification_file, "r") as f:
		f.readline()
		line=f.readline()
		while line!="":
			data=line.split(",")
			gene=data[0].strip('"')
			mitocore_id=data[1].strip('"').strip()
			#print(mitocore_id)
			fold_change=data[2].strip('"')
			multiplicateur=data[3].split('"')[1]
			if multiplicateur == "Hausse":
				multiplicateur=2.0
			if multiplicateur=="Baisse" :
				multiplicateur=0.5
			if multiplicateur=="Warning" :
				multiplicateur=0.5


			if mitocore_id!="None" :

				list_r_id=mitocore_id.split(";")

				for rid in list_r_id :
					list_modified_reactions.append(rid)
				#print(rid)
					r=model.reactions.get_by_id(rid)
					# print(r)
					#print(r.upper_bound)
					#print(r.lower_bound)
					#print(multiplicateur)

					# print(r.flux)
					if r.flux != 0.0 :
						#print(float(r.flux)*float(multiplicateur))
						print(r.flux)
						print(multiplicateur)
						r.upper_bound=max(float(r.flux)*float(multiplicateur),0)
						print()
					#print(r.upper_bound)
					count+=1

			line=f.readline()
	#print(count)
	solution=model.optimize()
	dict_voie_metabolique=dict()
	dict_voie_metabolique["glycolyse"]=abs(model.reactions.get_by_id(glycolyse_target_fonction).flux)
	dict_voie_metabolique["beta oxydation"]=abs(model.reactions.get_by_id(beta_oxydation_target_fonction).flux)
	dict_voie_metabolique["ketone bodies"]=abs(model.reactions.get_by_id(ketone_bodies_target_fonction).flux)
	dict_voie_metabolique["leucine degradation"]=abs(model.reactions.get_by_id(leucine_target_fonction).flux)
	dict_voie_metabolique["isoleucine degradation"]=abs(model.reactions.get_by_id(isoleucine_target_fonction).flux)
	#récupération des pourcentagest de l'acétyl-CoA produit par différentes voies métaboliques.
	dict_accoa=get_fluxes_from_mitocore_metabolite(model.metabolites.get_by_id(acetyl_CoA),model)
	#
	fva = flux_variability_analysis(model)
	
	return (dict_accoa,dict_voie_metabolique,model,solution,fva)

# def get_boundaries_and_run_model_calpain(model,glucose_upper_bound,FA_upper_bound,O2_upper_bound,Acetoacetate_upper_bound=0,hydroxybutyrate_upper_bound=0,
# 	lactate_upper_bound=0, atp_value=1000,
# 	glucose_id="GLCt1r",FA_id="HDCAtr",O2_id="O2t",Acetoacetate_id="ACACt2",hydroxybutyrate_id="BHBt",lactate_id="L_LACt2r",OF="OF_ATP_MitoCore",
# 	glycolyse_target_fonction="PGM",beta_oxydation_target_fonction="r0732",ketone_bodies_target_fonction="BDHm",leucine_target_fonction="LEUTAm",isoleucine_target_fonction="ILETAm",
# 	acetyl_CoA="accoa_m"):
# 	"""

# 	"""
# 	return (dict_accoa,dict_voie_metabolique,model,solution, dict_atp, fva)

def heatmap_flux(minimal_flux=0.05,save=False,**kwargs):
	"""
	Fonction qui plot un heat map en utilisant Seaborn
	la fonction prend n'importe quel nombre d'objet solutions produit par run_model_calpainopathie() ou define_boundary_and_run_model()
	"""
	list_df=list()
	for solution in kwargs["solutions"] :
		series=solution.fluxes[solution.fluxes!=0] #remove 0 values
		df=series.to_frame()
		df["name"]=df.index
		df = df.drop(df[df.fluxes < minimal_flux].index) #enlève les flux avec une valeur inférieur a minimal_flux
		list_df.append(df)
	#Quadraboucle qui s'assure que chaque DF contient les mêmes labels (que toutes les réactions soient présente dans chaque DF)
	#si ce n'est pas le cas rajoute une ligne avec pour valeur 0 dans les dataframe qui n'ont pas le label (Nom de réaction)
	for df in list_df :
		for i in df.index :
			for df2 in list_df :
				if i not in df2.index :
					df2.loc[i]=[0,i]
	# combine tous les DF dans un seul DF nommé results
	results=pd.DataFrame()
	compteur=0
	for name in kwargs["name_columns"] :
		results[name]=list_df[compteur]["fluxes"]
		compteur+=1

	#plot
	fig, ax = plt.subplots(figsize=(20,50))
	ax=sns.heatmap(results)
	#savefig
	if save==True :
		fig = ax.get_figure()
		fig.savefig(kwargs["name_plot"])
	return results