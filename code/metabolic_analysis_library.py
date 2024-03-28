"""
Library of functions used to build and analyse healthy and deseased models specifically during physical exercise.
This library contains functions to run modelsduring physical exercise, build models a dysregulated model, build data visualisation, and plot results.
"""

# import cobra
# from cobra.core import Metabolite, Reaction, Model
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd
import importlib
# import matplotlib.lines as mlines
import matplotlib
from cobra.flux_analysis import flux_variability_analysis
from cobra.flux_analysis.parsimonious import pfba
# import plotly.io as pio
import plotly.express as px
import d3flux
from jinja2 import Template
# from escher import Builder

############################ RUN MODEL - EXERCISE SIMULATION ################################

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

def define_boundary_and_run_model(model, glucose_upper_bound, FA_upper_bound,
	O2_upper_bound, Acetoacetate_upper_bound=0, hydroxybutyrate_upper_bound=0,
	lactate_upper_bound=0, atp_value=1000,
	glucose_id="GLCt1r", FA_id="HDCAtr", O2_id="O2t", Acetoacetate_id="ACACt2",hydroxybutyrate_id="BHBt", lactate_id="L_LACt2r", OF="OF_ATP_MitoCore",
	glycolyse_target_fonction="PGM", beta_oxydation_target_fonction="r0732", ketone_bodies_target_fonction="BDHm", leucine_target_fonction="LEUTAm", isoleucine_target_fonction="ILETAm",
	acetyl_CoA="accoa_m", FVA=False, pFBA=False):
	"""
    Runs FBA or pFBA with a given model and a set of boundary conditions.
    The model is considered as based on MitoCore. If it is not the case, the id of the reactions must be changed.
    FVA cans also be run if argument "FVA" is set to True.
	"""
	# Gets reactions
	
	glucose_reaction = model.reactions.get_by_id(glucose_id)
	FA_reaction = model.reactions.get_by_id(FA_id)
	O2_reaction = model.reactions.get_by_id(O2_id)
	Acetoacetate_reaction = model.reactions.get_by_id(Acetoacetate_id)
	hydroxybutyrate_reaction = model.reactions.get_by_id(hydroxybutyrate_id)
	lactate_reaction = model.reactions.get_by_id(lactate_id)
	OF_ATP_reaction = model.reactions.get_by_id(OF)

	# Objective function
	model.objective = model.reactions.get_by_id(OF)
	# Upper bounds changes
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
	
	# Runs FBA or pFBA depending on the argument pFBA
	if pFBA:
		solution = pfba(model)
	else:
		solution = model.optimize()
	
	# Gets flux values for each targeted pathway
	dict_voie_metabolique = dict()
	dict_voie_metabolique["glycolysis"]=abs(model.reactions.get_by_id(glycolyse_target_fonction).flux)
	dict_voie_metabolique["beta oxydation"]=abs(model.reactions.get_by_id(beta_oxydation_target_fonction).flux)
	dict_voie_metabolique["ketone bodies"]=abs(model.reactions.get_by_id(ketone_bodies_target_fonction).flux)
	dict_voie_metabolique["leucine degradation"]=abs(model.reactions.get_by_id(leucine_target_fonction).flux)
	dict_voie_metabolique["isoleucine degradation"]=abs(model.reactions.get_by_id(isoleucine_target_fonction).flux)

	dict_accoa = get_fluxes_from_mitocore_metabolite(model.metabolites.get_by_id(acetyl_CoA),model)

	if FVA:
		return (dict_accoa, dict_voie_metabolique, model, solution, fva)
	else:
		return (dict_accoa, dict_voie_metabolique, model, solution)


################################ VISUALISATION FUNCTIONS ####################################


def build_dataframe(model, soluce, threshold_down, threshold_up, 
                    hide_transport=True, from_list = False, 
                    pathways_list = ['TCA cycle']):
    """
    Builds a dataframe of reactions with flux values and subsystems from a FBA 
    solution (soluce) using a range of flux values to keep and optionnally a list of subsystems to keep
    This is used to build a treemap with plotly
    """
    df_percent = soluce.to_frame()
    subsystems = []

    # Gets subsystems for each reactions
    for index, row in df_percent.iterrows():
        notes = model.reactions.get_by_id(index).notes
        if notes != {}:
            subsystems.append(notes['SUBSYSTEM'])
        else:
            subsystems.append('unknown')

    # Adds columns subsystems and reactions to the dataframe
    df_percent['Subsystems'] = subsystems
    df_percent['Reactions'] = df_percent.index

    # Sorts dataframe from given parameters and removes unknown subsystems
    df_percent_sort = df_percent.loc[abs(df_percent['fluxes']) > threshold_down]
    df_percent_sort = df_percent_sort.loc[abs(df_percent_sort['fluxes']) < threshold_up]
    # if not from_list:
    df_percent_sort = df_percent_sort.loc[df_percent_sort['Subsystems'] != 'unknown']

    for index, row in df_percent_sort.iterrows():
        if row['fluxes'] < 0:
            df_percent_sort.loc[index, 'fluxes'] = -row['fluxes']
            df_percent_sort['sign'] = 'negative'
        else:
            df_percent_sort['sign'] = 'positive'

    # Hides some subsytems for better visualisation
    if hide_transport:
        row2drop = []
        subsystems_to_hide = ['Electron transport chain', 'Boundary conditions - core', 'Mitochondrial transport - diffusion / artificial', 'Mitochondrial transporters - characterised', 'Objective Function - ATP']
        # subsystems_to_hide = ['Electron transport chain', 'Boundary conditions - core', 'Mitochondrial transport - diffusion / artificial']
        
        for index, row in df_percent_sort.iterrows():
            if any(x in row['Subsystems'] for x in subsystems_to_hide):
                row2drop.append(index)

        df_percent_sort = df_percent_sort.drop(row2drop)
    
    # Sorts dataframe from a list of subsystem(s) to keep
    if from_list:
        row2drop = []
        for index, row in df_percent_sort.iterrows():
            if row['Subsystems'] not in pathways_list:
                row2drop.append(index)

        df_percent_sort = df_percent_sort.drop(row2drop)       

    return df_percent_sort

def display_treemap(dataframe, name):
  """
  Displays a treemap with plotly from a dataframe built with build_dataframe()
  """

  fig = px.treemap(dataframe,
    path = [px.Constant('Cell'),'Subsystems', 'Reactions'],
    values = 'fluxes',
    color = 'fluxes',
    color_continuous_scale=px.colors.diverging.RdBu[::-1],
    # color_continuous_scale=px.colors.diverging.curl,
    # color_continuous_scale=px.colors.cyclical.HSV,
    # width = 1400,
    # height = 900
    width = 1000,
    height = 900
  )
  fig.update_traces(textinfo = "label+value", textposition="middle center", hovertemplate='Pathway or enzyme=%{label}<br>Flux=%{value}<br>Percentage=.%{percent entry}<extra></extra>')

  # Changement du layout
  fig.update_layout(
    {'font_color': "white",
    'plot_bgcolor': 'rgba(0, 0, 0, 100)',
    'paper_bgcolor': '#2e2d2d'},
    title = f"Treemaps of metabolic pathways for {name}% physical intensity",
    title_x=0.5,
    title_y=0.97
    )

  fig.show()

def build_dataframe_visual(model, soluce, fva = False):
    """
    Builds a dataframe of reactions with flux values and subsystems from a FBA or FVA
    """

    if fva:
        df_percent = soluce
    else:
        df_percent = soluce.to_frame()

    subsystems = []

    for index, row in df_percent.iterrows():
        notes = model.reactions.get_by_id(index).notes
        if notes != {}:
            subsystems.append(notes['SUBSYSTEM'])
        else:
            subsystems.append('unknown')

    df_percent['Subsystems'] = subsystems
    df_percent['Reactions'] = df_percent.index

    if fva:
        df_percent = df_percent.sort_values(by=['Subsystems'])
    else:
        df_percent = df_percent.sort_values(by=['Subsystems', 'fluxes'])

    # print(df_percent.to_markdown())
    return df_percent


def print_reactions_from_met(model, metabolite, solution, notNull = True):
    """
    Prints all reactions and the corresponding flux involving a given metabolite
    """
    reactions = model.metabolites.get_by_id(metabolite).reactions
    for r in reactions:
        flux = solution.to_frame().loc[r.id]
        if notNull:
            if flux['fluxes'] != 0:
                print(r)
                print(solution.to_frame().loc[r.id])
                print('--------')
        else:
            print(r)
            print(solution.to_frame().loc[r.id])
            print('--------')


def get_alternative_pathways(fva_dt):
    """
    Returns all reactions from alternatives pathways (no essantial reactions)
    """
    essential_react = []
    alternative = []

    for index, row in fva_dt.iterrows():
        if row['minimum'] != row['maximum']:
            alternative.append(index)
        else: # if row['minimum'] > 0:
            essential_react.append(index)

    return alternative

def get_reactions_to_check(all_fluxes_dt, alternative, subsys):
    """
    Gets all alternative reactions to check given a subsystem
    """
    reaction_to_check = []

    for react_id in alternative:
        reaction = all_fluxes_dt.loc[react_id]
        if reaction['fluxes'] == 0 and reaction['Subsystems'] == subsys:
            # print(reaction['Subsystems'])
            reaction_to_check.append(react_id)

    return reaction_to_check

def print_alternative_reactions(list_id, fva_dt):
    for react_id in list_id:
        print(fva_dt.loc[react_id])
        print("---------------------------")


def get_key_fluxes(solution, list_enzyme = ["FBA", "PDHm", "LDH_L", "CSm", "r0722", "r0732", "OF_ATP_MitoCore", "O2t"]):
    list_key = []
    df_sol = solution.to_frame()

    for enzyme in list_enzyme:
        list_key.append(df_sol.loc[enzyme]['fluxes'])
    
    return list_key

def create_df_key_fluxes(dict_sol, list_enzyme = ["FBA", "PDHm", "LDH_L", "CSm", "r0722", "r0732", "OF_ATP_MitoCore", "O2t"]):
    """
    Builds a dataframe for the visualisation of key fluxes from a dictionary of solutions
    """
    df_key = pd.DataFrame(index=list_enzyme)
    for sol in dict_sol:
        list_key = get_key_fluxes(dict_sol[sol], list_enzyme)
        df_key[sol] = list_key
    
    return df_key

def show_network(model_json, solution):
    """
    Displays a network from a json file and a solution
    """
    css_custom = \
    """
    .node {
    stroke: #dedddd;
    stroke-width: 1px;
    r: 4px;
    }
    text.nodelabel {
    font-size:5pt;
    fill: #555;
    font-family: Arial;
    }

    text.nodelabel.cofactor {
    font-size: 4pt;
    fill: #778899;
    }

    text.nodelabel.rxn {
    fill: #A9A9A9;
    font-size: 7pt;
    font-style: italic;
    }
    """
    css = Template(css_custom)

    for reaction in model_json.reactions:
        if not reaction.notes['map_info']['hidden']:
            # print(solution100_calpain.fluxes[reaction.id] )
            reaction.notes['map_info']['flux'] = solution.fluxes[reaction.id] 
            if reaction.notes['map_info']['flux'] != 0:
                reaction.notes['map_info']['color'] = 'red'
            else:
                reaction.notes['map_info']['color'] = 'grey'

    map = d3flux.flux_map(model_json, figsize=(1900, 1600), flowLayout=True, custom_css=css_custom)

    return map

################################# MODEL BUILDING FUNCTIONS ######################################

def get_new_bounds(model, reaction, reg, ub85, ub100, upperbounds_before, upperbounds_after, lowerbounds_after, fva_min, fva_max, index, regulation):
    """
    Gets new bound to build the calpainopathy model
    """
    # default values for non numeric regulations
    dict_reg = {'Warning': 0.5, 'Baisse': 0.5, 'Hausse':2}
    index.append(reaction)
    upperbounds_before.append(model.reactions.get_by_id(reaction).upper_bound)
    if ub100.loc[reaction]['Subsystems'] == 'FA metabolism':
        if reg in dict_reg:
            factor = dict_reg[reg]
        else:
            factor = float(reg)
        # model.reactions.get_by_id(reaction).upper_bound = max((ub85.loc[reaction]['maximum']*factor), 0)
        bound = ub85.loc[reaction]['maximum']*factor
        if bound <= 0:
            bound = ub85.loc[reaction]['minimum']*factor
        fva_min.append(ub85.loc[reaction]['minimum'])
        fva_max.append(ub85.loc[reaction]['maximum'])
    else:
        if reg in dict_reg:
            factor = dict_reg[reg]
        else:
            factor = float(reg)
        # model.reactions.get_by_id(reaction).upper_bound = max((ub100.loc[reaction]['maximum']*factor), 0)
        bound = ub100.loc[reaction]['maximum']*factor
        if bound <= 0:
            bound = ub100.loc[reaction]['minimum']*factor
        fva_min.append(ub100.loc[reaction]['minimum'])
        fva_max.append(ub100.loc[reaction]['maximum'])
    regulation.append(factor)
    if bound <= 0:
        # print(reaction)
        # print(type(bound))
        model.reactions.get_by_id(reaction).lower_bound = bound
        lowerbounds_after.append(model.reactions.get_by_id(reaction).lower_bound)
        upperbounds_after.append(model.reactions.get_by_id(reaction).upper_bound)
    else:
        # print(reaction)
        # print(type(bound))
        model.reactions.get_by_id(reaction).upper_bound = bound
        lowerbounds_after.append(model.reactions.get_by_id(reaction).lower_bound)
        upperbounds_after.append(model.reactions.get_by_id(reaction).upper_bound)
    
    return model, upperbounds_before, upperbounds_after, lowerbounds_after, fva_min, fva_max, index, regulation


def create_model_calpain(model, ub100, ub85, csv_file):
    """
    Return a model with new bounds for reactions to regulate and simulate the calpainopathy or an other disease.
    """
    upperbounds_before = []
    upperbounds_after = []
    lowerbounds_after = []
    fva_min = []
    fva_max = []
    regulation = []
    index = []
    df_upperbounds = pd.DataFrame()

    csv = pd.read_csv(csv_file, sep=',')
    reaction_to_regulate = csv['id_reaction'].to_list()

    for reaction in reaction_to_regulate:
        reg = csv.loc[csv['id_reaction'] == reaction].iloc[0]['Regulation']
        if ';' in reaction:
            reaction_list = reaction.split(';')
            for r in reaction_list:
                model, upperbounds_before, upperbounds_after, lowerbounds_after, fva_min, fva_max, index, regulation = get_new_bounds(model, r, reg, ub85, ub100, upperbounds_before, upperbounds_after, lowerbounds_after, fva_min, fva_max, index, regulation)
                # print(model.reactions.get_by_id(r).bounds)
        else:
                model, upperbounds_before, upperbounds_after, lowerbounds_after, fva_min, fva_max, index, regulation = get_new_bounds(model, reaction, reg, ub85, ub100, upperbounds_before, upperbounds_after, lowerbounds_after, fva_min, fva_max, index, regulation)
                # print(model.reactions.get_by_id(reaction).bounds)

    df_upperbounds['reactions'] = index
    df_upperbounds['upperbounds_wt'] = upperbounds_before
    df_upperbounds['fva_min'] = fva_min
    df_upperbounds['fva_max'] = fva_max
    df_upperbounds['regulation'] = regulation
    df_upperbounds['lowerbounds_calpain'] = lowerbounds_after
    df_upperbounds['upperbounds_calpain'] = upperbounds_after
    
    # cobra.io.write_sbml_model(model, filename="Models/Mitocore_mouse_calpain.xml")

    return model, df_upperbounds



########################################### PLOTS ########################################### 

# dict_tot = {'25':0,'65':0,'85':0, '100':0} 
# dict_prop={'25':{},'65':{},'85':{}, '100':{}}
def build_df_proportion(df_pathways):
    """
    Builds a dataframe to plot proportion of pathways. It can be used to compare different conditions 
    for the same intensity.
    """
    
    df_prop = df_pathways.copy()

    for column in df_pathways.columns:
        df_prop[column] = df_pathways[column]/df_pathways[column].sum()

         
    # for intensity in dict_pathway:
    #     # print(intensity)
    #     for pathway in dict_pathway[intensity]:
    #         # print(pathway)
    #         dict_tot[intensity] += dict_pathway[intensity][pathway]
    #         # print(dict_tot[intensity])
    
    # for intensity in dict_pathway:
    #     for pathway in dict_pathway[intensity]:
    #         dict_prop[intensity][pathway] = dict_pathway[intensity][pathway] / dict_tot[intensity]
    
    return df_prop

def barplot_pathways(df_results, nom_graph, prop=False, keys = ['glycolysis', 'beta oxidation', 'Malate-Aspartate shuttle', 'Amino acid degradation'], figsize=(4,5), add_legend=False): #['glycolysis', 'beta oxydation', 'ketone bodies', 'leucine degradation', 'isoleucine degradation']
    # ancien plot_voie_metabolique_ou_accoa()
    """
    Display a barplot of the results from the function define_boundary_and_run_model
	"""
    # matplotlib.rcParams.update(mpl.rcParamsDefault)
    # list_voie = list()
    # list_intensite = []
    # df_results = pd.DataFrame(dict_results)

    # for intensité in dict_results:
	# 	# print(intensité)
    #     for voie in dict_results[intensité]:
    #         if voie not in list_voie:
    #             list_voie.append(voie)
    # for intensité in dict_results:
    #     for voie in list_voie:
    #         if voie not in dict_results[intensité].keys():
    #             dict_results[intensité][voie] = 0

    # lists_fluxes = {key: [] for key in keys}

    # for intensité in dict_results:
    #     for key in keys:
    #         try:
    #             lists_fluxes[key].append(dict_results[intensité][key])
    #         except KeyError:
    #             lists_fluxes[key] = [0, 0, 0, 0]
    #             pass
          
        # list_intensite.append(intensité)
          
    # lists_fluxes = {key: np.array(lists_fluxes[key]) for key in keys}
    
    plt.figure(figsize=figsize)
    plt.style.use("seaborn")
    # cmap = matplotlib.cm.get_cmap('Set3')
    labels = []
    bottoms = np.zeros(len(df_results.columns))
    for row in df_results.iterrows():
        try:
            plt.bar(df_results.columns.to_list(), row[1].to_list(), bottom=bottoms, label=row[0])
            bottoms += row[1]
            labels.append(row[0])
        except ValueError:
            pass

    # for i, key in enumerate(keys):
    #     try:
    #         # print(bottoms)
    #         plt.bar(list_intensite, lists_fluxes[key], bottom=bottoms, label=key) #DCCCA3 , color='#024547'
    #         bottoms += lists_fluxes[key]
    #         labels.append(key)
    #         # print(key)
    #     except ValueError:
    #         pass

    plt.xlabel("Exercise intensity (VO2max percentage)", fontsize=10)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)

    if add_legend:
        plt.legend(reversed(plt.legend().legendHandles), reversed(labels),bbox_to_anchor=(1.7, 0.5), loc='center right', frameon=False)
    
    if not prop:
        print(labels)
        # plt.legend(loc='upper left')
        plt.ylabel("Flux", fontsize=13)
    elif prop:
        plt.ylabel("ATP production pathways proportions", fontsize=13)
		# # plt.legend(reversed(plt.legend().legendHandles), ['leucine degradation', 'isoleucine degradation', 'ketone bodies', 'FA metabolism', 'glycolysis'])
		# plt.legend(loc=(1.04, 0))
		# plt.ylim((0, 1)) 
    else:
        plt.ylabel("Percentage", fontsize=13)
		# plt.legend(reversed(plt.legend().legendHandles), ['leucine degradation', 'isoleucine degradation', 'ketone bodies', 'FA metabolism', 'glycolysis'])
		# plt.legend(loc=(1.04, 0))
        plt.ylim = 1
	
    plt.savefig('Plots/'+nom_graph, bbox_inches='tight')
 
    plt.show()
    #return the final barplot object
    return plt



def create_df_flux_plot(dict_sol_fva, list_enzyme = ["FBA", "PDHm", "LDH_L", "CSm", "r0722", "r0732", "OF_ATP_MitoCore", "O2t"]):
    """
    Builds a dataframe of fluxes and FVA results for a list of key enzymes to be used by fva_plot()
    """
    df = pd.DataFrame(columns=['enzymes', 'intensity', 'fluxes', 'fva_min', 'fva_max'])
    # enzymes = []
    # intensities = []
    # fluxes = []
    # fva_min = []
    # fva_max = []
    row = 0

    for enzyme in list_enzyme:
        for intensity in dict_sol_fva:
            if '25' in intensity:
                percent = '25%'
            elif '65' in intensity:
                percent = '65%'
            elif '85' in intensity:
                percent = '85%'
            else:
                percent = '100%'
            df.loc[row] = [enzyme, percent, dict_sol_fva[intensity][0].to_frame().loc[enzyme]['fluxes'], dict_sol_fva[intensity][1].loc[enzyme]['minimum'], dict_sol_fva[intensity][1].loc[enzyme]['maximum']]
            row += 1
            # enzymes.append(enzyme)
            # intensities.append(intensity)
            # fluxes.append(dict_sol_fva[intensity][0].to_frame().loc[enzyme]['fluxes'])
            # fva_min.append(dict_sol_fva[intensity][1].loc[enzyme]['minimum'])
            # fva_max.append(dict_sol_fva[intensity][1].loc[enzyme]['maximum'])
    
    # df['enzymes'] = enzymes
    # df['intensities'] = intensities
    # df['fluxes'] = fluxes
    # df['fva_min'] = fva_min
    # df['fva_max'] = fva_max

    return df


def fva_plot (dataframe_calpain, dataframe_wt, dict_enzyme): 
    """ 
    Displays a plot of fluxes and the possible range of fluxes for each key enzyme defined in the model
    Compares the wild type and the calpainopathy model
    """
    cmap = matplotlib.cm.get_cmap('tab10')
    count = 0
    plt.style.use('seaborn')
    list_enzyme = list(dict_enzyme.keys())
    legend = dict_enzyme
    fig, axes = plt.subplots(len(dict_enzyme), 2, figsize=(12, 12), sharey='row')

    for enzyme in list_enzyme:
        df_enz = dataframe_calpain.loc[dataframe_calpain['enzymes'] == enzyme]
        df_enz_wt = dataframe_wt.loc[dataframe_wt['enzymes'] == enzyme]

        # y = 0.85-(0.1*count)
        # print(y)
        # fig.text(0.1, y, enzyme, fontsize=12)

        # fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 1.5), sharey=True)
        axes[count][0].scatter(df_enz_wt['intensity'], df_enz_wt['fluxes'], color=cmap(7))
        axes[count][0].plot(df_enz_wt['intensity'], df_enz_wt['fluxes'], color=cmap(7))
        axes[count][0].fill_between(df_enz_wt['intensity'], df_enz_wt['fva_min'], df_enz_wt['fva_max'], alpha=0.2, color=cmap(7), linewidth=0)
        axes[count][0].set_ylabel(legend[enzyme], fontsize=13)
        axes[count][0].tick_params(labelsize=13)

        axes[count][1].scatter(df_enz['intensity'], df_enz['fluxes'], color=cmap(3))
        axes[count][1].plot(df_enz['intensity'], df_enz['fluxes'], color=cmap(3))
        axes[count][1].fill_between(df_enz['intensity'], df_enz['fva_min'], df_enz['fva_max'], alpha=0.2, color=cmap(3), linewidth=0)
        axes[count][1].tick_params(labelsize=13)

        # leg.append(mlines.Line2D([], [], color=cmap(count), marker=markers[count], markersize=10, label=enzyme, ls='-'))
        # fig.legend(handles=[leg])
        count += 1

    # fig.legend(handles=leg, loc="center right")
    # fig.tight_layout()
    fig.align_ylabels()
    fig.text(0.22, 0.89, "CAPN3+ (control)", fontsize=15)
    fig.text(0.68, 0.89, "CAPN3-KO", fontsize=15)
    fig.text(0.5, 0.06, "Exercise intensity per conditions CAPN3+/CAPN3-KO (VO2max percentage)", ha='center', fontsize=15)
    fig.text(0.03, 0.5, 'Flux (mmol.gDW-1.h-1)', va='center', rotation='vertical', fontsize=15)
    # add three vertical bars to specify the pathways
    boundaries = [0, 1, 3, 6]  # The first color spans 2 units, the others 1 unit each
    colors = ['#937860', '#60ab71', '#587ab2']
    cmap = cmap= matplotlib.colors.ListedColormap(colors)
    norm = matplotlib.colors.BoundaryNorm(boundaries, cmap.N, clip=True)
    cbar_ax = fig.add_axes([0.92, 0.09, 0.005, 0.8])
    plt.colorbar(matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cbar_ax, boundaries=boundaries, ticks=[0.5, 2, 4.5], spacing='proportional')
    cbar_ax.set_yticklabels(['TCA cycle', '\u03B2-oxidation', 'Glycolysis'], fontsize=15)
    fig.savefig("Plots/FVA_wt_calpain.png")

def get_dict_pathways(solutions, pathways = {'glycolysis': 'PGM', "beta oxidation":'r0732', "Malate-Aspartate shuttle":'AKGMALtm', "Amino acid degradation":'GLUDxm'}):
    """
    Builds a dictionary of pathways with the corresponding fluxes. The final dictionnary is used to build a barplot
    """
    dict_pathway = dict()
    df_pathway = pd.DataFrame(columns=['25%', '65%', '85%', '100%'])

    for solution in solutions:
        for pathway, reaction in pathways.items():
            try:
                dict_pathway[pathway].append(abs(solution.fluxes[reaction]))
            except KeyError:
                dict_pathway[pathway] = [abs(solution.fluxes[reaction])]

    for pathway in dict_pathway:
        df_pathway.loc[pathway] = dict_pathway[pathway]

    # for pathway, reaction in pathways.items():
    #     dict_pathway[pathway] = abs(solutions.fluxes[reaction])
    # dict_pathway[list(pathways.keys())[0]] = abs(model.reactions.get_by_id(pathways['glycolysis']).flux)
    # dict_pathway[list(pathways.keys())[0]] = abs(model.reactions.get_by_id(pathways['beta oxidation']).flux)
    # dict_pathway[list(pathways.keys())[0]] = abs(model.reactions.get_by_id(pathways['Malate-Aspartate shuttle']).flux)
    # dict_pathway[list(pathways.keys())[0]] = abs(model.reactions.get_by_id(pathways['Amino acid degradation']).flux)

    return df_pathway

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

################################# ADDITIONAL FUNCTIONS ######################################

def get_reactions_from_subsys(model, dt_fluxes, list_of_subsys):
    
    for subsys in list_of_subsys:
        print(f"\n------ Reactions from {subsys} -------\n")
        df_subsys = dt_fluxes.loc[dt_fluxes['Subsystems'] == subsys]
        for index, row in df_subsys.iterrows():
            print(model.reactions.get_by_id(index))
            print(f"{row['fluxes']}\n")

def get_mitocore_respiratory_exchange_ratio(model):
	flux_model_o2 = model.reactions.O2t.flux
	flux_model_Co2 = abs(model.reactions.CO2t.flux)
	print("o2 : "+str(flux_model_o2))
	print("Co2 : "+str(flux_model_Co2))
	return flux_model_Co2/flux_model_o2

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