#%% Import packages
import time
start_time = time.time()

# Scientific and data packages
import numpy as np
import scipy as sp
from PIL import Image as im
from PIL import ImageDraw
import sympy as smp
import sys
import pandas as pd
pd.set_option('display.max_rows',100)
pd.set_option('display.min_rows',100)
pd.set_option('display.max_columns',100)
import seaborn as sns

# Plotting packages
import matplotlib as mpl
import matplotlib.cm as cm
from matplotlib import figure
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
plt.rcParams['figure.figsize'] = [25,25] # Default plot size
plt.rcParams['font.size']=48
plt.rcParams['lines.linewidth'] = 8

# Utility packages
import itertools as iter
from datetime import datetime,date
import os
home_dir = os.getcwd()
import re

# Set numpy print options for debugging
np.set_printoptions(threshold=100000)

# BioPython imports for sequence handling
from Bio.Data.IUPACData import protein_letters_3to1
from Bio.Data.IUPACData import protein_letters_1to3
from Bio import AlignIO,Align
from Bio.SeqIO import parse,to_dict

#%% Sequences
cell_time = time.time()

# List of species/clades to process
clade = np.array([
    'Atrichum undulatum','Marchantia polymorpha','Pteridium aquilinum','Platycerium superbum',
    'Cycas panzhihuaensis','Metasequoia glyptostroboides','Nymphaea alba','Agave victoriae',
    'Carpobrotus edulis','Echeveria elegans','Drosera capensis','Drosera venusta','Sarracenia flava',
    'Arabidopsis thaliana','Ceratophyllum demersum','Chenopodium album','Crithmum maritimum',
    'Dactylis glomerata','Flaveria pringlei','Flaveria floridana','Iris douglasiana',
    'Limonium latebracteatum','Limonium stenophyllum','Limonium virgatum','Nicotiana tabacum',
    'Pallenis maritima','Sideritis cretica','Spinacia oleracea','Teucrium heterophyllum',
    'Trachycarpus fortunei','Triticum aestivum','Flaveria bidentis','Zea mays',
    'Euglena gracilis','Chlamydomonas reinhardtii'
])

os.chdir(f'{home_dir}/input')
dateStamp = f"{datetime.today().strftime('%Y%m%d')}-normed"
spec_names = pd.read_csv('names.txt',names=['id','name'])
dataSet = '2016'
# Read multiple sequence alignment
fullAlignments = AlignIO.read('clustalo-I20210519-163022-0174-96769462-p1m.clustal_num','clustal')
clade_names = pd.DataFrame(clade,columns=['name'])
# Map clade names to IDs using spec_names
clade_names['id'] = [spec_names.id[np.where([nc in st for st in spec_names.name])[0][0]] for nc in clade]
cladeAlignments = Align.MultipleSeqAlignment([])
nonCladeAlignments = Align.MultipleSeqAlignment([])
# Separate alignments into clade and non-clade
for seq in fullAlignments:
    if any([ids in seq.description for ids in clade_names.id]):
        i = int(np.where([ids in seq.description for ids in clade_names.id])[0][0])
        seq.id = clade_names.id[i]
        seq.name = clade_names.name[i]
        cladeAlignments.append(seq)
    else:
        nonCladeAlignments.append(seq)
# Read kinetic rates data
rates = pd.read_excel('Rubisco Kinetics_20160530.xlsx',index_col=0,skiprows=11,usecols='A:J',engine='openpyxl')
rates.columns = np.concatenate([rates.columns[:3],['Kc', 'Ko', 'KcatC','Sco', 'KcatO','Reference']])
# Clean up rates columns
for stat in ['Kc', 'Ko', 'KcatC','Sco', 'KcatO']:
    rates[stat] = [float(x.split()[0]) if type(x)==str else float(x) if type(x)==int else x if type(x)==float else np.nan for x in rates[stat]]
rates = rates.loc[rates.index.intersection(clade_names.name)]

# Main loop over each wild type organism in clade
for wt_org in clade:
    print(wt_org)

    # Create wild type sequence and related variables
    alignments=cladeAlignments
    freqWT=wt_org.split(' ')[0][:3]+wt_org.split(' ')[1][:3]
    wtName = freqWT
    positionCount = []
    # Count amino acid frequency at each position
    for j in range(len(alignments[0])):
        positionCount.append(dict.fromkeys(alignments[:,j],0))
        for i in alignments[:,j]:
            positionCount[j][i] += 1
    # Determine wild type sequence
    if freqWT=='Freq':
        wt = ''
        for dic in positionCount:
            wt += max(dic, key=lambda key: dic[key])
    else:
        wt = str([rec.seq for rec in alignments if rec.name==wt_org][0])
    arabSeq = str([rec.seq for rec in alignments if rec.name==wt_org][0])
    arab_NoDash = [[i,j] for i,j in enumerate(arabSeq) if j!='-']
    wt_NoDash = [[i,j] for i,j in enumerate(wt) if j!='-']
    # DataFrame to store mutation information
    variantMatters = pd.DataFrame(np.full((len(alignments),len(wt)),-1),index = [pl.name for pl in alignments])
    # Map wild type to Arabidopsis sequence positions
    wt2Arab = pd.DataFrame([[[j+1 for resi in arab_NoDash if resi[0]==resj[0]][0],
                             [i+1 for i,resi in enumerate(arab_NoDash) if resi[0]==resj[0]][0],
                             [resi[1] for resi in arab_NoDash if resi[0]==resj[0]][0],
                             [resj[1] for resi in arab_NoDash if resi[0]==resj[0]][0]]
                            for j,resj in enumerate(wt_NoDash)],
                           columns=[f'{wtName}WTLoc','ArabLoc','ArabRes','wtRes'])
    wt2Arab[['nodashLoc','nodashRes']] = [[i+1,j] for i,j in enumerate(wt) if j!='-']
    proteinLength = len([i for i in wt if i!='-'])

    # Identify mutations for each sequence
    for pl in alignments:
        variantMatters.loc[pl.name] = np.array([
            np.nan if (pl.seq[i] == wt[i])|(pl.seq[i] == '-')|(wt[i]=='-')
            else wt[i]+str(np.where(np.array(wt_NoDash)==str(i))[0][0]+1)+pl.seq[i]
            for i in range(len(wt))
        ])
    # Collect unique mutation codes
    mutation_codes = []
    for i in np.arange(len(wt)):
        mutation_codes = mutation_codes + list(variantMatters[i].value_counts().index)
    mutation_codes = np.array(list(dict.fromkeys(mutation_codes)))
    mutation_codes = mutation_codes[mutation_codes!='nan']

    # Create dataframes grouping by mutation
    by_mutation = {}
    by_mutation_nogroup = []
    for mut in mutation_codes:
        orgs = []
        for i in np.arange(len(variantMatters)):
            if variantMatters.iloc[i].isin([mut]).any():
                orgs.append(variantMatters.index[i])
                by_mutation_nogroup.append([mut,variantMatters.index[i],int(mut[1:-1]),np.nan,np.nan,np.nan,np.nan,np.nan])
        by_mutation[mut] = orgs

    by_mutation_nogroupDF = pd.DataFrame(by_mutation_nogroup,columns = ['Mutation','Species','Location','Kc', 'Ko', 'KcatC','Sco', 'KcatO'])
    by_mutationDF = pd.DataFrame([[by_mutation[mut],int(mut[1:-1]),np.nan,np.nan,np.nan,np.nan,np.nan] for mut in by_mutation.keys()],
                                index=by_mutation.keys(),
                                columns = ['Species','Location','Kc', 'Ko', 'KcatC','Sco', 'KcatO'])

    os.chdir(home_dir)
    # Normalize kinetic rates for each variable
    rates['Kc_normed'] = [(x - np.min(rates.Kc))/(np.max(rates.Kc)-np.min(rates.Kc)) for x in rates.Kc]
    rates['Ko_normed'] = [(x - np.min(rates.Ko))/(np.max(rates.Ko)-np.min(rates.Ko)) for x in rates.Ko]
    rates['KcatC_normed'] = [(x - np.min(rates.KcatC))/(np.max(rates.KcatC)-np.min(rates.KcatC)) for x in rates.KcatC]
    rates['Sco_normed'] = [(x - np.min(rates.Sco))/(np.max(rates.Sco)-np.min(rates.Sco)) for x in rates.Sco]
    rates['KcatO_normed'] = [(x - np.min(rates.KcatO))/(np.max(rates.KcatO)-np.min(rates.KcatO)) for x in rates.KcatO]
    rates.index = [st.strip() for st in rates.index]

    # Assign normalized rates to mutation dataframes
    by_mutation_nogroupDF = by_mutation_nogroupDF.loc[by_mutation_nogroupDF.Species.isin(rates.index)]
    by_mutation_nogroupDF.reset_index(drop=True,inplace=True)
    for i,sp in enumerate(by_mutation_nogroupDF.Species):
        by_mutation_nogroupDF.loc[i,['Kc', 'Ko', 'KcatC','Sco', 'KcatO']] = rates.loc[sp,['Kc_normed', 'Ko_normed', 'KcatC_normed','Sco_normed', 'KcatO_normed']].values
    for mut in mutation_codes:
        spes = [li for li in by_mutation[mut] if li in rates.index]
        by_mutationDF.loc[mut,['Kc', 'Ko', 'KcatC','Sco', 'KcatO']] = rates.loc[spes][['Kc', 'Ko', 'KcatC','Sco', 'KcatO']].mean()
    by_mutationDF = by_mutationDF.dropna(axis=0,thresh=2,subset=['Kc','Ko','KcatC','KcatO','Sco'])
    # Calculate wild type position indices
    wtdash = [m.start() for m in re.finditer('-',wt)]
    arabdash = [m.start() for m in re.finditer('-',arabSeq)]
    wtPos = [i for i in range(len(wt)) if i not in wtdash]
    arabPos = [i for i in range(len(arabSeq)) if i not in arabdash]
    by_mutationDF[f'{wtName}WTLoc'] = (by_mutationDF.Location - 1)/(proteinLength - 1)
    by_mutation_nogroupDF[f'{wtName}WTLoc'] = by_mutation_nogroupDF.Location

    # Normalize mutation dataframe for kriging
    by_mutationDF_normed = by_mutationDF.copy()
    for col in by_mutationDF_normed.columns[2:-1]:
        by_mutationDF_normed[col] = (by_mutationDF_normed[col] - np.min(by_mutationDF_normed[col]))/(np.max(by_mutationDF_normed[col]) - np.min(by_mutationDF_normed[col]))
    by_mutationDF_normed[f'{wtName}WTLoc'] = (by_mutationDF_normed[f'{wtName}WTLoc']-1)/(len([i for i in wt if i!='-'])-1)
    by_mutation_nogroupDF[f'{wtName}WTLoc'] = (by_mutation_nogroupDF[f'{wtName}WTLoc']-1)/(len([i for i in wt if i!='-'])-1)

    print(f'---Cell run time: {round(np.floor((time.time() - cell_time)/60))} minutes, {round((time.time() - cell_time)%60)} seconds---')
    print(f'Cell completed: {datetime.now()}')

    # Generate outputs for kriging
    cell_time = time.time()

    def specific_mutation_kriging_output(innie,date,note,X,Y,Z,filterN=0,filter_UQ_LQ_NL_NU=0,norm=False,save=False,scatter=True,scatterSave=False,maha=False,yzscatter=True,yzscatterSave=False):
        """
        Generate kriging input/output files and plots for mutation data.
        innie: input dataframe
        date, note: output file naming
        X, Y, Z: variable names for axes
        filterN, filter_UQ_LQ_NL_NU: filtering options
        norm: normalize data
        save: save output files
        scatter: show scatter plot
        scatterSave: save scatter plot
        yzscatter: show Y vs Z scatter
        yzscatterSave: save Y vs Z scatter
        maha: unused
        """
        indir = os.getcwd()
        export_name = date+'-'+freqWT+'WTLoc-'+Y+'-'+Z+'-'+note
        df = innie.copy()
        # Optional filtering by group size
        if filter_UQ_LQ_NL_NU!=0:
            if filter_UQ_LQ_NL_NU==1:
                df = df[[len(spe)>=np.nanquantile([len(spec) for spec in df.Species],filterN) for spe in df.Species]]
            elif filter_UQ_LQ_NL_NU==2:
                df = df[[len(spe)<=np.nanquantile([len(spec) for spec in df.Species],filterN) for spe in df.Species]]
            elif filter_UQ_LQ_NL_NU==3:
                df = df[[len(spe)<=filterN for spe in df.Species]]
            elif filter_UQ_LQ_NL_NU==4:
                df = df[[len(spe)>=filterN for spe in df.Species]]
        # Select and rename columns for kriging
        df = df[[X,Y,Z]]
        df.columns=['x','y','z']
        df = df.dropna()
        # Y vs Z scatter plot
        if yzscatter:
            fig,ax = plt.subplots()
            ax.scatter(df.z,df.y,s=100)
            ax.set_xlabel(Z,fontweight='bold')
            ax.set_ylabel(Y,fontweight='bold')
            ax.grid(True)
            ax.set_xlim(-.01,1.01)
            ax.set_ylim(-.01,1.01)
            ax.set_aspect('equal')
            fig.tight_layout()
            plt.show()
        # Save output files if requested
        if save or yzscatterSave or scatterSave:
            data_directory = f'{home_dir}/output/kriging_runs/{export_name}'
            if not os.path.exists(data_directory):
                os.mkdir(data_directory)
            os.chdir(data_directory)
        if yzscatterSave:
            fig.savefig(export_name+'-yzscatter.png')
        # Normalize data if requested
        if norm:
            df.x = [(xx-1)/(proteinLength-1) for xx in df.x]
            df.y = [(yy - np.min(df.y))/(np.max(df.y)-np.min(df.y)) for yy in df.y]
            df.z = [(zz - np.min(df.z))/(np.max(df.z)-np.min(df.z)) for zz in df.z]
        if save:
            df.groupby(['x','y']).mean().to_csv(export_name+'.csv')
        # Main scatter plot
        if scatter:
            fig,ax = plt.subplots()
            ax.scatter(df.x,df.y,c=df.z,cmap='plasma',s=100,zorder=5)
            cb = fig.colorbar(plt.cm.ScalarMappable(cmap=cm.get_cmap('plasma'), norm=mpl.colors.Normalize(vmin=0, vmax=1)))
            ax.set_xlabel(X,fontweight='bold')
            ax.set_ylabel(Y,fontweight='bold')
            cb.set_label(Z,fontweight='bold')
            ax.set_xlim(-.01,1.01)
            ax.set_ylim(-.01,1.01)
            ax.set_aspect('equal')
            ax.grid(True)
            plt.show()
            if scatterSave:
                fig.savefig(export_name+'-scatter.png',bbox_inches='tight')
        os.chdir(indir)
        return df

    print(f'---Cell run time: {round(np.floor((time.time() - cell_time)/60))} minutes, {round((time.time() - cell_time)%60)} seconds---')
    print(f'Cell completed: {datetime.now()}')

    # Kriging execution cell

    cell_time = time.time()

    # Loop over all variable pairs for kriging analysis
    vars = ['Kc','Ko','KcatC','KcatO','Sco']
    for yy in vars:
        for zz in vars:
            if yy!=zz:
                runNote = f'{dataSet}-clade_only_alignment'
                specific_mutation_kriging_output(by_mutationDF_normed, dateStamp,runNote,f'{wtName}WTLoc',yy,zz,scatterSave=True,yzscatter=True,yzscatterSave=True,save=False)

    print(f'---Cell run time: {round(np.floor((time.time() - cell_time)/60))} minutes, {round((time.time() - cell_time)%60)} seconds---')
    print(f'Cell completed: {datetime.now()}')

    #%% Data save
    cell_time = time.time()

    # Save all output dataframes and variables for reproducibility
    if not os.path.exists(f'{home_dir}/output'):
        os.mkdir(f'{home_dir}/output')
    os.chdir(f'{home_dir}/output')
    if not os.path.exists(f'{dateStamp}-{wtName}WTLoc-data'):
        os.mkdir(f'{dateStamp}-{wtName}WTLoc-data')
    os.chdir(f'{dateStamp}-{wtName}WTLoc-data')

    by_mutationDF.to_csv(f'xyz_data.csv')
    by_mutationDF.to_pickle(f'{runNote}-by_mutationDF.pkl')
    by_mutationDF_normed.to_pickle(f'{runNote}-by_mutationDF_normed.pkl')
    by_mutation_nogroupDF.to_pickle(f'{runNote}-by_mutation_nogroupDF.pkl')
    rates.to_pickle('rates.pkl')
    np.save('xLength',proteinLength)
    np.save(f'clade_{dataSet}',clade)
    # Save metadata for downstream analysis
    with open('xLength.txt', 'w') as f:
        f.write(str(proteinLength))
    with open('wt_name.txt', 'w') as f:
        f.write(wt_org)
    with open('xName.txt', 'w') as f:
        f.write(f'{wtName}WTLoc')
    with open('y_vars.txt', 'w') as f:
        f.write('Kc Ko KcatC KcatO Sco')
    with open('z_vars.txt', 'w') as f:
        f.write('Kc Ko KcatC KcatO Sco')

    print(f'---Cell run time: {round(np.floor((time.time() - cell_time)/60))} minutes, {round((time.time() - cell_time)%60)} seconds---')
    print(f'Cell completed: {datetime.now()}')

#%% Total run time

print(f'---Total run time: {round(np.floor((time.time() - start_time)/3600))} hours {round(np.floor((time.time() - start_time)/60)%60)} minutes, {round((time.time() - start_time)%60)} seconds---')
print(f'Run completed: {datetime.now()}')