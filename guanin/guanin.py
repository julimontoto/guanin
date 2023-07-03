import os
import tempfile
import math
import numpy as np
import pandas as pd
import statistics
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.stats.mstats import gmean
from scipy import stats
# Monkeypatch matplotlib to avoid the failing from upsetplot
from matplotlib import tight_layout
tight_layout.get_renderer = ""
import logging
import argparse
from fpdf import FPDF
from sklearn.preprocessing import StandardScaler, quantile_transform
from sklearn.neighbors import KNeighborsClassifier
from mlxtend.feature_selection import SequentialFeatureSelector as SFS
import seaborn as sns
import time
import pathlib
import webbrowser
from ERgene import FindERG


def getfolderpath(folder):
    '''RCC path'''
    cwd = pathlib.Path(__file__).parent.absolute()
    path = cwd / folder
    return path

def loadrccs(args, start_time = 0):
    """ RCC loading to extract information"""
    columns = ['ID', 'Comments', 'FOV value', 'Binding Density', 'Background', 'Background2', 'Background3', 'Genes below backg %', 'nGenes', 'posGEOMEAN', 'Sum', 'Median', 'R2', 'limit of detection', '0,5fm']
    infolanes = pd.DataFrame(columns = columns)

    dfgenes = pd.DataFrame() #counts
    dfposneg = {} #posnegs

    geomeans = []
    logconc = [7, 5, 3, 1, -1]

    dfnegcount = pd.DataFrame()
    negnames = []
    dfhkecount = pd.DataFrame()
    hkenames = []

    a = 0 #loop count
    for file in os.listdir(getfolderpath(args.folder)):
        '''First data inspection'''
        if '.RCC' in file:
            df = pd.read_csv((getfolderpath(args.folder) / file), names=['CodeClass', 'Name', 'Accession', 'Count']) #count info dataframe
            df = df.dropna()
            df = df.replace('0','1')

            #separate dataframe for gene class
            dfpos = df[df.CodeClass =='Positive']
            dfpos1 = df[df.CodeClass =='Positive1']
            dfpos2 = df[df.CodeClass =='Positive2']
            dfpos = pd.concat([dfpos, dfpos1, dfpos2])
            dfneg = df[df.CodeClass =='Negative']
            dfposneg1 = pd.concat([dfpos,dfneg])
            dfend = df[df.CodeClass == 'Endogenous']
            dfend1 = df[df.CodeClass == 'Endogenous1']
            dfend2 = df[df.CodeClass == 'Endogenous2']
            dfend = pd.concat([dfend, dfend1,dfend2])
            dfhke = df[df.CodeClass =='Housekeeping']

            convert_dict = {'CodeClass': str, 'Name': str, 'Accession': str, 'Count': float}
            dfpos = dfpos.astype(convert_dict)
            dfpos = dfpos.sort_values(by = ['Count'], ascending=False)
            dfneg = dfneg.astype(convert_dict)
            dfend = dfend.astype(convert_dict)
            dfhke = dfhke.astype(convert_dict)
            dff = pd.concat([dfend, dfhke, dfneg, dfpos])

            names = ['parametro', 'valor']
            dinf = pd.read_csv(getfolderpath(args.folder) / file, names=names, nrows=30, on_bad_lines='skip') #dataframe with lane info
            dinf = dinf.dropna()
            thislane = [] #info for this lane

            #adds id from sample or from file name
            id1 = dinf.loc[dinf['parametro'] == 'ID']
            id1 = id1.loc[:,'valor']
            id1 = id1.iloc[0]
            id1 = str(id1)

            if args.modeid == 'filename':
                id1 = str(file)
            elif args.modeid == 'id+filename':
                id1 += str(file)

            args.current_state = str(' Loading RCC... ' + id1)
            print(args.current_state)
            logging.info(args.current_state)

            if args.autorename == 'on':
                id1 = 'Renamed'
                id1 = str(a) + id1
                id1 = id1.strip()

            comments = dinf.loc[dinf['parametro'] == 'Comments']
            comments =comments.loc[:,'valor']
            if len(comments) >= 1:
                comments = comments.iloc[0]
                comments = str(comments)
            elif len(comments) < 1:
                comments = 'no comments'

            thislane.append(id1)
            thislane.append(comments)

            fovcount = dinf[dinf.parametro == 'FovCount']
            fovcount = fovcount.valor
            fovcounted = dinf[dinf.parametro == 'FovCounted']
            fovcounted = fovcounted.valor
            fovvalue = float(fovcounted)/float(fovcount)

            thislane.append(fovvalue)

            bd = dinf.loc[dinf['parametro'] == 'BindingDensity']
            bd = bd.loc[:,'valor']
            bd = float(bd)
            thislane.append(bd)

            topneg = 3*(np.mean(dfneg['Count']))

            dfnegin = dfneg[dfneg.Count <= topneg]

            background = (np.mean(dfnegin.Count)) + 2*(np.std(dfnegin.Count))
            thislane.append(background)

            background2 = np.max(dfnegin['Count'])
            thislane.append(background2)

            background3 = np.mean(dfnegin['Count'])
            thislane.append(background3)

            negs = list(dfneg['Count'])
            negnames = dfneg['Name']
            negnames = list(negnames)

            maxout = 3*background3
            rownegs = []
            for i in negs:
                rownegs.append(i)
            rownegs.append(maxout)
            dfnegcount[id1] = rownegs

            hkes = list(dfhke['Count'])
            hkenames = dfhke['Name']
            hkenames = list(hkenames)
            rowhkes = []
            for i in hkes:
                rowhkes.append(i)
            dfhkecount[id1] = rowhkes

            dfin = dff[dff.Count >= background]
            dfout = dff[dff.Count < background]
            gbb = len(dfout)
            if a == 0:
                ngen = gbb + len(dfin)

            thislane.append(gbb*100/ngen)
            thislane.append(ngen)

            if a == 0:
                dfgenes['CodeClass'] = dff['CodeClass']
                dfgenes['Name'] = dff['Name']
                dfgenes['Accession'] = dff['Accession']

            dff.set_index('Name', inplace=True)

            diff = None
            if a != 0:
                if set(dff.index) != set(dfgenes.index):
                    diff = list(set(dff.index) - set(dfgenes.index))
                if diff != None:
                    diff.append(list(set(dfgenes.index) - set(dff.index)))
                common = list(set(list(dff.index)).intersection(list(dfgenes.index)))
                dff = dff.loc[common]
                dfgenes = dfgenes.loc[common]
                if diff != None:
                    logging.warning('Mismatch, genes not present in all samples: ' + str(diff))

            dfgenes.index = dff.index
            dfgenes[id1] = dff['Count']
            dfposneg[id1] = dfposneg1

            posgeomean = gmean(dfpos.Count)
            suma = dfpos.Count.sum()
            median = statistics.median(dfpos.Count)
            thislane.append(posgeomean)
            thislane.append(suma)
            thislane.append(median)
            geomeans.append(posgeomean)

            thislogsconc =[] #R2 calculation
            thiscountnolow = []

            #excluding lowest negative control
            y = 0
            while y < 5:
                thiscountnolow.append(dfpos['Count'].iloc[y])
                y = y+1

            #log2 counts
            for i in thiscountnolow:
                thislog = math.log(i,2)
                thislogsconc.append(thislog)

            #r2 score calculaion
            R2 = np.corrcoef(thislogsconc,logconc)
            R2 = R2[1,0]

            thislane.append(R2)

            if args.manualbackground != None:
                background = args.manualbackground
            elif args.background == 'Background2':
                background = background2
            elif args.background == 'Background3':
                background = background3

            lod = background >= dfpos.iloc[4].Count
            thislane.append(lod)

            cerocincofm = dfpos.iloc[4].Count
            thislane.append(cerocincofm)

            infolanes.loc[a] = thislane
            a = a+1

    ##CHECK FOR DUPLICATED IDS
    if all(infolanes.duplicated(subset=['ID']) == False):
        args.current_state = str('--> All ' +  str(len(infolanes['ID'])) + ' IDs are unique, proceeding with analysis. Elapsed %s seconds ' % (time.time() - args.start_time))
        logging.info(args.current_state)
        print(args.current_state)
    elif any(infolanes.duplicated(subset=['ID']) == True):
        args.current_state = str('WARNING ERROR! --> Duplicated IDs, rename samples with unique names or turn on autorename option')
        logging.warning(args.current_state)
        print(args.current_state)

    '''Adding calculated params to infolanes'''
    meangeomeans = np.mean(infolanes['posGEOMEAN'])
    scalingf = []

    if args.manualbackground != None:
        manualbglist = []
        for i in infolanes['ID']:
            manualbglist.append(args.manualbackground)
        infolanes['manual background'] = manualbglist

    for i in infolanes['posGEOMEAN']:
        scaling = meangeomeans/i
        scalingf.append(scaling)
    infolanes['scaling factor'] = scalingf

    dfnegcount = dfnegcount.T
    negnames.append('maxoutlier')
    dfnegcount.columns = negnames

    dfhkecount = dfhkecount.T
    dfhkecount.columns = hkenames

    infolanes.set_index('ID', inplace=True)
    dfgenes.drop('Name', axis=1, inplace=True)

    return infolanes, dfgenes, dfnegcount, dfhkecount, dfposneg

def createoutputfolder(args):
    pathout = str(args.outputfolder)
    pathlib.Path(pathout).mkdir(parents=True, exist_ok=True)
    pathoutimages = str(args.outputfolder) + '/images'
    pathlib.Path(pathoutimages).mkdir(parents=True, exist_ok=True)
    pathoutreports = str(args.outputfolder) + '/reports'
    pathlib.Path(pathoutreports).mkdir(parents=True, exist_ok=True)

def exportrawcounts(rawcounts, args):
    pathout = str(args.outputfolder)
    rawcounts2 = rawcounts
    pathraw = pathout + '/rawcounts.csv'
    rawcounts2.to_csv(pathraw, index=True)
    pathdfraw = pathout + '/dfgenes.csv'
    rawcounts2.to_csv(pathdfraw, index=True)

    rawcounts3 = rawcounts2.drop(['CodeClass', 'Accession'], axis=1)
    pathraw3 = pathout + '/rawcounts2.csv'
    rawcounts3.to_csv(pathraw3, index=True)

def exportdfgenes(dfgenes, args):
    pathout = str(args.outputfolder)
    pathdfgenes = pathout + '/dfgenes.csv'
    dfgenes.to_csv(pathdfgenes, index=True)

def exportrawinfolanes(infolanes, dfnegcount, dfhkecount, dfposneg, args):
    '''
    Exports raw infolanes and dfnegcount, dfhkecount, dfposneg
    '''
    pathout = str(args.outputfolder)
    pathinfolanes = pathout + '/rawinfolanes.csv'
    pathdfnegcount = pathout + '/dfnegcount.csv'
    pathdfhkecount = pathout + '/dfhkecount.csv'

    exportposneg(dfposneg, args)

    infolanes.to_csv(pathinfolanes)
    dfnegcount.to_csv(pathdfnegcount)
    dfhkecount.to_csv(pathdfhkecount)

def pathoutinfolanes(infolanes, args):
    pathout = str(args.outputfolder)
    pathinfolanes = pathout + '/infolanes.csv'
    infolanes.to_csv(pathinfolanes, index=True)

def pathoutrawsummary(rawsummary, args):
    pathout = str(args.outputfolder) + '/reports'
    pathrawsummary = pathout + '/rawsummary.csv'
    rawsummary.to_csv(pathrawsummary, index=True)

def pathoutsummary (summary, args):
    pathout = str(args.outputfolder) + '/reports'
    pathsummary = pathout + '/summary.csv'
    summary.to_csv(pathsummary, index=True)

def condformat_summary(val, top, bot, colorbien = '#a3c771', colorreg = '#f0e986', colormal = '#e3689b'):
    if top >= val >= bot:
        color = colorbien
    elif top*1.15 >= val >= bot*0.85:
        color = colorreg
    elif (bot*0.85 > val) | (val > top*1.15):
        color = colormal

    return 'background-color: {}'.format(color)

def summarizerawinfolanes(args):

    rawinfolanes = pd.read_csv(str(args.outputfolder) + '/rawinfolanes.csv', index_col='ID')

    rawinfofov = [np.min(rawinfolanes['FOV value']), np.max(rawinfolanes['FOV value']), np.mean(rawinfolanes['FOV value']), np.median(rawinfolanes['FOV value'])]
    rawinfobd = [np.min(rawinfolanes['Binding Density']), np.max(rawinfolanes['Binding Density']), np.mean(rawinfolanes['Binding Density']), np.median(rawinfolanes['Binding Density'])]
    rawinfolin = [np.min(rawinfolanes['R2']), np.max(rawinfolanes['R2']), np.mean(rawinfolanes['R2']), np.median(rawinfolanes['R2'])]
    rawinfobackg = [np.min(rawinfolanes['Background']), np.max(rawinfolanes['Background']), np.mean(rawinfolanes['Background']), np.median(rawinfolanes['Background'])]
    rawinfogbb = [np.min(rawinfolanes['Genes below backg %']), np.max(rawinfolanes['Genes below backg %']), np.mean(rawinfolanes['Genes below backg %']), np.median(rawinfolanes['Genes below backg %'])]
    rawinfopgm = [np.min(rawinfolanes['posGEOMEAN']), np.max(rawinfolanes['posGEOMEAN']), np.mean(rawinfolanes['posGEOMEAN']), np.median(rawinfolanes['posGEOMEAN'])]
    rawinfosum = [np.min(rawinfolanes['Sum']), np.max(rawinfolanes['Sum']), np.mean(rawinfolanes['Sum']), np.median(rawinfolanes['Sum'])]
    rawinfo05fm = [np.min(rawinfolanes['0,5fm']), np.max(rawinfolanes['0,5fm']), np.mean(rawinfolanes['0,5fm']), np.median(rawinfolanes['0,5fm'])]
    rawinfoscaf = [np.min(rawinfolanes['scaling factor']), np.max(rawinfolanes['scaling factor']), np.mean(rawinfolanes['scaling factor']), np.median(rawinfolanes['scaling factor'])]

    rawsummary = pd.DataFrame(columns= ['min', 'max', 'mean', 'Median'])
    rawsummary.loc['FOV'] = rawinfofov
    rawsummary.loc['Binding density'] = rawinfobd
    rawsummary.loc['R2'] = rawinfolin
    rawsummary.loc['Background'] = rawinfobackg
    rawsummary.loc['Genes below background'] = rawinfogbb
    rawsummary.loc['posGEOMEAN'] = rawinfopgm
    rawsummary.loc['Sum'] = rawinfosum
    rawsummary.loc['0,5 fm'] = rawinfo05fm
    rawsummary.loc['Scaling factor'] = rawinfoscaf

    pathoutrawsummary(rawsummary, args)
    rawsummary = rawsummary.T
    rawsummary = rawsummary.style.applymap(condformat_summary, top=args.maxfov, bot=args.minfov, subset='FOV')
    rawsummary = rawsummary.applymap(condformat_summary, top = args.maxbd, bot=args.minbd, subset='Binding density')
    rawsummary = rawsummary.applymap(condformat_summary, top= args.maxlin, bot=args.minlin,  subset='R2')
    rawsummary = rawsummary.applymap(condformat_summary, top= args.pbelowbackground, bot=0,  subset='Genes below background')
    rawsummary = rawsummary.applymap(condformat_summary, top= args.maxscalingfactor, bot=args.minscalingfactor,  subset='Scaling factor')


    rawsummary.to_html(str(args.outputfolder) + '/rawsummary.html')


    if args.showbrowserrawqc == True:
        webbrowser.open(str(args.outputfolder) + '/rawsummary.html')




def summarizeinfolanes(args):
    infolanes = pd.read_csv(str(args.outputfolder) + '/infolanes.csv', index_col='ID')
    infofov = [np.min(infolanes['FOV value']), np.max(infolanes['FOV value']), np.mean(infolanes['FOV value']), np.median(infolanes['FOV value'])]
    infobd = [np.min(infolanes['Binding Density']), np.max(infolanes['Binding Density']), np.mean(infolanes['Binding Density']), np.median(infolanes['Binding Density'])]
    infolin = [np.min(infolanes['R2']), np.max(infolanes['R2']), np.mean(infolanes['R2']), np.median(infolanes['R2'])]
    infobackg = [np.min(infolanes['Background']), np.max(infolanes['Background']), np.mean(infolanes['Background']), np.median(infolanes['Background'])]
    infogbb = [np.min(infolanes['Genes below backg %']), np.max(infolanes['Genes below backg %']), np.mean(infolanes['Genes below backg %']), np.median(infolanes['Genes below backg %'])]
    infopgm = [np.min(infolanes['posGEOMEAN']), np.max(infolanes['posGEOMEAN']), np.mean(infolanes['posGEOMEAN']), np.median(infolanes['posGEOMEAN'])]
    infosum = [np.min(infolanes['Sum']), np.max(infolanes['Sum']), np.mean(infolanes['Sum']), np.median(infolanes['Sum'])]
    info05fm = [np.min(infolanes['0,5fm']), np.max(infolanes['0,5fm']), np.mean(infolanes['0,5fm']), np.median(infolanes['0,5fm'])]
    infoscaf = [np.min(infolanes['scaling factor']), np.max(infolanes['scaling factor']), np.mean(infolanes['scaling factor']), np.median(infolanes['scaling factor'])]

    summary = pd.DataFrame(columns= ['min', 'max', 'mean', 'Median'])
    summary.loc['FOV'] = infofov
    summary.loc['Binding density'] = infobd
    summary.loc['R2'] = infolin
    summary.loc['Background'] = infobackg
    summary.loc['Genes below background'] = infogbb
    summary.loc['posGEOMEAN'] = infopgm
    summary.loc['Sum'] = infosum
    summary.loc['0,5 fm'] = info05fm
    summary.loc['Scaling factor'] = infoscaf

    summary = summary.T
    summary2view = summary.style.applymap(condformat_summary, top=args.maxfov, bot=args.minfov, subset='FOV')
    summary2view = summary2view.applymap(condformat_summary, top=args.maxbd, bot=args.minbd, subset='Binding density')
    summary2view = summary2view.applymap(condformat_summary, top=args.maxlin, bot=args.minlin, subset='R2')
    summary2view = summary2view.applymap(condformat_summary, top= args.pbelowbackground, bot=0,
                                     subset='Genes below background')
    summary2view = summary2view.applymap(condformat_summary, top=args.maxscalingfactor, bot=args.minscalingfactor,
                                     subset='Scaling factor')

    pathoutsummary(summary, args)
    summary2view.to_html(str(args.outputfolder) + '/Summary.html')

    if args.showbrowserqc == True:
        webbrowser.open(str(args.outputfolder) + '/Summary.html')

def exportposneg(dfposneg, args):
    posnegcounts = pd.DataFrame()

    for i in dfposneg.keys():
        a = dfposneg[i]
        posnegcounts['Name'] = a['Name']
        posnegcounts['CodeClass'] = a['CodeClass']
        posnegcounts[i] = a['Count']

    posnegcounts.set_index('Name', drop=True, inplace=True)

    pathout = str(args.outputfolder)
    pathposneg = pathout + '/posnegcounts.csv'
    posnegcounts.to_csv(pathposneg, index=True)


def plotfovvalue(args, infolanes):

    minfov = []
    maxfov = []
    for i in infolanes.index:
        minfov.append(args.minfov)
        maxfov.append(args.maxfov)
    plt.figure()
    plt.plot(infolanes.index, infolanes['FOV value'], 'bo')
    plt.plot(minfov, 'r', label='min')
    plt.plot(maxfov, 'g', label='optimal')
    plt.legend()
    plt.ylabel('fov value')
    plt.xlabel('samples')
    plt.title('IMAGE QC (FOV)')
    plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
    plt.grid(True)
    plt.savefig(str(args.outputfolder) + '/images/fovplot.png')

def plotbd(args, infolanes):
    minbd = []
    maxbd = []
    for i in infolanes.index:
        minbd.append(args.minbd)
        maxbd.append(args.maxbd)
    plt.figure()
    plt.plot(infolanes.index, infolanes['Binding Density'] ,'bo')
    plt.plot(infolanes.index, minbd, color='m', label='min')
    plt.plot(infolanes.index, maxbd, color='r', label = 'max')
    plt.title('Binding Density')
    plt.xlabel('samples')
    plt.ylabel('binding density')
    plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
    plt.legend()
    plt.grid(True)
    plt.savefig(str(args.outputfolder) + '/images/bdplot.png')

def plotgenbackground(args, infolanes):
    ngenlist = []
    ngen = infolanes['nGenes'][0]
    for i in infolanes.index:
        ngenlist.append(ngen)
    plt.figure()
    plt.bar(infolanes.index, infolanes['nGenes'] - infolanes['Genes below backg %'])
    plt.plot(infolanes.index, ngenlist, 'ro', label='total genes')
    plt.legend()
    plt.xticks(rotation=45)
    plt.xlabel('ID')
    plt.ylabel('genes')
    plt.title('Genes above background')
    plt.savefig(str(args.outputfolder) + '/images/genbackground.png')

def plotld(args, infolanes):
    plt.figure()
    plt.plot(infolanes.index, infolanes['0,5fm'], 'bo', label ='0,5fm')

    if args.manualbackground is not None:
        background = 'manual background'
    else:
        background = args.background
        if background == 'Backgroundalt':
            background = 'Background'

    plt.plot(infolanes.index, infolanes[background], 'r', label='Background')
    plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
    plt.title('Limit of detection')
    plt.xlabel('samples')
    plt.ylabel('0,5 fm')
    plt.legend()
    plt.grid(True)
    plt.savefig(str(args.outputfolder) + '/images/ldplot.png')

def plotocn(args, infolanes, dfnegcount):

    plt.figure()

    if args.manualbackground != None:
        background = 'manual background'
    else:
        background = args.background
        if background == 'Backgroundalt':
            background = 'Background'

    for i in dfnegcount.columns:
        plt.plot(dfnegcount.index, dfnegcount[i], 'o', label=i,)
    plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
    plt.plot(dfnegcount.index, dfnegcount['maxoutlier'])
    plt.plot(infolanes[background], label='Background')
    plt.xlabel('samples')
    plt.ylabel('counts')
    plt.legend(loc='lower right', ncol=4, mode='expand')
    plt.title('Outliers in neg_controls')
    plt.savefig(str(args.outputfolder) + '/images/ocnplot.png')

def plotlin(args, infolanes):
    minlin = []
    optlin = []
    for i in infolanes.index:
        minlin.append(args.minlin)
        optlin.append(args.maxlin)

    ngenlist = []
    ngen = infolanes['nGenes'][0]
    for i in infolanes.index:
        ngenlist.append(ngen)
    plt.figure()
    plt.bar(infolanes.index, (infolanes['nGenes'] - infolanes['Genes below backg %'])/infolanes['nGenes'], color='cyan')
    plt.plot(infolanes.index, infolanes['R2'], 'o' ,color='blue')
    plt.plot(infolanes.index, minlin, 'm')
    plt.plot(infolanes.index, optlin, 'g')
    plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
    plt.xlabel('samples')
    cyanbar = mpatches.Patch(color='cyan', label='% genes > background')
    purpleline = mpatches.Patch(color='m', label='min value')
    bluedot = mpatches.Patch(color='blue', label='R2 value')
    plt.legend(handles=[cyanbar, purpleline, bluedot], loc='lower right', ncol=3, mode='expand')
    plt.ylabel('genes')
    plt.xticks(rotation=45)
    plt.title('Linearity and genes above background')
    plt.savefig(str(args.outputfolder) + '/images/linplot.png')

def plothke(args, infolanes, dfhkecount):
    '''Housekeeping plot'''
    plt.figure()
    for i in dfhkecount.columns:
        plt.plot(dfhkecount.index,dfhkecount[i], 'o', label=i,)
    plt.xlabel('samples')
    plt.ylabel('counts')
    plt.legend(loc='upper left', ncol=3, mode='expand')
    plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
    plt.title('Housekeeping genes')
    plt.savefig(str(args.outputfolder) + '/images/hkeplot.png')


def plothkel(args, infolanes, dfhkecount):
    '''Closest housekeeping to background plot'''
    bb = np.mean(infolanes['Background'])

    bbmax = 6*bb
    bblist = []

    for i in infolanes.index:
        bblist.append(bb)

    hkelplot = plt.figure()
    ax1 = hkelplot.add_subplot(111)
    number_of_plots = len(dfhkecount.columns)
    colors = sns.color_palette("hls", number_of_plots)
    ax1.set_prop_cycle('color', colors)
    for i in dfhkecount.columns:
        ax1.plot(dfhkecount.index, dfhkecount[i], 'o', label=i,)
    plt.plot(infolanes.index, bblist, 'r')
    plt.xlabel('ID')
    plt.ylabel('counts')
    plt.ylim(0,2*bbmax)
    plt.legend(loc='upper left', ncol=3, mode='expand')
    plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
    plt.title('Housekeeping genes close to background')
    plt.savefig(str(args.outputfolder) + '/images/hkelplot.png')

def plotsca(args, infolanes):
    scalingflist = infolanes['scaling factor']
    slmin = []
    slmax = []
    for i in scalingflist:
        slmin.append(args.minscalingfactor)
        slmax.append(args.maxscalingfactor)

    plt.figure()
    plt.plot(infolanes.index, infolanes['scaling factor'], 'o')
    plt.plot(infolanes.index, slmin, 'm', label='min')
    plt.plot(infolanes.index, slmax, 'r', label='max')
    plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
    redline = mpatches.Patch(color='red', label='max scaling factor')
    purpleline = mpatches.Patch(color='m', label='min scaling factor')
    bluedot = mpatches.Patch(color='blue', label='sample scaling factor')
    plt.legend(handles=[redline, purpleline, bluedot], loc='lower right', ncol=2, mode='expand')
    plt.xlabel('samples')
    plt.ylabel('scaling factor')
    plt.title('scaling factor')
    plt.savefig(str(args.outputfolder) + '/images/scaplot.png')

def pdfreport(args):

    pdf = FPDF()
    pdf.add_page()
    pdf.set_font('Arial', 'B', 16)

    pdf.image(str(pathlib.Path(__file__).parent) + '/reports/images/qc_template_report.png',0,0,h=297)

    pdf.image(str(args.outputfolder) + '/images/ldplot.png', 12.5, 42, h=69)
    pdf.image(str(args.outputfolder) + '/images/bdplot.png', 110, 42, h=69)

    pdf.image(str(args.outputfolder) + '/images/fovplot.png', 10.5, 120, h=69)
    pdf.image(str(args.outputfolder) + '/images/linplot.png', 10.5, 200, h=69)
    pdf.image(str(args.outputfolder) + '/images/hkelplot.png', 110, 120, h=69)
    pdf.image(str(args.outputfolder) + '/images/ocnplot.png', 110, 200, h=69)

    pdf.output(str(args.outputfolder) + '/reports/QC_inspection.pdf', 'F')

    if args.showbrowserqc == True:
        os.system(str(args.outputfolder) + '/reports/QC_inspection.pdf')

def pdfreportnorm(args):
    pdf = FPDF()
    pdf.add_page()
    pdf.set_font('Arial', 'B', 16)

    pdf.image(str(pathlib.Path(__file__).parent) + '/reports/images/qc_template_report.png',0,0,h=297)

    pdf.image(str(args.outputfolder) + '/images/avgm.png', 12.5, 42, h=69)
    pdf.image(str(args.outputfolder) + '/images/uve.png', 110, 42, h=69)

    pdf.image(str(args.outputfolder) + '/images/rlerawplot.png', 10.5, 119, h=77.5)
    pdf.image(str(args.outputfolder) + '/images/rlenormplot.png', 10.5, 199, h=77.5)

    pdf.output(str(args.outputfolder) + '/reports/norm_report.pdf', 'F')

    if args.showbrowserqc == True:
        os.system(str(args.outputfolder) + '/reports/QC_inspection.pdf')

def flagqc(args):
    infolanes = pd.read_csv(str(args.outputfolder) + '/rawinfolanes.csv', index_col='ID')

    flagged = set([])
    if os.path.exists(str(args.outputfolder) + '/reports/QCflags.txt'):
        os.remove(str(args.outputfolder) + '/reports/QCflags.txt')

    f = open(str(args.outputfolder) + '/reports/QCflags.txt', 'a')
    args.current_state = '--> Starting QC flagging:'
    logging.info(args.current_state)
    print(args.current_state)

    if (
        all(i >= args.minfov for i in infolanes['FOV value']) and
        all( args.maxbd >= float(i) >= args.minbd for i in infolanes.loc[:,'Binding Density']) and
        all(i == False for i in infolanes.loc[:,'limit of detection']) and
        all(i >= j for i,j in zip(infolanes.loc[:,'0,5fm'],infolanes.loc[:,'Background'])) and
        all(args.maxscalingfactor > float(i) > args.minscalingfactor for i in infolanes.loc[:, 'scaling factor']) and
        all(i <= args.pbelowbackground for i in infolanes.loc[:,'Genes below backg %'])):
        info = 'Appropiate QC values for all samples. \n'
        args.current_state = info
        print(args.current_state)
        logging.info(args.current_state)
        f.writelines(info + '\n')
    else:
        args.current_state = 'Inappropiate QC values in some samples, revise QC report'
        logging.warning(args.current_state)
        print(args.current_state)
        for i in infolanes.index:
            thisFOV = infolanes.at[i,'FOV value']
            thisBD = infolanes.at[i, 'Binding Density']
            thisLOD = infolanes.at[i, 'limit of detection']
            thisBG = infolanes.at[i, 'Background']
            this05 = infolanes.at[i, '0,5fm']
            thisSF = infolanes.at[i,'scaling factor']
            thisgbb = infolanes.at[i,'Genes below backg %']
            if thisFOV < args.minfov:
                fovinfo = 'Low FOV value (' + str(thisFOV) +') in ' + i + '. Sample is flagged/discarded. \n'
                print(fovinfo)
                logging.warning(fovinfo)
                f.writelines(fovinfo)
                flagged.add(i)
            if thisBD > args.maxbd or thisBD < args.minbd:
                bdinfo = 'Wrong binding density value (' + str(thisBD) +') in ' + i + '. Sample is flagged/discarded.\n'
                print(bdinfo)
                logging.warning(bdinfo)
                f.writelines(bdinfo)
                flagged.add(i)
            if thisLOD == True:
                lodinfo = 'Wrong limit of detection value (' + str(thisLOD) +') in ' + i + '. Sample is flagged/discarded.\n'
                print(lodinfo)
                logging.warning(lodinfo)
                f.writelines(lodinfo)
                flagged.add(i)
            if thisBG > this05:
                bginfo = 'Wrong 0.5fm value (' + str(thisBG) +') in ' + i + '. Sample is flagged/discarded.\n'
                print(bginfo)
                logging.warning(bginfo)
                f.writelines(bginfo)
                flagged.add(i)
            if thisSF > 3 or thisSF < 0.3:
                sfinfo = 'Wrong scaling factor value (' + str(thisSF) +') in ' + i + '. Sample is flagged/discarded.\n'
                print(sfinfo)
                logging.warning(sfinfo)
                f.writelines(sfinfo)
                flagged.add(i)
            if thisgbb > args.pbelowbackground:
                gbbinfo = 'Wrong genes below background value (' + str(thisgbb) +') in ' + i + '. Sample is flagged/discarded.\n'
                print(gbbinfo)
                logging.warning(gbbinfo)
                f.writelines(gbbinfo)
                flagged.add(i)
    f.close()
    flaggeddf = pd.DataFrame(flagged, columns=['Flagged_samples'])
    flaggeddf.to_csv(str(args.outputfolder) + '/flagged.csv', index=False)

    if len(flagged) >= 3:
        args.badlanes = str(len(flagged)) + ' badlanes detected, check output/reports/QCflags.txt'
    elif 0 < len(flagged) < 3:
        args.badlanes = str(len(flagged)) + ' badlanes detected: ' + str(flagged)
    elif len(flagged) == 0:
        args.badlanes = 'No bad lanes detected from QC'

    return flagged

def removelanes(autoremove, args):
    infolanes = pd.read_csv(str(args.outputfolder) + '/rawinfolanes.csv', index_col='ID')
    dfgenes = pd.read_csv(str(args.outputfolder) + '/dfgenes.csv', index_col='Name')
    manualremove = args.remove
    if autoremove == None:
        autoremove = set(manualremove)
    if manualremove != None:
        if type(manualremove) == str:
            manualremove = set(manualremove.split())
        print('Se retiran manualmente las muestras:', manualremove, '.')
        logging.info('Se retiran manualmente las muestras:' + str(manualremove) + '.')
        autoremove.update(manualremove)
    dfgenes = dfgenes.T

    if set(dfgenes.index.intersection(autoremove)) == autoremove:
        dfgenes = dfgenes.T
        dfgenes11 = dfgenes.drop(autoremove, axis=1)
        infolanes = infolanes.drop(autoremove)
    else:
        dfgenes11 = dfgenes.T
        args.current_state = 'Lanes set to remove not present in analysis.'
        print('Error: ' + args.curent_state)
        logging.error(args.current_state)

    pathout = str(args.outputfolder)
    pathinfolanes = pathout + '/infolanes.csv'
    infolanes.to_csv(pathinfolanes, index=True)
    exportdfgenes(dfgenes, args)

    return dfgenes11, infolanes

def exportfilrawcounts(rawfcounts, args):
    pathout = str(args.outputfolder)
    rawfcounts2 = rawfcounts
    pathfraw = pathout + '/rawfcounts.csv'
    rawfcounts2.to_csv(pathfraw, index=True)

def rescalingfactor23(args):
    """Scaling factor needs to be recalculated after removing samples excluded by QC inspection"""
    infolanes = pd.read_csv(str(args.outputfolder) + '/infolanes.csv', index_col=0)

    if args.tecnormeth == 'posgeomean' or args.tecnormeth == 'regression':
        use = 'posGEOMEAN'
        if args.tecnormeth == 'regression':
            negs = pd.read_csv(str(args.outputfolder) + '/dfnegcount.csv', index_col=0)
            negs = negs.drop('maxoutlier', axis=1)
            corrected_negs = regretnegs(negs, args)
            backgr_regr = []
            for i in corrected_negs.index:
                thisbackg_regr = (np.mean(corrected_negs.loc[i])) + 2*(np.std(corrected_negs.loc[i]))
                backgr_regr.append(thisbackg_regr)

            infolanes['backgr_regr'] = backgr_regr

    elif args.tecnormeth == 'Sum':
        use = 'Sum'
    elif args.tecnormeth == 'Median':
        use = 'Median'
    ref = np.mean(infolanes[use])


    scalingf2 = []

    for i in infolanes[use]:
        scaling = ref/i
        scalingf2.append(scaling)

    infolanes['scaling factor2'] = scalingf2

    return infolanes

def reinfolanes(args):
    '''
    Whether background correction (transform low counts) or technorm happens first, background or scaling factor needs to be recalculated after.
    For background correction, new background shoult be calculated from dfgenes given by technorm
    For technorm, new scaling factor needs to be calculated from dfgenes given by transformlowcounts
    '''

    infolanes = findaltnegatives(args)
    fildfgenes = pd.read_csv(str(args.outputfolder) + '/dfgenes.csv', index_col='Name')
    rawdfgenes = pd.read_csv(str(args.outputfolder) + '/rawfcounts.csv', index_col='Name')

    if args.firsttransformlowcounts == True:
        dfgenes = rawdfgenes
    elif args.firsttransformlowcounts == False:
        dfgenes = fildfgenes


    dfneg = dfgenes[dfgenes['CodeClass'] == 'Negative'].drop(['CodeClass', 'Accession'], axis=1).T
    infolanes['Background'] = dfneg.mean(axis=1)+(2*(dfneg.std(axis=1)))
    infolanes['Background2'] = dfneg.max(axis=1)
    infolanes['Background3'] = dfneg.mean(axis=1)

    dfpos = fildfgenes[fildfgenes['CodeClass'] == 'Positive'].drop(['CodeClass', 'Accession'], axis=1).T
    infolanes['posGEOMEAN'] = gmean(dfpos, axis=1)
    infolanes['Sum'] = dfpos.sum(axis=1)
    infolanes['Median'] = np.median(dfpos, axis=1)

    dfgenes4mean = dfgenes.drop(['CodeClass', 'Accession'], axis = 1).T
    infolanes['meanexpr'] = gmean(dfgenes4mean, axis=1)

    pathoutinfolanes(infolanes, args)
    infolanes = rescalingfactor23(args)
    pathoutinfolanes(infolanes,args)

    dfgenes = dfgenes.drop(list(dfpos.columns))
    dfgenes = dfgenes.drop(list(dfneg.columns))
    exportdfgenes(dfgenes, args)


def regretnegs(negs, args):
    corrected_negs = pd.DataFrame()
    posneg = pd.read_csv(str(args.outputfolder) + '/posnegcounts.csv', index_col=0)
    for i in posneg.index:
        if 'NEG' in i:
            posneg = posneg.drop(i, axis=0)

    posneg2 = posneg.copy()
    posneg2 = posneg2.drop('CodeClass', axis=1)
    gmeans = []
    for i in posneg2.index:
        thisgmean = gmean(posneg2.loc[i])
        gmeans.append(thisgmean)

    posneg['mean'] = gmeans
    xmean = sorted(posneg['mean'])
    ymean = [31, 125, 500, 2000, 8000, 32000]

    meaneq_abc = np.polynomial.Polynomial.fit(xmean, ymean, 3)

    infolanes = pd.read_csv(str(args.outputfolder) + '/infolanes.csv')
    for i in infolanes['ID']:
        thisx = sorted(posneg[i])
        thiseq_abc = np.polynomial.Polynomial.fit(thisx, ymean, 3)

        def correct(thisseq_abc, meaneq_abc, counts):
            this = thisseq_abc(counts)
            solutions = (meaneq_abc - this).roots()
            corrected = abs(min(solutions, key=lambda x: abs(x - counts)))
            return corrected

        thissample = []

        newnegs = negs.T

        corrected_negs.index = newnegs.index

        for k in newnegs[str(i)]:
            corrected = correct(thisseq_abc=thiseq_abc, meaneq_abc=meaneq_abc, counts=k)

            thissample.append(corrected)

        corrected_negs[str(i)] = thissample

    return corrected_negs.T

def findaltnegatives(args):
    '''To find low and stably expressed genes through the endogenous, to use as negative controls
    in case native neg controls are not robust
    Generates new infolanes with background alt in it'''

    infolanes = pd.read_csv(str(args.outputfolder) + '/infolanes.csv', index_col=0)
    dfgenes = pd.read_csv(str(args.outputfolder) + '/rawfcounts.csv', index_col='Name')
    dfgenes.drop(['CodeClass', 'Accession'], inplace=True, axis=1)
    dfgenes = dfgenes.T

    genmean = dfgenes.mean()
    meangenmean = np.mean(genmean)
    genmean = genmean/meangenmean

    genmean.sort_values(inplace=True)

    genstd = dfgenes.std()
    meangenstd = np.mean(genstd)
    genstd = genstd/meangenmean
    genstd = genstd*2

    genstd.sort_values(inplace=True)

    genrank = genmean * genstd
    genrank.sort_values(inplace=True)

    bestaltnegs = genrank.head(10)
    bestaltnegsnames = bestaltnegs.index

    dfaltnegs = dfgenes[bestaltnegsnames]

    dfaltnegs = dfaltnegs.T
    backgroundaltmean = dfaltnegs.mean()
    backgroundalt2std = dfaltnegs.std()*2
    backgroundalt = backgroundaltmean + backgroundalt2std

    infolanes['Backgroundalt'] = backgroundalt

    return infolanes

def normtecnica(dfgenes, args):

    infolanes = pd.read_csv(str(args.outputfolder) + '/infolanes.csv')

    normgenes = pd.DataFrame()
    normgenes['CodeClass'] = dfgenes['CodeClass']
    normgenes['Name'] = dfgenes['Name']
    normgenes['Accession'] = dfgenes['Accession']

    for i in infolanes['ID']:
        j = float(infolanes.loc[infolanes['ID'] == str(i), 'scaling factor2'])

        this = []

        for k in dfgenes[str(i)]:
            m = j*k
            this.append(m)

        normgenes[str(i)] = this

    normgenes.set_index('Name', drop=True, inplace=True)

    return normgenes

def regresion(dfgenes, args):
    normgenes = pd.DataFrame()
    normgenes['CodeClass'] = dfgenes['CodeClass']
    normgenes['Name'] = dfgenes['Name']
    normgenes['Accession'] = dfgenes['Accession']

    posneg = pd.read_csv(str(args.outputfolder) + '/posnegcounts.csv', index_col=0)
    for i in posneg.index:
        if 'NEG' in i:
            posneg = posneg.drop(i, axis=0)

    posneg2 = posneg.copy()
    posneg2 = posneg2.drop('CodeClass', axis=1)
    gmeans = []
    for i in posneg2.index:

        thisgmean = gmean(posneg2.loc[i])
        gmeans.append(thisgmean)

    posneg['mean'] = gmeans

    xmean = sorted(posneg['mean'])
    ymean = [31, 125, 500, 2000, 8000, 32000]

    meaneq_abc = np.polynomial.Polynomial.fit(xmean, ymean, 3)

    infolanes = pd.read_csv(str(args.outputfolder) + '/infolanes.csv')
    for i in infolanes['ID']:
        thisx = sorted(posneg[i])
        thiseq_abc = np.polynomial.Polynomial.fit(thisx, ymean, 3)

        def correct(thisseq_abc, meaneq_abc, counts):
            this = thisseq_abc(counts)
            solutions = (meaneq_abc - this).roots()
            corrected = abs(min(solutions, key=lambda x: abs(x - counts)))
            return corrected

        thissample = []

        for k in dfgenes[str(i)]:
            corrected = correct(thisseq_abc=thiseq_abc, meaneq_abc=meaneq_abc, counts=k)

            thissample.append(corrected)

        normgenes[str(i)] = thissample

    normgenes.set_index('Name', drop=True, inplace=True)
    return normgenes


def transformlowcounts(args):

    dfgenes = pd.read_csv(str(args.outputfolder) + '/dfgenes.csv', index_col='Name')
    infolanes = pd.read_csv(str(args.outputfolder) + '/infolanes.csv')

    ilanes = infolanes.T

    ilanes.columns = ilanes.iloc[0]
    ilanes = ilanes.drop(['ID'])
    varback = args.lowcounts
    varbg = args.background

    # if args.tecnormeth == 'regression':
    #     varbg = 'backgr_regr'

    if args.manualbackground != None:
        mvarbg = args.manualbackground

        if varback == 'skip':
            pass

        elif varback == 'asim':
            for i in infolanes['ID']:
                estebg = mvarbg
                dfgenes.loc[dfgenes[i] <= estebg, i] = estebg
        elif varback == 'sustract':
            for i in infolanes['ID']:
                estebg = mvarbg
                dfgenes.loc[dfgenes[i] <= estebg, i] = 0
                dfgenes.loc[dfgenes[i] > estebg, i] = dfgenes[i] - estebg
                dfgenes.replace(0,1,inplace=True)

    else:
        if varback == 'skip':
            pass
        elif varback == 'asim':
            for i in infolanes['ID']:
                estebg = ilanes.at[varbg,i]
                dfgenes.loc[dfgenes[i] <= estebg, i] = estebg
        elif varback == 'sustract':
            for i in infolanes['ID']:
                estebg = ilanes.at[varbg,i]
                dfgenes.loc[dfgenes[i] <= estebg, i] = 0
                dfgenes.loc[dfgenes[i] > estebg, i] = dfgenes[i] - estebg
                dfgenes.replace(0,1,inplace=True)

    exportdfgenes(dfgenes, args)

    return dfgenes

def exporttnormgenes(normgenes, args):
    pathout = str(args.outputfolder)
    pathnormgenes = pathout + '/tnormcounts.csv'
    normgenes.to_csv(pathnormgenes, index=True)

def getallhkes(args):
    dfgenes = pd.read_csv(str(args.outputfolder) + '/dfgenes.csv', index_col='Name')

    allhkes = dfgenes.loc[dfgenes.loc[:,'CodeClass'] == 'Housekeeping']

    return allhkes

def filter50chkes(allhkes, args):
    '''Filters housekeeping genes with less than 50 counts'''
    infolanes = pd.read_csv(str(args.outputfolder) + '/infolanes.csv')
    selhkes = pd.DataFrame()

    for i in infolanes['ID']:
            selhkes = allhkes.loc[allhkes[i] >= args.mincounthkes, :]
            sumatorio = allhkes[i].sum()

    selhkes = selhkes.round(decimals=3)
    selhkes = selhkes.drop(['CodeClass', 'Accession'], axis=1)

    return selhkes

def findrefend(args, selhkes):
    '''Finds endogenous that can be used as reference genes'''

    dfgenes = pd.read_csv(str(args.outputfolder) + '/tnormcounts.csv')

    norm2end = dfgenes.loc[dfgenes['CodeClass'] == 'Endogenous']
    norm2end1 = dfgenes.loc[dfgenes['CodeClass'] == 'Endogenous1']
    norm2end = pd.concat([norm2end,norm2end1])
    norm2end = norm2end.drop(['CodeClass','Accession'], axis='columns')

    norm2end2 = norm2end.set_index('Name')

    if args.refendgenes == 'endhkes':

        endge = FindERG(norm2end2)

        bestend = endge[0:args.numend] #n best endogenous to include as reference genes
        logging.info('Most promising endogenous genes: ' +  str(bestend))
        print('Most promising endogenous genes: ', bestend)
    refgenes = selhkes
    # if 'CodeClass' in refgenes.columns:
    #     refgenes.drop('CodeClass', axis='columns', inplace=True)
    # if 'Accession' in refgenes.columns:
    #     refgenes.drop('Accession', axis='columns', inplace=True)

    if args.refendgenes == 'endhkes':
        for i in bestend:
            isbest = norm2end.loc[:,'Name'] == i
            refgen = norm2end.loc[isbest]
            refgen.set_index('Name', drop=True, inplace=True)
            refgenes = pd.concat([refgenes, refgen])

    return refgenes

def pathoutrefgenes(refgenes, args):
    pathout = str(args.outputfolder)
    pathrefgenesview = pathout + '/refgenesview.csv'
    refgenes.to_csv(pathrefgenesview, header=True, index=True)
    refgenes = refgenes.T
    pathrefgenes = pathout + '/refgenes.csv'
    refgenes.to_csv(pathrefgenes, header=True, index=True)

def getgroups(args):
    refgenes = pd.read_csv(str(args.outputfolder) + '/refgenes.csv', index_col=0)
    flagged = pd.read_csv(str(args.outputfolder) + '/flagged.csv')
    flagged = set(flagged['Flagged_samples'])
    dfgroups = pd.read_csv(args.groupsfile, header=0, index_col=0)
    if args.laneremover == 'yes':
        for i in flagged:
            if i in dfgroups.columns:
                dfgroups.drop(i, axis=0, inplace=True)
    groups = list(set(dfgroups['GROUP']))

    d = {}
    ddf = {}
    for i in groups:
        a = list(dfgroups[dfgroups['GROUP'] == i].index)
        d[i] = a
        newa = []
        for j in a:
            if j in refgenes.index:
                newa.append(j)
        ddf[i] = refgenes.loc[newa]

    return ddf

def calkruskal(*args):
    '''Kruskal wallis calculation
    Takes dfa-like dataframes, groups with samples at y and ref genes at x'''
    gencount = 0

    lk = {}
    for i in args[0]:
        la = []
        a = args[0].iloc[:,gencount]

        gcount = 0
        for j in args:
            b = args[gcount].iloc[:,gencount]
            la.append(b)
            gcount+=1
        try:
            krus = stats.kruskal(*la)
        except Exception:
            pass
        lk[a.name] = krus
        gencount +=1

    lk = pd.DataFrame.from_dict(lk)
    lk = lk.rename(index={0:'Result', 1: 'pvalue'})

    return lk

def calwilco(dfa,dfb):
    '''Calculates wilcoxon for every pair of groups'''
    count = 0
    lw = {}
    for i in dfa:
        a = dfa.iloc[:,count]
        b = dfb.iloc[:,count]
        k = stats.ranksums(a,b)
        lw[i] = k
        count += 1
    lw = pd.DataFrame.from_dict(lw)
    lw = lw.rename(index={0:'Result', 1: 'pvalue'})
    return lw

def calwilcopairs(*ddfc):
    name = ddfc[0][0]
    df = ddfc[0][1]


    lenargs = np.arange(0,len(ddfc))

    lw = {}
    for i in lenargs:
        for j in lenargs:
            if i < j:
                w = calwilco(ddfc[i][1], ddfc[j][1])
                pair = 'wilcox: ' + str(ddfc[i][0]) + '/' + str(ddfc[j][0])
                lw[pair] = w
    return lw

def flagkrus(reskrus):
    flaggedgenes = []
    for i in reskrus:
        if reskrus.loc['pvalue',i] < 0.05:
            flaggedgenes.append(i)
    return flaggedgenes

def filterkruskal(flaggedgenes, args):
    refgenes = pd.read_csv(str(args.outputfolder) + '/refgenes.csv', index_col=0)
    if args.filtergroupvariation == 'filterkrus':
        if (len(refgenes.columns) - len(flaggedgenes)) <=2:
            args.current_state = 'Too much genes to be removed from kruskal filtering, consider using another refgenes or change settings to "flagkrus".'
            print(args.current_state)
            logging.warning(args.current_state)
        else:
            refgenes = refgenes.drop(columns=flaggedgenes)
    elif args.filtergroupvariation == 'flagkrus':
        args.current_state = str('Genes not recommended as refgenes by kruskal: ' + str(flaggedgenes) + '.')
        logging.warning(args.current_state)
        print(args.current_state)
    pathout = str(args.outputfolder)
    pathrefgenes = pathout + '/refgenes.csv'
    refgenes.to_csv(pathrefgenes, header=True, index=True)
    return refgenes

def flagwilcox(reswilcopairs):
    flaggedwilcox = []
    for i in reswilcopairs.values():
        for j in i:
            if i.loc['pvalue',j] < 0.05:
                flaggedwilcox.append(j)
    flaggedwilcox = set(flaggedwilcox)
    return flaggedwilcox

def filterwilcox(flaggedwilcox, args):
    refgenes = pd.read_csv(str(args.outputfolder) + '/refgenes.csv', index_col=0)
    if args.filtergroupvariation == 'filterwilcox':
        if len(flaggedwilcox) < len(refgenes.columns) and (len(refgenes.columns) - len(flaggedwilcox)) > 2:
            refgenes = refgenes.drop(columns=flaggedwilcox)
        else:
            args.current_state = 'Too many ref genes filtered by wilcox, skipping filtering. Consider re-designing reference/housekeeping genes.'
            print(args.current_state)
            logging.warning(args.current_state)
    elif args.filtergroupvariation == 'flagwilcox':
        args.current_state = 'Genes not recommended as refgenes by wilcoxon: ' + str(flaggedwilcox) + '.'
        print(args.current_state)
        logging.warning(args.current_state)
    pathout = str(args.outputfolder)
    pathrefgenes = pathout + '/refgenes.csv'
    refgenes.to_csv(pathrefgenes, header=True, index=True)
    return refgenes

def measureM(df, ctVal=False):
    if ctVal == True:
        for column in df:
            minimo = min(df[column])
            df[column] = df[column] - minimo
        df = 2**-df
    else:
        for column in df:
            maximo = max(df[column])
            df[column] = df[column] / maximo
    m=list(df.index)

    n=list(df.columns)

    M_a = pd.DataFrame()
    for j in n:
        Vjk= []
        for k in n:
            Ajk=[]
            for i in m:
                Ajk.append(np.log2(df.at[i,j]/df.at[i,k]))
            Vjk.append(np.std(Ajk, ddof=1))

        M_a.at[j,'Genes'] = j

        M_a.at[j, 'M'] = np.sum(Vjk)/(len(n)-1)
    M_a = M_a.sort_values('M', ascending=False)

    return M_a

def geNorm(df, avgm=pd.DataFrame()):
    result = measureM(df)
    n = len(df.columns)
    if n <= 2:
        bestgen = result.iat[0,0]
        newrow = pd.DataFrame([[bestgen, result['M'].mean()]])
        avgm = pd.concat([avgm, newrow])
        lastrow = pd.DataFrame([[result.iat[1,0], result.iat[1,1]]])
        avgm = pd.concat([avgm, lastrow])
        newindex2 = np.arange(start=1, stop = len(avgm[0])+1)
        avgm.index = newindex2
        return avgm
    else:
        bestgen = result.iat[0,0]
        newrow = pd.DataFrame([[bestgen, result['M'].mean()]])
        avgm = pd.concat([avgm, newrow])
        newdf = df.drop(bestgen,axis=1)
        return geNorm(newdf, avgm)

def pairwiseV(datarefgenes):
    Vs = pd.DataFrame()

    buf = geNorm(datarefgenes)
    n = len(buf[0])
    m = np.arange(start=0, stop=n-2)

    a = 2
    for i in m:
        mas = np.arange(start=1, stop=a+1)
        genes2 = []
        genes3 = []
        for z in mas:
            gen = buf.iloc[n-z,0]
            genes2.append(gen)
            genes3.append(gen)
        lmas = len(mas)+1
        genes3.append(buf.iloc[n-lmas,0])

        df2 = datarefgenes.loc[:,genes2]
        df3 = datarefgenes.loc[:,genes3]

        medias2 = []
        medias3 = []
        for j in df2.index:
            media2 = gmean(df2.loc[j])
            medias2.append(media2)
            media3 = gmean(df3.loc[j])
            medias3.append(media3)

        logmedias = []
        for k,l in zip(medias2,medias3):
            n1 = k/l
            An_n1 = np.log2(n1)
            logmedias.append(An_n1)

        Vn_n1 = np.std(logmedias, ddof=1)

        newvn = pd.DataFrame([[f'V{a}/V{a+1}', Vn_n1]])
        Vs = pd.concat([Vs, newvn])

        newindex = np.arange(start=1, stop = len(Vs[0])+1)
        Vs.index = newindex

        a = a+1
    return Vs

def ploteme(eme, args):
    plt.figure()
    plt.plot(eme['Genes'], eme['M'], 'o')
    plt.xticks(rotation=45)
    plt.xlabel('refgenes')
    plt.ylabel('measured M')
    plt.xticks(rotation=45)
    plt.title('measure M')
    plt.savefig(str(args.outputfolder) + '/images/eme.png')

def plotavgm(genorm, args):
    plt.figure()
    plt.plot(genorm[0], genorm[1], 'o')
    plt.xticks(rotation=45)
    plt.xlabel('refgenes')
    plt.ylabel('Avg. M')
    plt.xticks(rotation=45)
    plt.title('Genorm result')
    plt.savefig(str(args.outputfolder) + '/images/avgm.png')

def plotuve(uve, args):
    plt.figure()
    plt.bar(uve[0], uve[1])
    plt.xticks(rotation=45)
    plt.xlabel('gen pairs')
    plt.ylabel('pairwise variation')
    plt.xticks(rotation=45)
    plt.title('Pairwise variation')
    plt.savefig(str(args.outputfolder) + '/images/uve.png')

def getnrefgenes(uvedf):
    minuve = uvedf[1].min()
    iminuve = uvedf[uvedf[1]==minuve].index
    iminuve = iminuve[0] +2
    return iminuve

def getnamesrefgenes(uvedf, genorm, args):
    if args.nrefgenes != None:
        n = args.nrefgenes
    elif args.nrefgenes == None:
        n = getnrefgenes(uvedf)

    l = np.arange(start=1, stop=n+1)

    names = []
    for i in l:
        a = genorm.loc[i, 0]
        names.append(a)
    return names

def takerefgenes(names, args):
    datarefgenes = pd.read_csv(str(args.outputfolder) + '/refgenes.csv')
    datarefgenes = datarefgenes.rename(columns={'Unnamed: 0': 'Name'})
    datarefgenes = datarefgenes.set_index('Name')
    bestrefgenes = datarefgenes[names]
    bestrefgenes.to_csv(str(args.outputfolder) + '/bestrefgenes.csv', index=True)

    return bestrefgenes

def rankfeaturegenes(data, targets, args, verbose=0):
    '''
    'data' can be refgenes (usual, fast exploration) or all genes (very large analysis, several hours) using all genes to further visualization"
    'targets' must be single column sample-class association
    'num_neighbors' can simplify analysis, default 5, shouldnt be lower than 3
    '''
    num_neighbors_neighbors = args.featureselectionneighbors
    knn = KNeighborsClassifier(n_neighbors=num_neighbors_neighbors)
    targets.set_index('SAMPLE', inplace=True)

    stargets = set(targets.index)
    sdata = set(data.index)

    if len(stargets) != len(sdata):
        if len(stargets) > len(sdata):
            notboth = stargets - sdata
            targets.drop(notboth, axis=0, inplace=True)
        elif len(sdata) > len(stargets):
            notboth = sdata - stargets
            data.drop(notboth, axis=0, inplace=True)

    targets = targets.loc[data.index,:]
    data.index = targets.index

    X = data
    y = targets

    sfs1 = SFS(knn, k_features=1, forward=False, floating=False, verbose=0, scoring='balanced_accuracy', cv=2, n_jobs=-1)
    sfs1 = sfs1.fit(X, y.values.ravel())
    subsets = sfs1.subsets_

    subsets = pd.DataFrame.from_dict(subsets)

    feature_names = subsets.iloc[3,:]

    metrics = pd.DataFrame.from_dict(sfs1.get_metric_dict()).T
    return metrics

def rankstatsrefgenes(metrics, reskrus, reswilcopairs):
     '''Takes info from kruskal and wilcoxon to create a summary dataframe of confidence for refgenes'''
     ranking = pd.DataFrame()
     ranking['Genes'] = list(reskrus.columns)
     ranking.set_index('Genes', inplace=True)
     ranking['Kruskal p-value'] = reskrus.loc['pvalue',:]
     ranking.sort_values('Kruskal p-value', ascending=True, inplace=True)

     for i,j in reswilcopairs.items():
         ranking[i] = j.loc['pvalue',:]

     return ranking

def getallgenesdf(args):
     df = pd.read_csv(str(args.outputfolder) + '/tnormcounts.csv')
     df.drop(['CodeClass', 'Accession'], axis=1, inplace=True)
     df.set_index('Name', inplace=True)
     df = df.T
     return df

def gettopngenesdf(args):
    df = pd.read_csv(str(args.outputfolder) + '/tnormcounts.csv')
    df.drop(['CodeClass', 'Accession'], axis=1, inplace=True)
    df.set_index('Name', inplace=True)
    df['mean'] = df.mean(axis=1)
    df.sort_values(by = ['mean'], ascending=False, inplace=True)
    df = df[:args.topngenestocontnorm]
    df.drop('mean', inplace=True, axis=1)
    df = df.T
    return df

def getnormfactor(refgenesdf, eme, args):

    infolanes = pd.read_csv(str(args.outputfolder) + '/infolanes.csv')
    geomeans1 = {}
    refgenesdf.set_index(infolanes['ID'], inplace=True)

    if args.contnorm == 'ponderaterefgenes':
        eme.set_index('Genes', drop=True, inplace=True)
        eme2 = eme

        eme2['M'] = eme2['M'] * len(eme['M']) / sum(eme['M'])

        for i in refgenesdf.columns:
            refgenesdf[i] = refgenesdf[i] * float(eme2.loc[i])

    for i in refgenesdf.index:
        igeomean = gmean(refgenesdf.loc[i])
        geomeans1[i] = igeomean
    prenormfactor = np.mean(list((geomeans1.values())))
    normfactor = {}
    for i,j in geomeans1.items():
        nfactor = prenormfactor/j
        normfactor[i] = nfactor

    return normfactor

def refnorm(normfactor, args):

    df = pd.read_csv(str(args.outputfolder) + '/tnormcounts.csv')

    rnormgenes = pd.DataFrame()
    rnormgenes['Name'] = df.loc[:,'Name']
    thisnormgenes = df.drop(['CodeClass','Name','Accession'], axis=1)
    for i in thisnormgenes:
        nfactor = normfactor[i]
        thisnormed = []
        for j in thisnormgenes[i]:
            thiscell = nfactor*j
            thisnormed.append(thiscell)
        rnormgenes[i] = thisnormed
    rnormgenes.set_index('Name', inplace=True)
    return rnormgenes

def grouprnormgenes(args, *dfs):
    if args.groupsinrnormgenes == 'yes':
        count = 0
        ddf2 = []
        for i in dfs:
            listkeys = list(ddf.keys())
            group = listkeys[count]
            i['group'] = group
            count += 1
            ddf2.append(i)
    else:
        ddf2 = []
        for i in dfs:
            ddf2.append(i)
    return ddf2

def pathoutrnormgenes(df, args):
    pathout = str(args.outputfolder)
    pathrnorm = pathout + '/rnormcounts.csv'
    df.to_csv(pathrnorm)
    path2rnorm = pathout + '/rnormcounts2.csv'
    df2 = df
    df2.to_csv(path2rnorm, index=False, header=False)

def pathoutadnormgenes(df, args):
    pathout = str(args.outputfolder)
    pathadnorm = pathout + '/adnormcounts.csv'
    df.to_csv(pathadnorm)
    path2adnorm = pathout + '/adnormcounts2.csv'
    df2 = df
    df2.to_csv(path2adnorm, index=False, header=False)

def adnormalization(df, args, rnormgenes):

    if args.adnormalization == 'no':
        pathoutrnormgenes(df, args)
        return df
        df.set_index('Name', inplace=True)

    df = df.T

    if args.adnormalization == 'standarization':
        if args.groupsinrnormgenes == 'yes':
            df.drop('group', axis=1, inplace=True)
        scaler = StandardScaler()
        stdnormgenes = scaler.fit_transform(df)
        stdnormgenes = pd.DataFrame(stdnormgenes)
        stdnormgenes.index = rnormgenes.columns
        stdnormgenes.columns = rnormgenes.index
        stdnormgenes = stdnormgenes.T
        pathoutadnormgenes(stdnormgenes, args)
        return stdnormgenes

    elif args.adnormalization == 'quantile':
        if args.groupsinrnormgenes == 'yes':
            df.drop('group', axis=1, inplace=True)
        df.applymap(lambda x: float(x))

        df = df.T
        qnormgenes = quantile_transform(df)
        qnormgenes = pd.DataFrame(qnormgenes)

        qnormgenes.index = df.index
        qnormgenes.columns = df.columns
        qnormgenes = qnormgenes.T
        pathoutadnormgenes(qnormgenes, args)

        return qnormgenes

def logarizeoutput(rnormgenes, args):
    if args.logarizedoutput != 'no':
        if 'group' in rnormgenes.index:
            rnormgenes2 = rnormgenes.drop('group', axis=0, inplace=False)
        else:
            rnormgenes2 = rnormgenes

        if args.logarizedoutput == '2':
            logarizedgenes = rnormgenes2.applymap(lambda x: np.log2(x))

        if args.logarizedoutput == '10':
            logarizedgenes = rnormgenes2.applymap(lambda x: np.log10(x))

        pathout = str(args.outputfolder)
        pathlogarized = pathout + '/logarized_rnormcounts.csv'
        logarizedgenes.to_csv(pathlogarized)
        return logarizedgenes

def logarizegroupedcounts(rnormgenesgroups, args):
    if args.logarizedoutput != 'no':
        rngg = rnormgenesgroups.drop('group', axis=0)
        if args.logarizeforeval == '2':
            rngg = rngg.applymap(lambda x: np.log2(x))
        if args.logarizeforeval == '10':
            rngg = rngg.applymap(lambda x: np.log10(x))

        rngg.loc['group'] = rnormgenesgroups.loc['group']
        pathout = str(args.outputfolder)
        pathlogarizedgrouped = pathout + '/logarized_grouped_rnormcounts.csv'
        rngg.to_csv(pathlogarizedgrouped)

        return rngg

    if args.logarizedoutput == 'no':
        return rnormgenesgroups

def RLEcal(rnormgenes, args):
    if args.groupsinrnormgenes == 'yes' and 'group' in rnormgenes.index:
            rlegenes = rnormgenes.drop('group', axis=0)
    else:
        rlegenes = rnormgenes

    rlegenes = rlegenes.T

    for i in rlegenes.columns:
        median = np.median(rlegenes[i])
        rlegenes[i] = rlegenes[i] - median
    rlegenes = rlegenes.T
    if 'group' in rnormgenes.index:
        rlegenes.loc['group'] = rnormgenes.loc['group']
    return rlegenes

def getmeaniqr(rlegenes):
    iqrlist = []
    iqrlist2 = []
    if 'group' in rlegenes.index:
        rlegenesng = rlegenes.drop('group', axis=0)
    else:
        rlegenesng = rlegenes
    for i in rlegenesng:
        iqr = stats.iqr(rlegenesng[i])
        iqrlist.append(iqr)
    for i in rlegenesng:
        iqr2 = stats.iqr(rlegenesng[i], rng=(10,90))
        iqrlist2.append(iqr2)
    meaniqr1 = np.mean(iqrlist)
    meaniqr2 = np.mean(iqrlist2)
    meaniqr = np.mean([meaniqr1,meaniqr2])

    meaniqr = meaniqr * 100
    return meaniqr

def logarize(x):
    x = float(x)
    x  = np.log10(x)
    return x

def plotevalnorm(matrix, what, meaniqr, args):
    matrix = matrix.apply(lambda x: x+1)
    matrix = matrix.applymap(lambda x: np.log10(x))

    matrix = matrix.T
    for i in matrix.columns:
        median = np.median(matrix[i])
        matrix[i] = matrix[i] - median

    matrix = matrix.T

    estoo = matrix

    plt.figure(figsize=(30,12))
    sns.boxplot(data=matrix, showfliers=False, showmeans=True)
    plt.title(what + '. IQR: ' + str(meaniqr), fontsize=24)
    plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
    plt.ylim(-1,1)
    plt.ylabel('RLE', fontsize=24)
    plt.xlabel('Samples', fontsize=24)
    sns.stripplot(data=matrix, size=2, palette='dark:black')
    plt.savefig(str(args.outputfolder) + '/images/rlenormplot.png')
    plt.savefig(str(args.outputfolder) + '/images/rlenormplot2.png', dpi=15)


    return estoo

def plotevalraw(matrix, what, meaniqrraw, args):
    matrix = matrix.apply(lambda x: x+1)
    matrix = matrix.applymap(lambda x: np.log10(x))

    matrix = matrix.T
    for i in matrix.columns:
        median = np.median(matrix[i])
        matrix[i] = matrix[i] - median

    matrix = matrix.T

    estoo = matrix

    plt.figure(figsize=(30,12))
    sns.boxplot(data=matrix, showfliers=False, showmeans=True)
    plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
    plt.title(what + '. IQR: ' + str(meaniqrraw), fontsize=24)
    plt.ylabel('RLE', fontsize=24)
    plt.xlabel('Samples', fontsize=24)
    plt.ylim(-1,1)
    sns.stripplot(data=matrix, size=2, palette='dark:black')
    plt.savefig(str(args.outputfolder) + '/images/rlerawplot.png')
    plt.savefig(str(args.outputfolder) + '/images/rlerawplot2.png', dpi=15)
    return estoo


def argParser():
    parser = argparse.ArgumentParser(description="Nanostring quality control analysis")
    parser.add_argument('-f', '--folder', type=str, default= pathlib.Path.cwd() / '../examples/d1_COV_GSE183071', help='relative folder where RCC set is located. Default: /data')
    parser.add_argument('-minf', '--minfov', type=float, default=0.75, help='set manually min fov for QC')
    parser.add_argument('-maxf', '--maxfov', type=float, default=1, help='set manually max fov for QC')
    parser.add_argument('-minbd', '--minbd', type=float, default=0.1, help='set manually min binding density for QC')
    parser.add_argument('-maxbd', '--maxbd', type=float, default=1.8, help='set manually max binding density for QC')
    parser.add_argument('-minlin', '--minlin', type=float, default=0.75, help='set manually min linearity for QC')
    parser.add_argument('-maxlin', '--maxlin', type=float, default=1, help='set manually max linearity for QC')
    parser.add_argument('-minscaf', '--minscalingfactor', type=float, default=0.3, help='set manually min scaling factor for QC')
    parser.add_argument('-maxscaf', '--maxscalingfactor', type=float, default=3, help='set manually max scaling factor for QC')
    parser.add_argument('-swbrrq', '--showbrowserrawqc', type=bool, default=False, help='pops up infolanes and qc summary')
    parser.add_argument('-swbrq', '--showbrowserqc', type=bool, default=False, help='pops up infolanes and qc summary')
    parser.add_argument('-swbrcn', '--showbrowsercnorm', type=bool, default=False, help='pops up infolanes and qc summary')
    parser.add_argument('-lc', '--lowcounts', type=str, default='sustract', choices=['skip', 'asim', 'sustract'],  help='what to do with counts below background?')
    parser.add_argument('-mi', '--modeid', type=str, default='filename', choices=['sampleID','filename', 'id+filename'], help='choose sample identifier. sampleID: optimal if assigned in rccs. filenames: easier to be unique. id+filename: care with group assignment coherence')
    parser.add_argument('-mv', '--modeview', type=str, default='view', choices=['justrun', 'view'], help='choose if plot graphs or just run calculations')
    parser.add_argument('-tnm', '--tecnormeth', type=str, default='regression', choices=['posgeomean','Sum', 'Median', 'regression'], help='choose method for technical normalization')
    parser.add_argument('-reg', '--refendgenes', type=str, default= 'endhkes', choices=['hkes', 'endhkes'], help='choose refgenes, housekeeping, or hkes and endogenous')
    parser.add_argument('-re', '--remove', type=str, nargs='+', default=None, help='lanes to be removed from the analysis')
    parser.add_argument('-bg', '--background', type=str, default= 'Background', choices=['Background', 'Background2', 'Background3', 'Backgroundalt'], help='choose background: b1=meancneg+(2*std), b2=maxcneg, b3=meancneg, balt=')
    parser.add_argument('-pbb', '--pbelowbackground', type=int, default=85, help='if more than %bb genes are below background, sample gets removed from analysis')
    parser.add_argument('-mbg', '--manualbackground', type=float, default=None, help='set manually background')
    parser.add_argument('-crg', '--chooserefgenes', type=str, nargs='+', default = None, help = 'list of strings like. choose manualy reference genes to use over decided-by-program ones')
    parser.add_argument('-fgv', '--filtergroupvariation', type=str, default='filterkrus', choices=['filterkrus', 'filterwilcox', 'flagkrus', 'flagwilcox', 'nofilter'], help='¿filter or flag preselected ref genes by significative group-driven differences? needs groups to be declared')
    parser.add_argument('-fsn', '--featureselectionneighbors', type=float, default=4, help='number of neighbors for feature selection analysis of refgenes. recommended 3-6')
    parser.add_argument('-g', '--groups', type=str, default='yes', choices=['yes','no'], help='defining groups for kruskal/wilcox/fs analysis?')
    parser.add_argument('-ne', '--numend', type=int, default=6, help='number of endogenous tofind by ERgene to include in analysis to check viability as refgenes')
    parser.add_argument('-ar', '--autorename', type=str, default='off', choices=['on', 'off'], help='turn on when sample IDs are not unique, be careful on sample identification detail')
    parser.add_argument('-cn', '--contnorm', type=str, default='refgenes', choices=['ponderaterefgenes', 'refgenes', 'all', 'topn'])
    parser.add_argument('-an', '--adnormalization', type=str, default='no', choices=['no', 'standarization', 'quantile'], help='perform additional normalization? standarization and quantile normalization available')
    parser.add_argument('-tn', '--topngenestocontnorm', type=int, default=100, help='set n genes to compute for calculating norm factor from top n expressed endogenous genes')
    parser.add_argument('-mch', '--mincounthkes', type=int, default=80, help='set n min counts to filter hkes candidate as refgenes')
    parser.add_argument('-nrg', '--nrefgenes', type=int, default=None, help='set n refgenes to use, overwriting geNorm calculation')
    parser.add_argument('-lr', '--laneremover', type=str, default='yes', choices=['yes', 'no'], help='option to perform analysis with all lanes if set to no')
    parser.add_argument('-grn', '--groupsinrnormgenes', type=str, default='no', choices=['yes', 'no'], help='want groups to be specified in last column of rnormgenes dataframe?')
    parser.add_argument('-lo', '--logarizedoutput', type=str, default='10', choices=['2', '10', 'no'], help='want normed output to be logarized? in what logbase?')
    parser.add_argument('-le', '--logarizeforeval', type=str, default='10', choices=['2', '10', 'no'], help= 'logarithm base for RLE calculations')
    parser.add_argument('-gf', '--groupsfile', type=str, default='../examples/groups_d1_COV_GSE183071.csv', help='enter file name where groups are defined')
    parser.add_argument('-st', '--start_time', type=float, default = time.time())
    parser.add_argument('-cs', '--current_state', type=str, default='Ready')
    parser.add_argument('-ftl', '--firsttransformlowcounts', type=bool, default=True)
    parser.add_argument('-of', '--outputfolder', type=str, default= tempfile.gettempdir() + '/guanin_output')
    parser.add_argument('-sll', '--showlastlog', type=bool, default = False)
    return parser.parse_args()

#################BIG BLOCKS -- BUTTONS
def showinfolanes(args):
    '''
    Load RCCS and show infolanes
    '''

    infolanes, dfgenes, dfnegcount, dfhkecount, dfposneg = loadrccs(args)

    createoutputfolder(args)

    exportrawinfolanes(infolanes, dfnegcount, dfhkecount, dfposneg, args)

    exportrawcounts(dfgenes, args)

    summarizerawinfolanes(args)

    def condformat_infolanes(val, top, bot, colorbien = '#a3c771', colorreg = '#f0e986', colormal = '#e3689b'):
        if top >= val >= bot:
            color = colorbien
        elif top*1.15 >= val >= bot*0.85:
            color = colorreg
        elif (bot*0.85 > val) | (val > top*1.15):
            color = colormal

        return 'background-color: {}'.format(color)

    def condformat_LOD(val, colorbien = '#a3c771', colormal = '#e3689b'):
        if val == True:
            color = colormal
        elif val == False:
            color = colorbien

        return 'background-color: {}'.format(color)


    infolanes = infolanes.style.applymap(condformat_infolanes, top=args.maxfov, bot=args.minfov, subset='FOV value')
    infolanes = infolanes.applymap(condformat_infolanes, top = args.maxbd, bot=args.minbd, subset='Binding Density')
    infolanes = infolanes.applymap(condformat_LOD, subset='limit of detection')
    infolanes = infolanes.applymap(condformat_infolanes, top= args.maxlin, bot=args.minlin,  subset='R2')
    infolanes = infolanes.applymap(condformat_infolanes, top= args.pbelowbackground, bot=0,  subset='Genes below backg %')
    infolanes = infolanes.applymap(condformat_infolanes, top= args.maxscalingfactor, bot=args.minscalingfactor,  subset='scaling factor')


    infolanes.to_html(str(args.outputfolder) + '/rawinfolanes.html')

    if args.showbrowserrawqc == True:
        webbrowser.open(str(args.outputfolder) + '/rawinfolanes.html')

def plotandreport(args, whatinfolanes = 'rawinfolanes'):
    '''
    Plots QC from infolanes or raw info lanes (file)

    :param whatinfolanes: chooses file rawinfolanes.csv or infolanes.csv
    :return:
    '''
    if whatinfolanes == 'rawinfolanes':
        infolanes = pd.read_csv(str(args.outputfolder) + '/rawinfolanes.csv', index_col='ID')
    elif whatinfolanes == 'infolanes':
        infolanes = pd.read_csv(str(args.outputfolder) + '/infolanes.csv', index_col=0)

    dfnegcount = pd.read_csv(str(args.outputfolder) + '/dfnegcount.csv', index_col=0)
    dfhkecount = pd.read_csv(str(args.outputfolder) + '/dfhkecount.csv', index_col=0)

    if args.modeview != 'justrun':
        plotfovvalue(args, infolanes)
        plotbd(args, infolanes)
        plotgenbackground(args, infolanes)
        plotld(args, infolanes)
        plotocn(args, infolanes, dfnegcount)
        plotlin(args, infolanes)
        plothke(args, infolanes, dfhkecount)
        plothkel(args, infolanes, dfhkecount)
        plotsca(args, infolanes)

    args.current_state = '--> Generating pdf report'
    print(args.current_state)
    logging.info(args.current_state)
    pdfreport(args)


def runQCview(args):
    try:
        showinfolanes(args)
        state = 'All RCCs loaded succesfully'
        args.current_state = state
        print(args.current_state)
        logging.info(state)
    except Exception as e:
        state = 'Something went wrong loading files, check input folder. Error: ' + str(e)
        args.current_state = state
        logging.error(args.current_state)
        print(args.current_state)
        if args.showbrowserrawqc == True:
            webbrowser.open(str(pathlib.Path.cwd()) + '/guanin_analysis_description.log')
        return
    try:
        plotandreport(args)
        state = 'Data loaded succesfuly, preliminary analysis and plots ready to inspect'
        args.current_state = state
        print(args.current_state)
        logging.info(state)


    except Exception as e:
        state = 'Something went wrong with preliminary analysis and/or plotting. Error: ' + str(e)
        args.current_state = state
        logging.error(args.current_state)
        print(args.current_state)

    if args.showbrowserrawqc == True:
        webbrowser.open(str(pathlib.Path.cwd()) + '/guanin_analysis_description.log')


def runQCfilterpre(args):
    '''
    This func assumes you have run runQCview and selected some filtering and configuration parameters
    '''
    dfgenes = pd.read_csv(str(args.outputfolder) + '/dfgenes.csv')
    flagged = flagqc(args)

    if args.laneremover == 'yes':
        dfgenes, infolanes = removelanes(flagged, args)
        exportfilrawcounts(dfgenes, args)
    elif args.laneremover =='no':
        dfgenes.reset_index(inplace=True)
    exportdfgenes(dfgenes, args)

    infolanes = rescalingfactor23(args)
    pathoutinfolanes(infolanes, args)

    infolanes = findaltnegatives(args)
    pathoutinfolanes(infolanes, args)

    def condformat_infolanes(val, top, bot, colorbien = '#a3c771', colorreg = '#f0e986', colormal = '#e3689b'):
        if top >= val >= bot:
            color = colorbien
        elif top*1.15 >= val >= bot*0.85:
            color = colorreg
        elif (bot*0.85 > val) | (val > top*1.15):
            color = colormal

        return 'background-color: {}'.format(color)

    def condformat_LOD(val, colorbien = '#a3c771', colormal = '#e3689b'):
        if val == True:
            color = colormal
        elif val == False:
            color = colorbien

        return 'background-color: {}'.format(color)


    infolanes = infolanes.style.applymap(condformat_infolanes, top=args.maxfov, bot=args.minfov, subset='FOV value')
    infolanes = infolanes.applymap(condformat_infolanes, top = args.maxbd, bot=args.minbd, subset='Binding Density')
    infolanes = infolanes.applymap(condformat_LOD, subset='limit of detection')
    infolanes = infolanes.applymap(condformat_infolanes, top= args.maxlin, bot=args.minlin,  subset='R2')
    infolanes = infolanes.applymap(condformat_infolanes, top= args.pbelowbackground, bot=0,  subset='Genes below backg %')
    infolanes = infolanes.applymap(condformat_infolanes, top= args.maxscalingfactor, bot=args.minscalingfactor,  subset='scaling factor')

    infolanes.to_html(str(args.outputfolder) + '/infolanes.html')

    if args.showbrowserqc == True:
        webbrowser.open(str(args.outputfolder) + '/infolanes.html')

    summarizeinfolanes(args)

    return flagged

def runQCfilter(args):
    try:
        runQCfilterpre(args)
        args.current_state = 'QC filter applied, ready to perform technical normalization'
        print(args.current_state)
        logging.info(args.current_state)


    except Exception as e:
        args.current_state = 'Unknown error while QC filtering, check input data and parameters. Error: ' + str(e)
        print(args.current_state)
        logging.info(args.current_state)

    if args.showbrowserqc == True:
        webbrowser.open(str(pathlib.Path.cwd()) + '/guanin_analysis_description.log')


def technorm(args):


    if args.firsttransformlowcounts == True:
        transformlowcounts(args)
        reinfolanes(args)
        dfgenes = pd.read_csv(str(args.outputfolder) + '/dfgenes.csv')
        if args.tecnormeth != 'regression':
            normgenes = normtecnica(dfgenes, args)
        elif args.tecnormeth == 'regression':
            normgenes = regresion(dfgenes, args)

    elif args.firsttransformlowcounts == False:
        dfgenes = pd.read_csv(str(args.outputfolder) + '/dfgenes.csv')
        if args.tecnormeth != 'regression':
            normgenes = normtecnica(dfgenes, args)
        elif args.tecnormeth == 'regression':
            normgenes = regresion(dfgenes, args)
        reinfolanes(args)
        transformlowcounts(args)

    exporttnormgenes(normgenes, args)
    exportdfgenes(normgenes, args)

    args.current_state = 'Technical normalization done'
    # print(args.current_state)
    # logging.info(args.current_state)

    # except Exception as e:
    #     args.current_state = 'Failed technical normalization'
    #     print(args.current_state, e)
    #     logging.error(args.current_state)


def contnorm(args):
    logging.info('Starting content normalization')

    if os.path.isfile(args.groupsfile):
        targets = pd.read_csv(args.groupsfile)
        groups = set(targets['GROUP'])
        if len(groups)>1:
            args.groups = 'yes'

    allhkes = getallhkes(args)
    args.current_state = '--> Selecting refgenes. Elapsed %s seconds ' + str((time.time() - args.start_time))
    print(args.current_state)
    logging.info(args.current_state)
    print('Housekeeping genes present in analysis: ', list(allhkes.index))

    selhkes = filter50chkes(allhkes, args)
    if len(selhkes.index) <= 2:
        selhkes = allhkes
        args.current_state = 'All or almost all housekeeping genes are low expressed. Consider re-design experiment'
        print(args.current_state)
        logging.warning(args.current_state)
    else:
        args.current_state = 'Housekeeping genes with more than 50 counts for all lanes: ' + str(list(selhkes.index))
        print(args.current_state)
        logging.info(args.current_state)

    try:
        refgenes = findrefend(args, selhkes)
        args.current_state = 'Refgenes in analysis including housekeepings + best endogenous selected: ' +  str(list(refgenes.index))
        print(args.current_state)
        logging.info(args.current_state)
    except Exception as e:
        logging.warning('Unable to retrieve candidate ref genes from endogenous, ERROR: ', e)
        refgenes = selhkes

    pathoutrefgenes(refgenes, args)

    args.current_state = str(('--> Group-driven refining candidate reference genes selection through kruskal, wilcoxon and feature selection. Elapsed %s seconds ' % (time.time() - args.start_time)))
    print(args.current_state)
    logging.info(args.current_state)

    if args.groups == 'yes':
        args.current_state = '--> Performing kruskal-wallis analysis'
        print(args.current_state)
        logging.info(args.current_state)
        ddf = getgroups(args)
        ddfc = list(ddf.items())
        ddfb = list(ddf.values())
        reskrus = calkruskal(*ddfb)

        args.current_state = '--> Performing wilcoxon analysis'
        print(args.current_state)
        logging.info(args.current_state)
        reswilcopairs = calwilcopairs(*ddfc)

        flaggedgenes = flagkrus(reskrus)
        flaggedwilcox = flagwilcox(reswilcopairs)

        flaggedboth = set(flaggedgenes).intersection(set(flaggedwilcox))

        print('Flagged genes by kruskal-wallis and/or wilcoxon: ' + str(flaggedboth))
        logging.info('Flagged genes by kruskal-wallis and/or wilcoxon: ' + str(flaggedboth))

        if args.filtergroupvariation == 'filterkrus' or args.filtergroupvariation == 'flagkrus':
            refgenes = filterkruskal(flaggedgenes, args)
        elif args.filtergroupvariation == 'filterwilcox' or args.filtergroupvariation == 'flagwilcox':
            refgenes = filterwilcox(flaggedwilcox, args)

        print('Ref genes present in analysis after applying kruskal-wallis or wilcoxon filtering: ',
              list(refgenes.columns))

    elif args.groups == 'no':
        pass

    print(
        '--> Applying genorm to select best ranking selection of refgenes from candidate refgenes. Elapsed %s seconds ' % (
                    time.time() - args.start_time))
    datarefgenes = pd.read_csv(str(args.outputfolder) + '/refgenes.csv', index_col=0)

    eme = measureM(datarefgenes)

    genorm = geNorm(datarefgenes)

    uve = pairwiseV(datarefgenes)

    print('--> Performing geNorm calculations')

    ploteme(eme, args)

    plotavgm(genorm, args)

    plotuve(uve, args)

    args.current_state = '--> Calculating optimal number of reference genes'
    print('--> Calculating optimal number of reference genes')

    if args.chooserefgenes == None:
        names = getnamesrefgenes(uve, genorm, args)
        if args.nrefgenes == None:
            args.current_state = 'Ref. genes selected (auto): ' + str(names)
            args.refgenessel = str(names)
            print(args.current_state)
        elif args.nrefgenes != None:
            args.current_state = 'Ref. genes selected (n_manual): ' + str(names)
            args.refgenessel = str(names)
            print(args.current_state)
    else:
        if type(args.chooserefgenes) == str:
            names = args.chooserefgenes.split()
        else:
            names = args.chooserefgenes

        args.current_state = 'Ref. genes selected (manual): ' + str(names)
        args.refgenessel = str(names)
        print(args.current_state)

    bestrefgenes = takerefgenes(names, args)

    dataref = pd.read_csv(str(args.outputfolder) + '/refgenes.csv', index_col=0)


    print('--> Performing feature selection for refgenes evaluation and control.')
    print(args.groups)
    if (args.groups == 'yes') | (os.path.exists(args.groupsfile) == True):
        print('2')
        metrics = rankfeaturegenes(dataref, targets, args)
        print(metrics)
        ranking = rankstatsrefgenes(metrics, reskrus, reswilcopairs)
        print(ranking)

        def condformat_ranking(val, top, bot, colorbien='#a3c771', colorreg='#f0e986', colormal='#e3689b'):
            if top >= val >= bot:
                color = colorbien
            elif top * 1.15 >= val >= bot * 0.85:
                color = colorreg
            elif (bot * 0.85 > val) | (val > top * 1.15):
                color = colormal

            return 'background-color: {}'.format(color)

        ranking2 = ranking.style.applymap(condformat_ranking, top = 1, bot=0.05, subset='Kruskal p-value')
        for i in ranking2.columns:
            ranking2 = ranking2.applymap(condformat_ranking, top=1, bot=0.05, subset = i)

        ranking2.to_html(str(args.outputfolder) + '/ranking_kruskal_wilcox.html')
        print('done')
        ranking.to_csv(str(args.outputfolder) + '/reports/ranking_kruskal_wilcox.csv')

    if args.showbrowsercnorm == True:
            webbrowser.open(str(args.outputfolder) + '/ranking_kruskal_wilcox.html')

    if args.groups == 'yes':
        def condformat_metrics(val, top, bot, colorbien = '#a3c771', colorreg = '#f0e986', colormal = '#e3689b'):
            if top >= val >= bot:
                color = colorbien
            elif top*1.15 >= val >= bot*0.85:
                color = colorreg
            elif (bot*0.85 > val) | (val > top*1.15):
                color = colormal

            return 'background-color: {}'.format(color)

        metrics2 = metrics.style.applymap(condformat_metrics, top=1.5/len(groups), bot=0.5/len(groups), subset='avg_score')

        metrics2.to_html(str(args.outputfolder) + '/metrics_reverse_feature_selection.html')
        metrics.to_csv(str(args.outputfolder) + '/reports/metrics_reverse_feature_selection.csv')

    if (args.showbrowsercnorm == True) and (args.groups == 'yes'):
        webbrowser.open(str(args.outputfolder) + '/metrics_reverse_feature_selection.html')


    print('--> Getting lane-specific normfactor and applying content normalization. Elapsed %s seconds ' % (
                time.time() - args.start_time))

    allgenes = getallgenesdf(args)

    if args.contnorm == 'refgenes' or args.contnorm == 'ponderaterefgenes':
        normfactor = getnormfactor(bestrefgenes, eme, args)
    elif args.contnorm == 'all':
        normfactor = getnormfactor(allgenes, eme, args)
    elif args.contnorm == 'topn':
        topngenes = gettopngenesdf(args)
        normfactor = getnormfactor(topngenes, eme, args)

    rnormgenes = refnorm(normfactor, args)
    pathoutrnormgenes(rnormgenes, args)

    print('--> Performing additional normalization. Elapsed %s seconds ' % (time.time() - args.start_time))

    if args.groupsinrnormgenes == 'yes' and args.groups == 'yes':
        adnormgenes = adnormalization(rnormgenesgroups, args, rnormgenes)
    elif args.groupsinrnormgenes == 'no' or args.groups == 'no':
        adnormgenes = adnormalization(rnormgenes, args, rnormgenes)

    print('--> Exporting normalization results. Elapsed %s seconds ' % (time.time() - args.start_time))
    pathoutadnormgenes(adnormgenes, args)

    if args.groups == 'yes' and args.groupsinrnormgenes == 'yes':
        rngg = logarizegroupedcounts(rnormgenesgroups, args)
    elif args.groups == 'yes' and args.groupsinrnormgenes == 'no':
        rngg = logarizeoutput(rnormgenes, args)
    else:
        rngg = logarizeoutput(adnormgenes, args)

    rngg.to_csv(str(args.outputfolder) + '/rngg.csv', index=True)
    return rngg, names

def evalnorm(args):
    args.current_state = '--> Evaluating and plotting normalization results. Elapsed %s seconds ' + str((time.time() - args.start_time))
    print(args.current_state)
    logging.info(args.current_state)
    rngg = pd.read_csv(str(args.outputfolder) + '/rngg.csv', index_col = 'Name')
    rawcounts = pd.read_csv(str(args.outputfolder) + '/rawcounts2.csv', index_col=0)
    rlegenes = RLEcal(rngg, args)
    rleraw = RLEcal(logarizeoutput(rawcounts, args), args)

    meaniqr = getmeaniqr(rlegenes)
    meaniqrraw = getmeaniqr(rleraw)

    rawcounts = pd.read_csv(str(args.outputfolder) + '/rawcounts.csv', index_col='Name')
    rawcounts.drop(['CodeClass', 'Accession'], inplace=True, axis=1)

    rawfcounts = pd.read_csv(str(args.outputfolder) + '/rawfcounts.csv', index_col='Name')
    rawfcounts.drop(['CodeClass', 'Accession'], inplace=True, axis=1)

    tnormcounts = pd.read_csv(str(args.outputfolder) + '/tnormcounts.csv', index_col='Name')
    tnormcounts.drop(['CodeClass', 'Accession'], inplace=True, axis=1)

    rnormcounts = pd.read_csv(str(args.outputfolder) + '/rnormcounts.csv', index_col='Name')

    print('Plotting raw RLE plot...')
    plotevalraw(rawcounts, 'RAW counts', meaniqrraw, args)
    logging.info('Plotted raw RLE plot')

    print('Plotting normalized RLE plot...')
    plotevalnorm(rnormcounts, 'fully normalized counts', meaniqr, args)
    logging.info('Plotted normalized RLE plot')

    pdfreportnorm(args)

    args.current_state = '--> Finished. Elapsed %s seconds ' + str((time.time() - args.start_time))
    print(args.current_state)
    logging.info(args.current_state)

    if args.showlastlog == True:
        webbrowser.open(str(pathlib.Path.cwd()) + '/guanin_analysis_description.log')

    return (meaniqrraw, meaniqr)


if __name__ == '__main__':
    args = argParser()
    runQCview(args)
    runQCfilter(args)
    technorm(args)
    contnorm(args)
    evalnorm(args)