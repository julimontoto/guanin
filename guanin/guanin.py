import math
import os
from pathlib import Path
import pathlib
import statistics
import tempfile
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
from scipy.stats.mstats import gmean
from scipy import stats
from scipy.linalg import svd
# Monkeypatch matplotlib to avoid the failing from upsetplot
# from matplotlib.pyplot import tight_layout
# tight_layout.get_renderer = ""
import logging
import argparse
from fpdf import FPDF
from sklearn.decomposition import PCA
from mlxtend.feature_selection import SequentialFeatureSelector as SFS
import seaborn as sns
import time
import webbrowser
from sklearn.preprocessing import MinMaxScaler, quantile_transform, scale, PowerTransformer
from sklearn.neighbors import KNeighborsClassifier
from ERgene import FindERG
from pydeseq2.preprocessing import deseq2_norm
import traceback
import sys
import warnings


class Setup:
    def __init__(self, cmd_line_args, *args, **kwargs):
        self.set_output(cmd_line_args.output)

    def set_output(self, output):
        """Build the most probable default dir.

        - If we detect NT (Windows), defaults to User %LOCALAPPDATA% / guanin
        - If we detect other (Posix), defaults to $XDG_DATA_HOME / guanin
        - If the variables %LOCALAPPDATA% or XDG_DATA_HOME are not set, the
          base is at $TMP / userlogin / guanin
        """
        if output:
            self.output = output
        else:
            base = Path(tempfile.gettempdir()) / os.getlogin()
            if os.name == "posix":
                base = Path(os.getenv("XDG_DATA_HOME", base))
            else:
                base = Path(os.getenv("LOCALAPPDATA", base))

            self.output = base / "guanin"

        Path(self.output).mkdir(parents=True, exist_ok=True)


def default_output(make_dir=True):
    """Build the most probable default dir.

    - If we detect NT (Windows), defaults to User %LOCALAPPDATA% / guanin
    - If we detect other (Posix), defaults to $XDG_DATA_HOME / guanin
    - If the variables %LOCALAPPDATA% or XDG_DATA_HOME are not set, the base
      is at $TMP / userlogin / guanin

    """
    apps_dir = Path(tempfile.gettempdir()) / os.getlogin()
    if os.name == "nt":
        apps_dir = Path(os.getenv("LOCALAPPDATA", apps_dir))
    else:  # posix
        apps_dir = Path(os.getenv("XDG_DATA_HOME", apps_dir))
    app_dir = Path(apps_dir / "guanin")
    if make_dir:
        app_dir.mkdir(parents=True, exist_ok=True)

    return app_dir


def getfolderpath(folder):
    """RCC path."""
    cwd = Path(__file__).parent.absolute()
    path = cwd / folder
    return path


def loadrccs(args, start_time = 0):
    """RCC loading to extract information."""
    columns = ['ID', 'Comments', 'FOV value', 'Binding Density', 'Background',
               'Background2', 'Background3', 'Genes below backg %', 'nGenes',
               'posGEOMEAN', 'Sum', 'Median', 'R2', 'limit of detection',
               '0,5fm']
    infolanes = pd.DataFrame(columns = columns)

    dfgenes = pd.DataFrame() # counts
    dfposneg = {} # posnegs

    geomeans = []
    logconc = [7, 5, 3, 1, -1]

    dfnegcount = pd.DataFrame()
    negnames = []
    dfhkecount = pd.DataFrame()
    hkenames = []

    a = 0
    for file in os.listdir(getfolderpath(args.folder)):
        '''First data inspection'''
        if '.RCC' in file:
            # Count info dataframe
            df = pd.read_csv((getfolderpath(args.folder) / file),
                             names=['CodeClass', 'Name', 'Accession', 'Count'])
            df = df.dropna()

            #separate dataframe for gene class
            posvals = ["Positive", "Positive1", "Positive2"]
            dfpos = df[df["CodeClass"].isin(posvals)]
            negvals = ["Negative", "Negative1", "Negative2"]
            dfneg = df[df["CodeClass"].isin(negvals)]
            dfposneg1 = pd.concat([dfpos,dfneg])
            endvals = ["Endogenous", "Endogenous1", "Endogenous2"]
            dfend = df[df["CodeClass"].isin(endvals)]
            dfhke = df[df.CodeClass =='Housekeeping']


            convert_dict = {'CodeClass': str,
                            'Name': str,
                            'Accession': str,
                            'Count': float}
            dfpos = dfpos.astype(convert_dict)
            dfpos = dfpos.sort_values(by=['Count'], ascending=False)
            dfneg = dfneg.astype(convert_dict)
            dfend = dfend.astype(convert_dict)
            dfhke = dfhke.astype(convert_dict)
            dff = pd.concat([dfend, dfhke, dfneg, dfpos])

            names = ['parametro', 'valor']
            dinf = pd.read_csv(getfolderpath(args.folder) / file,
                               names=names,
                               nrows=30,
                               on_bad_lines='skip')  # dataframe with lane info
            dinf = dinf.dropna()
            thislane = []  # info for this lane

            # adds id from sample or from file name
            id1 = dinf.loc[dinf['parametro'] == 'ID']
            id1 = id1.loc[:, 'valor']
            id1 = id1.iloc[0]
            id1 = str(id1)

            if args.modeid == 'filename':
                id1 = str(file)
            elif args.modeid == 'id+filename':
                id1 += str(file)

            args.current_state = f"Loading RCC... {id1}"
            logging.info(args.current_state)

            if args.autorename:
                id1 = f"{a} Renamed".strip()

            comments = dinf.loc[dinf['parametro'] == 'Comments']
            comments = comments.loc[:, 'valor']
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
            fovvalue = float(fovcounted.iloc[0])/float(fovcount.iloc[0])

            thislane.append(fovvalue)

            bd = dinf.loc[dinf['parametro'] == 'BindingDensity']
            bd = bd.loc[:, 'valor']
            bd = float(bd.iloc[0])
            thislane.append(bd)

            topneg = 3*(np.mean(dfneg['Count']))

            dfnegin = dfneg[dfneg.Count <= topneg]

            background = (np.mean(dfnegin.Count)) + 2 * (np.std(dfnegin.Count))
            thislane.append(background)

            background2 = np.max(dfnegin['Count'])
            thislane.append(background2)

            background3 = np.mean(dfnegin['Count'])
            thislane.append(background3)

            negs = list(dfneg['Count'])
            negnames = dfneg['Name']
            negnames = list(negnames)

            maxout = 3 * background3
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

            thislane.append(gbb * 100 / ngen)
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
                if diff is not None:
                    diff.append(list(set(dfgenes.index) - set(dff.index)))
                common = list(set(list(dff.index)).intersection(list(dfgenes.index)))
                dff = dff.loc[common]
                dfgenes = dfgenes.loc[common]
                if diff is not None:
                    logging.warning(
                        'Mismatch, genes not present in all samples: {diff}')

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

            thislogsconc = []  # R2 calculation
            thiscountnolow = []

            # excluding lowest negative control
            y = 0
            while y < 5:
                thiscountnolow.append(dfpos['Count'].iloc[y])
                y = y + 1

            # log2 counts
            for i in thiscountnolow:
                thislog = math.log(i, 2)
                thislogsconc.append(thislog)

            # r2 score calculation
            R2 = np.corrcoef(thislogsconc, logconc)
            R2 = R2[1, 0]

            thislane.append(R2)

            if args.background == 'Manual':
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
            a = a + 1

    # CHECK FOR DUPLICATED IDS
    if not all(infolanes.duplicated(subset=['ID'])):
        args.current_state = f"--> All {len(infolanes['ID'])} IDs are unique, "
        args.current_state += "proceeding with analysis. "
        logging.info(args.current_state)
    elif any(infolanes.duplicated(subset=['ID'])):
        args.current_state = \
            "--> Duplicated IDs, rename samples with unique names " +\
            "or turn on autorename option"
        logging.warning(args.current_state)


    # Adding calculated params to infolanes
    meangeomeans = np.mean(infolanes['posGEOMEAN'])
    scalingf = []

    manualbglist = []
    for i in infolanes['ID']:
        manualbglist.append(args.manualbackground)
    infolanes['manual background'] = manualbglist

    for i in infolanes['posGEOMEAN']:
        scaling = meangeomeans / i
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
    args.outputfolder = Path(args.outputfolder)
    pathout = args.outputfolder
    pathout.mkdir(parents=True, exist_ok=True)
    pathoutimages = args.outputfolder / 'images'
    pathoutimages.mkdir(parents=True, exist_ok=True)
    pathoutotherfiles = args.outputfolder / 'otherfiles'
    pathoutotherfiles.mkdir(parents=True, exist_ok=True)
    pathoutinfo = args.outputfolder / 'info'
    pathoutinfo.mkdir(parents=True, exist_ok=True)
    pathoutresults = args.outputfolder / 'results'
    pathoutresults.mkdir(parents=True, exist_ok=True)
    pathoutreports = args.outputfolder / 'reports'
    pathoutreports.mkdir(parents=True, exist_ok=True)


def exportrawcounts(rawcounts, args):
    pathraw = args.outputfolder / 'otherfiles' / 'rawcounts.csv'
    rawcounts.to_csv(pathraw, index=True)
    pathdfraw = args.outputfolder / 'otherfiles' / 'dfgenes_raw.csv'
    rawcounts.drop(['CodeClass', 'Accession'], axis=1)
    rawcounts.to_csv(pathdfraw, index=True)
    pathdfgenes = args.outputfolder / 'otherfiles' / 'dfgenes.csv'
    rawcounts.to_csv(pathdfgenes, index=True)

    rawcounts3 = rawcounts.drop(['CodeClass', 'Accession'], axis=1)
    pathraw3 = args.outputfolder / 'otherfiles' / 'rawcounts2.csv'
    rawcounts3.to_csv(pathraw3, index=True)

def exportdfgenes(dfgenes, args):
    pathdfgenes = args.outputfolder / 'otherfiles' /'dfgenes.csv'
    dfgenes.to_csv(pathdfgenes, index=True)

def exportdfgenes_qc(dfgenes, args):
    pathdfqc = args.outputfolder / 'otherfiles' / 'dfgenes_qc.csv'
    dfgenes.to_csv(pathdfqc, index=True)

def exportrawinfolanes(infolanes, dfnegcount, dfhkecount, dfposneg, args):
    """Exports raw infolanes and dfnegcount, dfhkecount, dfposneg"""
    pathinfolanes = args.outputfolder / 'info' / 'rawinfolanes.csv'
    pathdfnegcount = args.outputfolder / 'otherfiles' / 'dfnegcount.csv'
    pathdfhkecount = args.outputfolder / 'otherfiles' / 'dfhkecount.csv'

    exportposneg(dfposneg, args)
    infolanes.to_csv(pathinfolanes)
    dfnegcount.to_csv(pathdfnegcount)
    dfhkecount.to_csv(pathdfhkecount)

def pathoutinfolanes(infolanes, args):
    pathinfolanes = args.outputfolder / 'info' / 'infolanes.csv'
    infolanes.to_csv(pathinfolanes, index=True)

def pathoutrawsummary(rawsummary, args):
    pathrawsummary = args.outputfolder / 'info' / 'rawsummary.csv'
    rawsummary.to_csv(pathrawsummary, index=True)


def pathoutsummary (summary, args):
    pathsummary = args.outputfolder / 'info' / 'summary.csv'
    summary.to_csv(pathsummary, index=True)

def format_dfgenesall_to_namematrix(args, dfgenes):
    dfgenes = dfgenes.drop(['CodeClass', 'Accession'], axis=1)
    dfgenes.set_index('Name', drop=True, inplace=True)
    return dfgenes

def format_namematrix_to_dfgenes(args, namematrix):
    dfwithgeneinfo = pd.read_csv(args.outputfolder / 'otherfiles' / 'dfgenes_raw.csv', index_col='Name', usecols=['CodeClass', 'Accession'])
    namematrix.merge(dfwithgeneinfo, on='Name', how='left')
    return namematrix


def condformat(val,
                       top,
                       bot,
                       colorbien='#a3c771',
                       colorreg='#f0e986',
                       colormal='#e3689b'):
    if top >= val >= bot:
        color = colorbien
    elif top * 1.15 >= val >= bot * 0.85:
        color = colorreg
    elif (bot * 0.85 > val) | (val > top * 1.15):
        color = colormal

    return 'background-color: {}'.format(color)

def condformat_LOD(val, colorbien='#a3c771', colormal='#e3689b'):
    color = colormal if val else colorbien

    return f"background-color: {color}"


def summarizerawinfolanes(args):

    rawinfolanes = pd.read_csv(
        args.outputfolder / 'info' / 'rawinfolanes.csv',
        index_col='ID')

    rawinfofov = [np.min(rawinfolanes['FOV value']),
                  np.max(rawinfolanes['FOV value']),
                  np.mean(rawinfolanes['FOV value']),
                  np.median(rawinfolanes['FOV value'])]
    rawinfobd = [np.min(rawinfolanes['Binding Density']),
                 np.max(rawinfolanes['Binding Density']),
                 np.mean(rawinfolanes['Binding Density']),
                 np.median(rawinfolanes['Binding Density'])]
    rawinfolin = [np.min(rawinfolanes['R2']),
                  np.max(rawinfolanes['R2']),
                  np.mean(rawinfolanes['R2']),
                  np.median(rawinfolanes['R2'])]
    rawinfobackg = [np.min(rawinfolanes['Background']),
                    np.max(rawinfolanes['Background']),
                    np.mean(rawinfolanes['Background']),
                    np.median(rawinfolanes['Background'])]
    rawinfogbb = [np.min(rawinfolanes['Genes below backg %']),
                  np.max(rawinfolanes['Genes below backg %']),
                  np.mean(rawinfolanes['Genes below backg %']),
                  np.median(rawinfolanes['Genes below backg %'])]
    rawinfopgm = [np.min(rawinfolanes['posGEOMEAN']),
                  np.max(rawinfolanes['posGEOMEAN']),
                  np.mean(rawinfolanes['posGEOMEAN']),
                  np.median(rawinfolanes['posGEOMEAN'])]
    rawinfosum = [np.min(rawinfolanes['Sum']),
                  np.max(rawinfolanes['Sum']),
                  np.mean(rawinfolanes['Sum']),
                  np.median(rawinfolanes['Sum'])]
    rawinfo05fm = [np.min(rawinfolanes['0,5fm']),
                   np.max(rawinfolanes['0,5fm']),
                   np.mean(rawinfolanes['0,5fm']),
                   np.median(rawinfolanes['0,5fm'])]
    rawinfoscaf = [np.min(rawinfolanes['scaling factor']),
                   np.max(rawinfolanes['scaling factor']),
                   np.mean(rawinfolanes['scaling factor']),
                   np.median(rawinfolanes['scaling factor'])]

    rawsummary = pd.DataFrame(columns=['min', 'max', 'mean', 'Median'])
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
    rawsummary = rawsummary.style.applymap(condformat, top=args.maxfov, bot=args.minfov, subset='FOV')
    rawsummary = rawsummary.applymap(condformat, top = args.maxbd, bot=args.minbd, subset='Binding density')
    rawsummary = rawsummary.applymap(condformat, top= args.maxlin, bot=args.minlin,  subset='R2')
    rawsummary = rawsummary.applymap(condformat, top= args.pbelowbackground, bot=0,  subset='Genes below background')
    rawsummary = rawsummary.applymap(condformat, top= args.maxscalingfactor, bot=args.minscalingfactor,  subset='Scaling factor')


    rawsummary.to_html(str(args.outputfolder / 'info' / 'rawsummary.html'))

def summarizeinfolanes(args):
    infolanes = pd.read_csv(args.outputfolder / 'info' / "infolanes.csv",
        index_col='ID')
    infofov = [
        np.min(infolanes['FOV value']),
        np.max(infolanes['FOV value']),
        np.mean(infolanes['FOV value']),
        np.median(infolanes['FOV value'])
    ]
    infobd = [
        np.min(infolanes['Binding Density']),
        np.max(infolanes['Binding Density']),
        np.mean(infolanes['Binding Density']),
        np.median(infolanes['Binding Density'])
    ]
    infolin = [
        np.min(infolanes['R2']),
        np.max(infolanes['R2']),
        np.mean(infolanes['R2']),
        np.median(infolanes['R2'])
    ]
    infobackg = [
        np.min(infolanes['Background']),
        np.max(infolanes['Background']),
        np.mean(infolanes['Background']),
        np.median(infolanes['Background'])
    ]
    infogbb = [
        np.min(infolanes['Genes below backg %']),
        np.max(infolanes['Genes below backg %']),
        np.mean(infolanes['Genes below backg %']),
        np.median(infolanes['Genes below backg %'])
    ]
    infopgm = [
        np.min(infolanes['posGEOMEAN']),
        np.max(infolanes['posGEOMEAN']),
        np.mean(infolanes['posGEOMEAN']),
        np.median(infolanes['posGEOMEAN'])
    ]
    infosum = [
        np.min(infolanes['Sum']),
        np.max(infolanes['Sum']),
        np.mean(infolanes['Sum']),
        np.median(infolanes['Sum'])
    ]
    info05fm = [
        np.min(infolanes['0,5fm']),
        np.max(infolanes['0,5fm']),
        np.mean(infolanes['0,5fm']),
        np.median(infolanes['0,5fm'])
    ]
    infoscaf = [
        np.min(infolanes['scaling factor']),
        np.max(infolanes['scaling factor']),
        np.mean(infolanes['scaling factor']),
        np.median(infolanes['scaling factor'])
    ]

    summary = pd.DataFrame(columns=['min', 'max', 'mean', 'Median'])
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
    summary2view = summary.style.applymap(
        condformat,
        top=args.maxfov,
        bot=args.minfov,
        subset='FOV')
    summary2view = summary2view.applymap(
        condformat,
        top=args.maxbd,
        bot=args.minbd,
        subset='Binding density')
    summary2view = summary2view.applymap(
        condformat,
        top=args.maxlin,
        bot=args.minlin,
        subset='R2')
    summary2view = summary2view.applymap(
        condformat,
        top=args.pbelowbackground,
        bot=0,
        subset='Genes below background')
    summary2view = summary2view.applymap(
        condformat,
        top=args.maxscalingfactor,
        bot=args.minscalingfactor,
        subset='Scaling factor')

    pathoutsummary(summary, args)
    summary2view.to_html(str(args.outputfolder / "info" / "Summary.html"))

    if args.showbrowserqc:
        webbrowser.open(str(args.outputfolder / "info" / "Summary.html"))

def exportposneg(dfposneg, args):
    posnegcounts = pd.DataFrame()

    for i in dfposneg.keys():
        a = dfposneg[i]
        posnegcounts['Name'] = a['Name']
        posnegcounts['CodeClass'] = a['CodeClass']
        posnegcounts[i] = a['Count']

    posnegcounts.set_index('Name', drop=True, inplace=True)

    pathposneg = args.outputfolder / 'otherfiles' / 'posnegcounts.csv'
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
    plt.savefig(str(args.outputfolder / 'images' / 'fovplot.png'))
    plt.close()


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
    plt.savefig(str(args.outputfolder / 'images' / 'bdplot.png'))
    plt.close()


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
    plt.savefig(str(args.outputfolder / 'images' /'genbackground.png'))
    plt.close()


def plotld(args, infolanes):
    plt.figure()
    plt.plot(infolanes.index, infolanes['0,5fm'], 'bo', label='0,5fm')

    if args.background == 'Manual':
        background = 'manual background'
    else:
        background = args.background
        if background == 'Backgroundalt':
            background = 'Background'

    plt.plot(infolanes.index, infolanes[background], 'r', label='Background')
    plt.tick_params(
        axis='x',
        which='both',
        bottom=False,
        top=False,
        labelbottom=False)
    plt.title('Limit of detection')
    plt.xlabel('samples')
    plt.ylabel('0,5 fm')
    plt.legend()
    plt.grid(True)
    plt.savefig(str(args.outputfolder / 'images' / 'ldplot.png'))
    plt.close()


def plotocn(args, infolanes, dfnegcount):
    plt.figure()

    if args.background == 'Manual':
        background = 'manual background'
    else:
        background = args.background
        if background == 'Backgroundalt':
            background = 'Background'

    for i in dfnegcount.columns:
        plt.plot(dfnegcount.index, dfnegcount[i], 'o', label=i,)
    plt.tick_params(
        axis='x',
        which='both',
        bottom=False,
        top=False,
        labelbottom=False)
    plt.plot(dfnegcount.index, dfnegcount['maxoutlier'])
    plt.plot(infolanes[background], label='Background')
    plt.xlabel('samples')
    plt.ylabel('counts')
    plt.legend(loc='lower right', ncol=4, mode='expand')
    plt.title('Outliers in neg_controls')
    plt.savefig(str(args.outputfolder / 'images' / 'ocnplot.png'))
    plt.close()


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
    plt.bar(
        infolanes.index,
        ((infolanes['nGenes'] - infolanes['Genes below backg %']) /
         infolanes['nGenes']),
        color='cyan')
    plt.plot(infolanes.index, infolanes['R2'], 'o', color='blue')
    plt.plot(infolanes.index, minlin, 'm')
    plt.plot(infolanes.index, optlin, 'g')
    plt.tick_params(
        axis='x',
        which='both',
        bottom=False,
        top=False,
        labelbottom=False)
    plt.xlabel('samples')
    cyanbar = mpatches.Patch(color='cyan', label='% genes > background')
    purpleline = mpatches.Patch(color='m', label='min value')
    bluedot = mpatches.Patch(color='blue', label='R2 value')
    plt.legend(
        handles=[cyanbar, purpleline, bluedot],
        loc='lower right',
        ncol=3,
        mode='expand')
    plt.ylabel('genes')
    plt.xticks(rotation=45)
    plt.title('Linearity and genes above background')
    plt.savefig(str(args.outputfolder / 'images' / 'linplot.png'))
    plt.close()

def plothke(args, infolanes, dfhkecount):
    """Housekeeping plot."""
    plt.figure()
    for i in dfhkecount.columns:
        plt.plot(dfhkecount.index, dfhkecount[i], 'o', label=i)
    plt.xlabel('samples')
    plt.ylabel('counts')
    plt.legend(loc='upper left', ncol=3, mode='expand')
    plt.tick_params(
        axis='x',
        which='both',
        bottom=False,
        top=False,
        labelbottom=False)
    plt.title('Housekeeping genes')
    plt.savefig(str(args.outputfolder / 'images' / 'hkeplot.png'))
    plt.close()


def plothkel(args, infolanes, dfhkecount):
    """Closest housekeeping to background plot."""
    bb = np.mean(infolanes['Background'])

    bbmax = 6 * bb
    bblist = []

    for i in infolanes.index:
        bblist.append(bb)

    hkelplot = plt.figure()
    ax1 = hkelplot.add_subplot(111)
    number_of_plots = len(dfhkecount.columns)
    colors = sns.color_palette("hls", number_of_plots)
    ax1.set_prop_cycle('color', colors)
    for i in dfhkecount.columns:
        ax1.plot(dfhkecount.index, dfhkecount[i], 'o', label=i)
    plt.plot(infolanes.index, bblist, 'r')
    plt.xlabel('ID')
    plt.ylabel('counts')
    plt.ylim(0, 2 * bbmax)
    plt.legend(loc='upper left', ncol=3, mode='expand')
    plt.tick_params(
        axis='x',
        which='both',
        bottom=False,
        top=False,
        labelbottom=False)
    plt.title('Housekeeping genes close to background')
    plt.savefig(str(args.outputfolder / 'images' / 'hkelplot.png'))
    plt.close()


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
    plt.tick_params(
        axis='x',
        which='both',
        bottom=False,
        top=False,
        labelbottom=False)
    redline = mpatches.Patch(color='red', label='max scaling factor')
    purpleline = mpatches.Patch(color='m', label='min scaling factor')
    bluedot = mpatches.Patch(color='blue', label='sample scaling factor')
    plt.legend(
        handles=[redline, purpleline, bluedot],
        loc='lower right',
        ncol=2,
        mode='expand')
    plt.xlabel('samples')
    plt.ylabel('scaling factor')
    plt.title('scaling factor')
    plt.savefig(args.outputfolder / 'images' / 'scaplot.png')
    plt.close()

def pdfreport(args, rawreport=True):
    pdf = FPDF()
    pdf.add_page()
    pdf.set_font('Arial', 'B', 16)

    pdf.image(
        str(Path(__file__).parent
            / "reports" / "images" / "qc_template_report.png"),
        0,
        0,
        h=297)

    pdf.image(str(args.outputfolder / "images" / "ldplot.png"), 12.5, 42, h=69)
    pdf.image(str(args.outputfolder / "images" / "bdplot.png"), 110, 42, h=69)

    pdf.image(str(args.outputfolder / "images" / "fovplot.png"),
              10.5,
              120,
              h=69)
    pdf.image(str(args.outputfolder / "images" / "linplot.png"),
              10.5,
              200,
              h=69)
    pdf.image(str(args.outputfolder / "images" / "hkelplot.png"),
              110,
              120,
              h=69)
    pdf.image(str(args.outputfolder / "images" / "ocnplot.png"),
              110,
              200,
              h=69)
    if rawreport:
        pdf.output(str(args.outputfolder / "reports" / "QC_inspection.pdf"), 'F')
    elif not rawreport:
        pdf.output(str(args.outputfolder / "reports" / "QC_inspection_filtered.pdf"), 'F')

    if args.showbrowserqc and rawreport:
        os.system(str(args.outputfolder / "reports" / "QC_inspection.pdf"))
    elif args.showbrowserqc and not rawreport:
        os.system(str(args.outputfolder / "reports" / "QC_inspection_filtered.pdf"))

def pdfreportnorm(args):
    pdf = FPDF()
    pdf.add_page()
    pdf.set_font('Arial', 'B', 16)

    pdf.image(
        str(Path(__file__).parent /
            "reports" / "images" / "norm_template_report.png"),
        0,
        0,
        h=297)

    pdf.image(str(args.outputfolder / "images" / "avgm.png"), 12.5, 42, h=69)
    pdf.image(str(args.outputfolder / "images" / "uve.png"), 110, 42, h=69)
    pdf.image(str(args.outputfolder / "images" / "rlerawplot.png"),
              10.5,
              119,
              h=77.5)
    pdf.image(str(args.outputfolder / "images" / "rlenormplot.png"),
              10.5,
              199,
              h=77.5)
    pdf.output(str(args.outputfolder / "reports" / "norm_report.pdf"), 'F')

    if args.showbrowserqc:
        os.system(str(args.outputfolder / "reports" / "QC_inspection.pdf"))


def flagqc(args):
    infolanes = pd.read_csv(args.outputfolder / "info" / "rawinfolanes.csv",
        index_col='ID')

    flagged = set([])
    ### TODO ASK XABI ABOUT THIS LINE
    # (args.outputfolder / "reports" / "QCflags.txt").unlink(missing_ok=False)

    qc_flags_report = args.outputfolder / "reports" / "QCflags.txt"

    with open(str(qc_flags_report), "a") as f:
        args.current_state = '--> Starting QC flagging:'
        logging.info(args.current_state)

        if (
            all(_ >= args.minfov for _ in infolanes['FOV value']) and
            all(args.maxbd >= float(_) >= args.minbd
                for _ in infolanes.loc[:, 'Binding Density']) and
            all(_ is False for _ in infolanes.loc[:, 'limit of detection']) and
            all(i >= j for i, j in zip(infolanes.loc[:, '0,5fm'],
                                       infolanes.loc[:, 'Background'])) and
            all(args.maxscalingfactor > float(_) > args.minscalingfactor
                for _ in infolanes.loc[:, 'scaling factor']) and
            all(_ <= args.pbelowbackground
                for _ in infolanes.loc[:, 'Genes below backg %'])):
            info = 'Appropiate QC values for all samples. \n'
            args.current_state = info
            logging.info(args.current_state)
            f.writelines(info + '\n')
        else:
            args.current_state =\
                'Inappropiate QC values in some samples, revise QC report'
            logging.warning(args.current_state)
            for i in infolanes.index:
                thisFOV = infolanes.at[i, 'FOV value']
                thisBD = infolanes.at[i, 'Binding Density']
                thisLOD = infolanes.at[i, 'limit of detection']
                thisBG = infolanes.at[i, 'Background']
                this05 = infolanes.at[i, '0,5fm']
                thisSF = infolanes.at[i, 'scaling factor']
                thisgbb = infolanes.at[i, 'Genes below backg %']
                if thisFOV < args.minfov:
                    fovinfo = f"Low FOV value ({thisFOV}) in {i}. " +\
                        "Sample is flagged/discarded. \n"
                    logging.warning(fovinfo)
                    f.writelines(fovinfo)
                    flagged.add(i)
                if thisBD > args.maxbd or thisBD < args.minbd:
                    bdinfo = "Wrong binding density value " +\
                        f"({thisBD}) in {i}. " +\
                        "Sample is flagged/discarded."
                    logging.warning(bdinfo)
                    f.writelines(bdinfo)
                    flagged.add(i)
                if thisLOD:
                    lodinfo = "Wrong limit of detection value " +\
                        f"({thisLOD}) in {i}. " +\
                        "Sample is flagged/discarded."
                    logging.warning(lodinfo)
                    f.writelines(lodinfo)
                    flagged.add(i)
                if thisBG > this05:
                    bginfo = f"Wrong 0.5fm value ({thisBG}) in {i}. " +\
                        "Sample is flagged/discarded."
                    logging.warning(bginfo)
                    f.writelines(bginfo)
                    flagged.add(i)
                if thisSF > 3 or thisSF < 0.3:
                    sfinfo = "Wrong scaling factor value " +\
                        f"({thisSF}) in {i}" +\
                        "Sample is flagged/discarded."
                    logging.warning(sfinfo)
                    f.writelines(sfinfo)
                    flagged.add(i)
                if thisgbb > args.pbelowbackground:
                    gbbinfo = "Wrong genes below background value " +\
                        f"({thisgbb}) in {i}. " +\
                        "Sample is flagged/discarded."
                    logging.warning(gbbinfo)
                    f.writelines(gbbinfo)
                    flagged.add(i)

    flaggeddf = pd.DataFrame(flagged, columns=['Flagged_samples'])
    flaggeddf.to_csv(str(args.outputfolder / 'otherfiles' / 'flagged.csv'), index=False)

    if len(flagged) >= 3:
        args.nbadlanes = f"{len(flagged)} badlanes detected, " +\
            "check {qc_flags_report}."
    elif 0 < len(flagged) < 3:
        args.nbadlanes = f"{len(flagged)} badlanes detected: {flagged}"
    elif len(flagged) == 0:
        args.nbadlanes = 'No bad lanes detected from QC'

    args.badlanes = flagged

    return flagged


def removelanes(autoremove, args):
    infolanes = pd.read_csv(
        str(args.outputfolder / 'info' / 'rawinfolanes.csv'),
        index_col='ID')
    dfgenes = pd.read_csv(
        str(args.outputfolder / 'otherfiles' / 'dfgenes.csv'),
        index_col=0).T
    manualremove = args.remove

    if autoremove is None:
        autoremove = set(manualremove)
    if args.manual_remove == True:
        if type(manualremove) == str:
            manualremove = set(manualremove.split())
        print(f"Se retiran manualmente las muestras: {manualremove}.")
        logging.info(f"Se retiran manualmente las muestras: {manualremove}.")

    if set(dfgenes.index.intersection(autoremove)) == autoremove:
        dfgenes11 = dfgenes.drop(list(autoremove), axis=0)
        infolanes = infolanes.drop(autoremove)
    else:
        dfgenes11 = dfgenes
        if args.laneremover or args.manual_remove:
            args.current_state = 'Lanes set to remove not present in analysis.'
            print('Error: ' + args.current_state)
            logging.error(args.current_state)

    dfgenes11 = dfgenes11.T

    pathinfolanes = args.outputfolder / 'info' / 'infolanes.csv'
    infolanes.to_csv(pathinfolanes, index=True)
    exportdfgenes(dfgenes11, args)

    return dfgenes11, infolanes


def exportfilrawcounts(rawfcounts, args):
    pathfraw = args.outputfolder / 'otherfiles' / 'rawfcounts.csv'
    rawfcounts.to_csv(pathfraw, index=True)


def rescalingfactor23(args):
    """Recalculate scaling factor.

    Scaling factor needs to be recalculated after removing samples excluded
    by QC inspection or by applying background correction
    """
    infolanes = pd.read_csv(
        args.outputfolder / 'info' / 'infolanes.csv',
        index_col=0)

    if args.tecnormeth in ["posgeomean", "regression"]:
        use = 'posGEOMEAN'
        if args.tecnormeth == 'regression':
            negs = pd.read_csv(
                str(args.outputfolder / 'otherfiles' / 'dfnegcount.csv'), index_col=0)
            negs = negs.drop('maxoutlier', axis=1)
            corrected_negs = regretnegs(negs, args)
            backgr_regr = []
            for i in corrected_negs.index:
                backgr_regr.append(
                    (np.mean(corrected_negs.loc[i])) +
                    2 * (np.std(corrected_negs.loc[i]))
                )

            infolanes['backgr_regr'] = backgr_regr

    elif args.tecnormeth == 'Sum':
        use = 'Sum'
    elif args.tecnormeth == 'Median':
        use = 'Median'
    ref = np.mean(infolanes[use])

    scalingf2 = []

    for i in infolanes[use]:
        scaling = ref / i
        scalingf2.append(scaling)

    infolanes['scaling factor2'] = scalingf2

    return infolanes


def reinfolanes(args):
    '''
    Whether background correction (transform low counts) or technorm happens first, background or scaling factor needs to be recalculated after.
    For background correction, new background should be calculated from dfgenes given by technorm
    For technorm, new scaling factor needs to be calculated from dfgenes given by transformlowcounts
    '''

    dfgenes = pd.read_csv(args.outputfolder / 'otherfiles' / 'dfgenes.csv', index_col=0)
    infolanes = pd.read_csv(args.outputfolder / 'info' / 'infolanes.csv', index_col=0)
    dfneg = pd.read_csv(
        args.outputfolder / 'otherfiles' / 'dfnegcount.csv',
        index_col=0)
    dfpos = pd.read_csv(
        args.outputfolder / 'otherfiles' / 'posnegcounts.csv',
        index_col=0)
    dfneg = dfneg.drop(['maxoutlier'], axis=1)

    negnames = list(dfneg.columns)
    dfneg = dfgenes.loc[negnames]

    for i in infolanes.index:
        infolanes.at[i, 'Background'] = dfneg[i].mean()+(2*(dfneg[i].std()))
        infolanes.at[i, 'Background2'] = dfneg[i].max()
        infolanes.at[i, 'Background3'] = dfneg[i].mean()

    posnames = list(dfpos[dfpos['CodeClass'] == 'Positive'].index)
    dfpos = dfgenes.loc[posnames]

    for i in infolanes.index:
        infolanes.at[i, 'posGEOMEAN'] = gmean(dfpos[i])
        infolanes.at[i, 'Sum'] = dfpos[i].sum()
        infolanes.at[i, 'Median'] = np.median(dfpos[i])
        infolanes.at[i, 'meanexpr'] = gmean(dfgenes[i])

    infolanes = rescalingfactor23(args)

    pathoutinfolanes(infolanes, args)

def regretnegs(negs, args):
    corrected_negs = pd.DataFrame()
    posneg = pd.read_csv(
        args.outputfolder / 'otherfiles' / 'posnegcounts.csv',
        index_col=0)
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

    infolanes = pd.read_csv(args.outputfolder / 'info' / 'infolanes.csv')
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
            corrected = correct(
                thisseq_abc=thiseq_abc,
                meaneq_abc=meaneq_abc,
                counts=k)

            thissample.append(corrected)

        corrected_negs[str(i)] = thissample

    return corrected_negs.T



    # infolanes = pd.read_csv(str(args.outputfolder) + '/info/infolanes.csv', index_col=0)
    #
    # dfgenes = pd.read_csv(str(args.outputfolder) + '/otherfiles/rawfcounts.csv', index_col='Name')

def findaltnegatives(args):
    """Find alt negatives.

    To find low and stably expressed genes through the endogenous,
      to use as negative controls.
    In case native neg controls are not robust generates new infolanes with
    background alt in it.
    """

    infolanes = pd.read_csv(
        args.outputfolder / 'info' / 'infolanes.csv',
        index_col=0)
    dfgenes = pd.read_csv(
        args.outputfolder / 'otherfiles' / 'rawfcounts.csv',
        index_col=0).T
    genmean = dfgenes.mean()
    meangenmean = np.mean(genmean)
    genmean = genmean/meangenmean

    genmean.sort_values(inplace=True)

    genstd = dfgenes.std()
    genstd = genstd / meangenmean
    genstd = genstd * 2

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

    infolanes = pd.read_csv(args.outputfolder / 'info' / 'infolanes.csv', index_col=0)
    for i in infolanes.index:
        j = float(infolanes.at[i, 'scaling factor2'])
        dfgenes[i] = dfgenes[i] * j

    return dfgenes


def regresion(dfgenes, args):

    posneg = pd.read_csv(
        args.outputfolder / 'otherfiles' / 'posnegcounts.csv',
        index_col=0)
    for i in posneg.index:
        if 'NEG' in i:
            posneg = posneg.drop(i, axis=0)

    posneg = posneg.drop(['CodeClass'], axis=1)

    posneg['mean'] = posneg.mean(axis=1)

    xmean = sorted(posneg['mean'])
    ymean = [31, 125, 500, 2000, 8000, 32000]

    meaneq_abc = np.polynomial.Polynomial.fit(xmean, ymean, 3)

    infolanes = pd.read_csv(args.outputfolder / 'info' / 'infolanes.csv', index_col=0)
    for i in infolanes.index:
        thisx = sorted(posneg[i])
        thiseq_abc = np.polynomial.Polynomial.fit(thisx, ymean, 3)

        def correct(thisseq_abc, meaneq_abc, counts):
            this = thisseq_abc(counts)
            solutions = (meaneq_abc - this).roots()
            corrected = abs(min(solutions, key=lambda x: abs(x - counts)))
            return corrected

        for k,n in zip(dfgenes[i],dfgenes.index):
            corrected = correct(thisseq_abc=thiseq_abc, meaneq_abc=meaneq_abc, counts=k)
            dfgenes.at[n,i] = corrected

    return dfgenes


def transformlowcounts(dfgenes, args):

    infolanes = pd.read_csv(args.outputfolder / 'info' / 'infolanes.csv', index_col=0)
    varback = args.lowcounts
    varbg = args.background

    if args.background == 'Manual':
        mvarbg = args.manualbackground

        if varback == 'skip':
            pass

        elif varback == 'asim':
            dfgenes[dfgenes <= mvarbg] = mvarbg
        elif varback == 'subtract':
            dfgenes = dfgenes.sub(mvarbg)
            dfgenes[dfgenes <= 0] = 0.01


    else:
        if varback == 'skip':
            pass
        elif varback == 'asim':
            for i in infolanes.index:
                estebg = infolanes.at[i, varbg]
                dfgenes[i][dfgenes[i] <= estebg] = estebg
        elif varback == 'subtract':
            for i in infolanes.index:
                estebg = infolanes.at[i, varbg]
                dfgenes[i] = dfgenes[i].sub(estebg)
            dfgenes[dfgenes <= 0] = 0.01

    exportdfgenes(dfgenes, args)
    return dfgenes


def exporttnormgenes(normgenes, args):
    pathnormgenes = args.outputfolder / 'otherfiles' / 'tnormcounts.csv'
    normgenes.to_csv(pathnormgenes, index=True)


def getallhkes(args):
    dfhkegenes = pd.read_csv(
        args.outputfolder / 'otherfiles' / 'dfhkecount.csv',
        index_col=0)

    allhkes = dfhkegenes.T

    return allhkes


def filter50chkes(allhkes, args):
    """Filter housekeeping genes with less than 50 counts."""
    infolanes = pd.read_csv(args.outputfolder / 'info' / 'infolanes.csv')
    selhkes = pd.DataFrame()

    for i in infolanes['ID']:
        selhkes = allhkes.loc[allhkes[i] >= args.mincounthkes, :]

    selhkes = selhkes.round(decimals=3)
    # selhkes = selhkes.drop(['CodeClass', 'Accession'], axis=1)
    return selhkes

def findrefend(args, selhkes):
    """Find endogenous that can be used as reference genes."""
    dfgenes = pd.read_csv(args.outputfolder / 'otherfiles' / 'tnormcounts.csv', index_col=0)

    dfnamescodeclass = pd.read_csv(args.outputfolder / 'otherfiles' / 'dfgenes_raw.csv', usecols=['Name', 'CodeClass'])
    endogenousnames = dfnamescodeclass['Name'][dfnamescodeclass['CodeClass'].str.contains('Endogenous')]
    endogenousnames = endogenousnames.to_list()
    norm2end2 = dfgenes.loc[endogenousnames]
    if args.refendgenes == 'endhkes':
        endge = FindERG(norm2end2)
        # n best endogenous to include as reference genes
        bestend = list(endge[0:args.numend])
        logging.info(f"Most promising endogenous genes: {bestend}")
        print('Most promising endogenous genes: ', bestend)

        refgenesnames = selhkes.index.tolist() + bestend
        refgenes = dfgenes.loc[refgenesnames]
        refgenes = filter50chkes(refgenes, args)
    else:
        refgenesnames = selhkes.index.tolist()
        refgenes = dfgenes.loc[selhkes.index.tolist()]

    return refgenes


def pathoutrefgenes(refgenes, args):
    pathrefgenesview = args.outputfolder / 'otherfiles' / 'refgenesview.csv'
    refgenes.to_csv(pathrefgenesview, header=True, index=True)
    pathrefgenes = args.outputfolder / 'otherfiles' / 'refgenes.csv'
    refgenes.to_csv(pathrefgenes, header=True, index=True)


def getgroups(args):
    refgenes = pd.read_csv(
        args.outputfolder / 'otherfiles' / 'refgenes.csv',
        index_col=0).T
    flagged = pd.read_csv(args.outputfolder / 'otherfiles' / 'flagged.csv')
    flagged = set(flagged['Flagged_samples'])
    dfgroups = pd.read_csv(args.groupsfile, header=0, index_col=0)
    if args.laneremover:
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
    """Kruskal wallis calculation.
    Takes dfa-like dataframes, groups with samples at y and ref genes at x.
    A list of dataframes, a dataframe for each condition
    """
    lk = {}
    for i in args[0].columns: #for each candidate reference gene, column of the first dataframe
        la = []
        for j in args: #for each dataframe
            b = j[i] #we take the column referred to the gene of interest
            la.append(b)
        try:
            krus = stats.kruskal(*la)
        except Exception:
            pass
        lk[i] = krus

    lk = pd.DataFrame.from_dict(lk)
    lk = lk.rename(index={0: 'Result', 1: 'pvalue'})

    return lk


def calwilco(dfa, dfb):
    """Calculate Wilcoxon for every pair of groups."""
    lw = {}
    for i in dfa: #for each dataframe of the pair
        a = dfa[i]
        b = dfb[i]
        k = stats.ranksums(a, b)
        lw[i] = k
    lw = pd.DataFrame.from_dict(lw)
    lw = lw.rename(index={0: 'Result', 1: 'pvalue'})
    return lw


def calwilcopairs(*ddfc):  # perform wilcoxon calculations for each par of conditions
    lenargs = np.arange(0,len(ddfc))

    lw = {}
    for i in lenargs:
        for j in lenargs:
            if i < j:
                w = calwilco(ddfc[i][1], ddfc[j][1])
                pair = f"wilcox: {ddfc[i][0]} / {ddfc[j][0]}"
                lw[pair] = w
    return lw


def flagkrus(reskrus):
    flaggedgenes = []
    for i in reskrus:
        if reskrus.loc['pvalue', i] < 0.05:
            flaggedgenes.append(i)
    return flaggedgenes

def filterkruskal(flaggedgenes, args):
    refgenes = pd.read_csv(
        args.outputfolder / 'otherfiles' / 'refgenes.csv',
        index_col=0).T
    if args.filtergroupvariation == 'filterkrus':
        if (len(refgenes.columns) - len(flaggedgenes)) <= 2:
            args.current_state = \
                "Too many genes to be removed from kruskal filtering, " +\
                "consider using another refgenes or change settings to " +\
                "'flagkrus'."
            logging.warning(args.current_state)
        else:
            if len(flaggedgenes)>=1:
                refgenes = refgenes.drop(columns=flaggedgenes)
            else:
                refgenes = refgenes
    elif args.filtergroupvariation == 'flagkrus':
        args.current_state = \
            f"Genes not recommended as refgenes by kruskal: {flaggedgenes}."
        logging.warning(args.current_state)
    pathrefgenes = args.outputfolder / 'otherfiles' / 'refgenes.csv'
    refgenes.T.to_csv(pathrefgenes, header=True, index=True)

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
    refgenes = pd.read_csv(
        args.outputfolder / 'otherfiles' / 'refgenes.csv',
        index_col=0).T
    if args.filtergroupvariation == 'filterwilcox':
        if (len(flaggedwilcox) < len(refgenes.columns) and
            (len(refgenes.columns) - len(flaggedwilcox)) > 2):
            refgenes = refgenes.drop(columns=flaggedwilcox)
        else:
            args.current_state = \
                "Too many ref genes filtered by wilcox, skipping filtering." +\
                " Consider re-designing reference/housekeeping genes."
            logging.warning(args.current_state)
    elif args.filtergroupvariation == 'flagwilcox':
        args.current_state = \
            f"Genes not recommended as refgenes by wilcoxon: {flaggedwilcox}."
        logging.warning(args.current_state)
    pathrefgenes = args.outputfolder / 'otherfiles' / 'refgenes.csv'
    refgenes.T.to_csv(pathrefgenes, header=True, index=True)

    return refgenes


def measureM(df, ctVal=False):
    if ctVal:
        for column in df:
            minimo = min(df[column])
            df[column] = df[column] - minimo
        df = 2**-df
    else:
        for column in df:
            maximo = max(df[column])
            df[column] = df[column] / maximo
    # XXX The naming here is insane
    m = list(df.index)

    n = list(df.columns)

    M_a = pd.DataFrame()
    for j in n:
        Vjk = []
        for k in n:
            Ajk = []
            for i in m:
                if df.at[i,k] == 0:
                    df.at[i,k] = 0.000001
                Ajk.append(np.log2(df.at[i, j] / df.at[i, k]))
            Vjk.append(np.std(Ajk, ddof=1))


        M_a.at[j, 'Genes'] = j

        M_a.at[j, 'M'] = np.sum(Vjk) / (len(n) - 1)
    M_a = M_a.sort_values('M', ascending=False)
    return M_a


def geNorm(df, avgm=pd.DataFrame()):
    result = measureM(df)
    n = len(df.columns)
    if n <= 2:
        bestgen = result.iat[0, 0]
        newrow = pd.DataFrame([[bestgen, result['M'].mean()]])
        avgm = pd.concat([avgm, newrow])
        lastrow = pd.DataFrame([[result.iat[1, 0], result.iat[1, 1]]])
        avgm = pd.concat([avgm, lastrow])
        newindex2 = np.arange(start=1, stop=len(avgm[0]) + 1)
        avgm.index = newindex2
        return avgm
    else:
        bestgen = result.iat[0, 0]
        newrow = pd.DataFrame([[bestgen, result['M'].mean()]])
        avgm = pd.concat([avgm, newrow])
        newdf = df.drop(bestgen, axis=1)
        return geNorm(newdf, avgm)


def pairwiseV(datarefgenes):
    Vs = pd.DataFrame()

    buf = geNorm(datarefgenes)
    n = len(buf[0])
    m = np.arange(start=0, stop=n - 2)

    a = 2
    for i in m:
        mas = np.arange(start=1, stop=a + 1)
        genes2 = []
        genes3 = []
        for z in mas:
            gen = buf.iloc[n-z, 0]
            genes2.append(gen)
            genes3.append(gen)
        lmas = len(mas)+1
        genes3.append(buf.iloc[n-lmas, 0])

        df2 = datarefgenes.loc[:, genes2]
        df3 = datarefgenes.loc[:, genes3]

        medias2 = []
        medias3 = []
        for j in df2.index:
            media2 = gmean(df2.loc[j])
            medias2.append(media2)
            media3 = gmean(df3.loc[j])
            medias3.append(media3)

        logmedias = []
        for k, m in zip(medias2, medias3):
            n1 = k / m
            An_n1 = np.log2(n1)
            logmedias.append(An_n1)

        Vn_n1 = np.std(logmedias, ddof=1)

        newvn = pd.DataFrame([[f'V{a}/V{a+1}', Vn_n1]])
        Vs = pd.concat([Vs, newvn])

        newindex = np.arange(start=1, stop=len(Vs[0]) + 1)
        Vs.index = newindex

        a += 1

    return Vs


def ploteme(eme, args):
    plt.figure()
    plt.plot(eme['Genes'], eme['M'], 'o')
    plt.xticks(rotation=45)
    plt.xlabel('refgenes')
    plt.ylabel('measured M')
    plt.xticks(rotation=45)
    plt.title('measure M')
    plt.savefig(args.outputfolder / 'images' / 'eme.png')
    plt.close()


def plotavgm(genorm, args):
    plt.figure()
    plt.plot(genorm[0], genorm[1], 'o')
    plt.xticks(rotation=45)
    plt.xlabel('refgenes')
    plt.ylabel('Avg. M')
    plt.xticks(rotation=45)
    plt.title('Genorm result')
    plt.savefig(args.outputfolder / 'images' / 'avgm.png')
    plt.close()


def plotuve(uve, args):
    plt.figure()
    plt.bar(uve[0], uve[1])
    plt.xticks(rotation=45)
    plt.xlabel('gen pairs')
    plt.ylabel('pairwise variation')
    plt.xticks(rotation=45)
    plt.title('Pairwise variation')
    plt.savefig(args.outputfolder / 'images' / 'uve.png')
    plt.close()


def getnrefgenes(uvedf):
    minuve = uvedf[1].min()
    iminuve = uvedf[uvedf[1] == minuve].index
    iminuve = iminuve[0] + 2

    return iminuve


def getnamesrefgenes(uvedf, genorm, args):
    if args.nrefgenes is not None:
        n = args.nrefgenes
    elif args.nrefgenes is None:
        n = getnrefgenes(uvedf)

    names = []
    for i in np.arange(start=1, stop=n + 1):
        a = genorm.loc[i, 0]
        names.append(a)
    return names


def takerefgenes(names, args):
    datarefgenes = pd.read_csv(args.outputfolder / 'otherfiles' / 'refgenes.csv', index_col=0)
    bestrefgenes = datarefgenes.loc[names]
    bestrefgenes.to_csv(args.outputfolder / 'otherfiles' / 'bestrefgenes.csv', index=True)

    return bestrefgenes


def rankfeaturegenes(data, targets, args, verbose=0):
    """Rank feature genes.

    'data' can be refgenes (usual, fast exploration)
      or all genes (very large analysis, several hours)
      using all genes to further visualization
    'targets' must be single column sample-class association
    'num_neighbors' can simplify analysis, default 5, shouldnt be lower than 3
    """
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

    targets = targets.loc[data.index, :]
    data.index = targets.index

    X = data
    y = targets

    sfs1 = SFS(knn,
               k_features=1,
               forward=False,
               floating=False,
               verbose=0,
               scoring='balanced_accuracy',
               cv=2,
               n_jobs=-1)
    sfs1 = sfs1.fit(X, y.values.ravel())

    metrics = pd.DataFrame.from_dict(sfs1.get_metric_dict()).T

    return metrics


def rankstatsrefgenes(reskrus, reswilcopairs):
    """Create a summary dataframe of confidence.

    Takes info from kruskal and wilcoxon to create a summary dataframe
    of confidence for refgenes
    """
    ranking = pd.DataFrame()
    ranking['Genes'] = list(reskrus.columns)
    ranking.set_index('Genes', inplace=True)
    ranking['Kruskal p-value'] = reskrus.loc['pvalue', :]
    ranking.sort_values('Kruskal p-value', ascending=True, inplace=True)

    for i, j in reswilcopairs.items():
        ranking[i] = j.loc['pvalue', :]

    return ranking


def getallgenesdf(args):
    df = pd.read_csv(args.outputfolder / 'otherfiles' / 'tnormcounts.csv', index_col=0)
    return df


def gettopngenesdf(args):
    df = pd.read_csv(args.outputfolder / 'otherfiles' / 'tnormcounts.csv', index_col=0)
    df['mean'] = df.mean(axis=1)
    df.sort_values(by=['mean'], ascending=False, inplace=True)
    df = df[:args.topngenestocontnorm]
    df.drop('mean', inplace=True, axis=1)

    return df

def getnormfactor(refgenesdf, eme, args):
    geomeans1 = {}

    if args.contnorm == 'ponderaterefgenes':
        eme2 = eme

        eme2['M'] = eme2['M'] * len(eme['M']) / sum(eme['M'])
        
        for i in refgenesdf:
            
            a = (refgenesdf[i] * eme2['M']).sum() / eme2['M'].sum()
            geomeans1[i] = a

    else:
        for i in refgenesdf.columns:
            igeomean = gmean(refgenesdf[i])
            geomeans1[i] = igeomean
    
    prenormfactor = np.mean(list((geomeans1.values())))
    normfactor = {}
    for i, j in geomeans1.items():
        nfactor = prenormfactor / j
        normfactor[i] = nfactor

    return normfactor


def refnorm(normfactor, args):
    df = pd.read_csv(args.outputfolder / 'otherfiles' / 'tnormcounts.csv', index_col=0)

    for i in df.columns:
        nfactor = normfactor[i]
        df[i] = df[i]*nfactor

    return df


def pathoutrnormgenes(df, args):
    pathrnorm = args.outputfolder / 'results' / 'rnormcounts.csv'
    df.to_csv(pathrnorm)


def pathoutadnormgenes(df, args):
    pathadnorm = args.outputfolder / 'otherfiles' / 'adnormcounts.csv'
    df.to_csv(pathadnorm)


def adnormalization(df, args):
    scaler = MinMaxScaler()
    stdnormgenes = scaler.fit_transform(df)
    stdnormgenes = pd.DataFrame(stdnormgenes)
    stdnormgenes.index = df.index
    stdnormgenes.columns = df.columns
    pathoutadnormgenes(stdnormgenes, args)
    return stdnormgenes


def logarizeoutput(rnormgenes, args):
    if 'group' in rnormgenes.index:
        rnormgenes.drop('group', axis=0, inplace=True)
    with np.errstate(divide = 'ignore'):
        if args.logarizedoutput == '2':
            logarizedgenes = rnormgenes.applymap(lambda x: np.log2(x))
        if args.logarizedoutput == '10':
            logarizedgenes = rnormgenes.applymap(lambda x: np.log10(x))

    pathlogarized = args.outputfolder / 'otherfiles' / 'logarized_rnormcounts.csv'
    logarizedgenes.to_csv(pathlogarized)

    return logarizedgenes

def RLEcal(rnormgenes, args):

    rlegenes = rnormgenes

    rlegenes = rlegenes.T

    median = rlegenes.median()
    rlegenes = rlegenes - median
    rlegenes = rlegenes.T
    if 'group' in rnormgenes.index:
        rlegenes.loc['group'] = rnormgenes.loc['group']
    return rlegenes


def getmeaniqr(rlegenes):
    iqrlist = []
    if 'group' in rlegenes.index:
        rlegenesng = rlegenes.drop('group', axis=0)
    else:
        rlegenesng = rlegenes
    for i in rlegenesng:
        iqr = stats.iqr(rlegenesng[i])
        iqrlist.append(iqr)
    meaniqr = np.mean(iqrlist)
    
    return meaniqr


def plotevalnorm(matrix, what, meaniqr, args):
    matrix = matrix.apply(lambda x: x + 1)
    matrix = matrix.applymap(lambda x: np.log10(x))

    matrix = matrix.T
    for i in matrix.columns:
        median = np.median(matrix[i])
        matrix[i] = matrix[i] - median

    matrix = matrix.T

    estoo = matrix

    plt.figure(figsize=(30, 12))
    sns.boxplot(data=matrix, showfliers=False, showmeans=True)
    plt.title(what + '. IQR: ' + str(meaniqr), fontsize=40)
    plt.tick_params(
        axis='x',
        which='both',
        bottom=False,
        top=False,
        labelbottom=False)
    plt.ylim(-1, 1)
    plt.ylabel('RLE', fontsize=36)
    plt.xlabel('Samples', fontsize=40)
    sns.stripplot(data=matrix, size=2, palette='dark:black')
    plt.savefig(args.outputfolder / 'images' / 'rlenormplot.png')
    plt.savefig(args.outputfolder / 'images' / 'rlenormplot2.png', dpi=17)
    plt.close()


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
    plt.tick_params(
        axis='x',
        which='both',
        bottom=False,
        top=False,
        labelbottom=False)
    plt.title(f"{what}. IQR: {meaniqrraw}", fontsize=40)
    plt.ylabel('RLE', fontsize=36)
    plt.xlabel('Samples', fontsize=40)
    plt.ylim(-1, 1)
    sns.stripplot(data=matrix, size=2, palette='dark:black')

    plt.savefig(args.outputfolder / 'images' / 'rlerawplot.png')
    plt.savefig(args.outputfolder / 'images' / 'rlerawplot2.png', dpi=17)
    plt.close()

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
    parser.add_argument('-lc', '--lowcounts', type=str, default='skip', choices=['skip', 'asim', 'subtract'],  help='what to do with counts below background?')
    parser.add_argument('-mi', '--modeid', type=str, default='filename', choices=['sampleID','filename', 'id+filename'], help='choose sample identifier. sampleID: optimal if assigned in rccs. filenames: easier to be unique. id+filename: care with group assignment coherence')
    parser.add_argument('-mv', '--modeview', type=str, default='view', choices=['justrun', 'view'], help='choose if plot graphs or just run calculations')
    parser.add_argument('-tnm', '--tecnormeth', type=str, default='posgeomean', choices=['posgeomean','Sum', 'Median', 'regression'], help='choose method for technical normalization')
    parser.add_argument('-reg', '--refendgenes', type=str, default= 'endhkes', choices=['hkes', 'endhkes'], help='choose refgenes, housekeeping, or hkes and endogenous')
    parser.add_argument('-re', '--remove', type=str, nargs='+', default=None, help='lanes to be removed from the analysis')
    parser.add_argument('-bg', '--background', type=str, default= 'Background', choices=['Background', 'Background2', 'Background3', 'Backgroundalt', 'Manual'], help='choose background: b1=meancneg+(2*std), b2=maxcneg, b3=meancneg, balt=uses alternative subset of negative controls')
    parser.add_argument('-pbb', '--pbelowbackground', type=int, default=85, help='if more than %bb genes are below background, sample gets removed from analysis')
    parser.add_argument('-mbg', '--manualbackground', type=float, default=None, help='set manually background')
    parser.add_argument('-crg', '--chooserefgenes', type=list, nargs='+', default = None, help = 'list of strings like. choose manualy reference genes to use over decided-by-program ones')
    parser.add_argument('-fgv', '--filtergroupvariation', type=str, default='filterkrus', choices=['filterkrus', 'filterwilcox', 'flagkrus', 'flagwilcox', 'nofilter'], help='¿filter or flag preselected ref genes by significative group-driven differences? needs groups to be declared')
    parser.add_argument('-fsn', '--featureselectionneighbors', type=float, default=4, help='number of neighbors for feature selection analysis of refgenes. recommended 3-6')
    parser.add_argument('-g', '--groups', type=str, default='yes', choices=['yes','no'], help='defining groups for kruskal/wilcox/fs analysis?')
    parser.add_argument('-ne', '--numend', type=int, default=6, help='number of endogenous tofind by ERgene to include in analysis to check viability as refgenes')
    parser.add_argument('-ar', '--autorename', type=bool, default=False, help='turn on when sample IDs are not unique, be careful on sample identification detail')
    parser.add_argument('-cn', '--contnorm', type=str, default='refgenes', choices=['ponderaterefgenes', 'refgenes', 'all', 'topn'])
    parser.add_argument('-an', '--adnormalization', type=str, default='no', choices=['no', 'standarization'], help='perform additional normalization? standarization available')
    parser.add_argument('-tn', '--topngenestocontnorm', type=int, default=100, help='set n genes to compute for calculating norm factor from top n expressed endogenous genes')
    parser.add_argument('-mch', '--mincounthkes', type=int, default=80, help='set n min counts to filter hkes candidate as refgenes')
    parser.add_argument('-nrg', '--nrefgenes', type=int, default=None, help='set n refgenes to use, overwriting geNorm calculation')
    parser.add_argument('-lr', '--laneremover', type=bool, default=True, choices=[True, False], help='option to perform analysis with all lanes if set to no')
    parser.add_argument('-lo', '--logarizedoutput', type=str, default='no', choices=['2', '10', 'no'], help='want normed output to be logarized? in what logbase?')
    parser.add_argument('-le', '--logarizeforeval', type=str, default='10', choices=['2', '10', 'no'], help= 'logarithm base for RLE calculations')
    parser.add_argument('-gf', '--groupsfile', type=str, default='../examples/groups_d1_COV_GSE183071.csv', help='enter file name where groups are defined')
    parser.add_argument('-st', '--start_time', type=float, default = time.time())
    parser.add_argument('-cs', '--current_state', type=str, default='Ready')
    parser.add_argument('-ftl', '--tnormbeforebackgcorr', type=int, default=1, help='0= False, 1= True, 2= ruvg')
    parser.add_argument('-of', '--outputfolder', type=str, default=tempfile.gettempdir() + '/guanin_output')
    parser.add_argument('-sll', '--showlastlog', type=bool, default = False)
    parser.add_argument('-rgs', '--refgenessel', type=list, default = [])
    parser.add_argument('-im2', '--indexmethod2', type=int, default=0)
    parser.add_argument('-k', '--kvalue', type=int, default=3)
    parser.add_argument('-pip', '--pipeline', type=str, default='ruvgnorm', choices=['ruvgnorm', 'scalingfactors'])
    parser.add_argument('-wrg', '--whatrefgenes', type=list, default=[])
    parser.add_argument('-m', '--eme', type=object, default=None)
    parser.add_argument('-gpca', '--grouppca', type=str, default='GROUP')
    parser.add_argument('-dm', '--deseq2_mor', type=bool, default=True)
    parser.add_argument('-mr', '--manual_remove', type=bool, default=False)
    parser.add_argument('-nbl', '--nbadlanes', type=str, default='No badlanes detected')
    parser.add_argument('-bl', '--badlanes', type=set, default=set())
    parser.add_argument('-e', '--elapsed', type=float, default=0.0)
    return parser.parse_args()


#################BIG BLOCKS -- BUTTONS
def showinfolanes(args):
    """Load RCCS and show infolanes."""

    infolanes, dfgenes, dfnegcount, dfhkecount, dfposneg = loadrccs(args)

    createoutputfolder(args)
    exportrawcounts(dfgenes, args)

    dfgenes2 = dfgenes.drop(['CodeClass', 'Accession'], axis=1)
    dfgenes2 +=1

    exportdfgenes(dfgenes2, args)
    exportrawinfolanes(infolanes, dfnegcount, dfhkecount, dfposneg, args)
    summarizerawinfolanes(args)


    infolanes = infolanes.style.applymap(
        condformat,
        top=args.maxfov,
        bot=args.minfov,
        subset='FOV value')
    infolanes = infolanes.applymap(
        condformat,
        top=args.maxbd,
        bot=args.minbd,
        subset='Binding Density')
    infolanes = infolanes.applymap(
        condformat_LOD,
        subset='limit of detection')
    infolanes = infolanes.applymap(
        condformat,
        top=args.maxlin,
        bot=args.minlin,
        subset='R2')
    infolanes = infolanes.applymap(
        condformat,
        top=args.pbelowbackground,
        bot=0,
        subset='Genes below backg %')
    infolanes = infolanes.applymap(
        condformat,
        top=args.maxscalingfactor,
        bot=args.minscalingfactor,
        subset='scaling factor')

    infolanes.to_html(args.outputfolder / 'info' / 'rawinfolanes.html')

    if args.showbrowserrawqc:
        ##DEBUG
        #FOR SOME REASON IT POPS 2 TABS
        webbrowser.open(str(args.outputfolder / 'info' / 'rawinfolanes.html'))

def plotandreport(args, whatinfolanes="rawinfolanes", rawreport=True):
    """Plot QC from infolanes or raw info lanes (file).

    :param whatinfolanes: chooses file rawinfolanes.csv or infolanes.csv
    :return:
    """
    if whatinfolanes == 'rawinfolanes':
        infolanes = pd.read_csv(
            args.outputfolder / 'info' / 'rawinfolanes.csv',
            index_col='ID')
    elif whatinfolanes == 'infolanes':
        infolanes = pd.read_csv(
            args.outputfolder / 'info' / 'infolanes.csv',
            index_col=0)

    dfnegcount = pd.read_csv(
        args.outputfolder / 'otherfiles' / 'dfnegcount.csv',
        index_col=0)
    dfhkecount = pd.read_csv(
        args.outputfolder / 'otherfiles'  / 'dfhkecount.csv',
        index_col=0)

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
    logging.info(args.current_state)
    pdfreport(args, rawreport=rawreport)


def runQCview(args):
    args.start_time = time.time()
    try:
        showinfolanes(args)
        state = 'All RCCs loaded succesfully'
        args.current_state = state
        logging.info(state)
    except Exception as e:
        state = "Something went wrong loading files, check input folder. " +\
            f"Error: {e}"
        traceback.print_exception()
        args.current_state = state
        logging.error(args.current_state)
        return
    try:
        plotandreport(args)
        state = "Data loaded succesfuly, preliminary analysis and plots " +\
            "ready to inspect in reports output folder"
        args.current_state = state
        logging.info(state)
    except Exception as e:
        state = "Something went wrong with preliminary analysis and/or " +\
            f"plotting. Error: {e}"
        args.current_state = state
        logging.error(args.current_state)
    args.elapsed = time.time() - args.start_time

    args.current_state = f"Elapsed loading RCCs {args.elapsed} seconds"


def runQCfilterpre(args):
    """Filter the result of runQC.
    """
    args.start_time = time.time()
    dfgenes = pd.read_csv(args.outputfolder / 'otherfiles' / 'dfgenes_raw.csv', index_col=0)
    flagged = flagqc(args)

    if args.laneremover:
        dfgenes, infolanes = removelanes(flagged, args)
        exportfilrawcounts(dfgenes, args)

    elif not args.laneremover:
        dfgenes, infolanes = removelanes([], args)
        exportfilrawcounts(dfgenes, args)

    exportdfgenes(dfgenes, args)
    exportdfgenes_qc(dfgenes, args)

    infolanes = rescalingfactor23(args)
    pathoutinfolanes(infolanes, args)

    infolanes = findaltnegatives(args)
    pathoutinfolanes(infolanes, args)

    html_infolanes(args, infolanes)

    summarizeinfolanes(args)
    plotandreport(args, whatinfolanes='infolanes', rawreport=False)

    args.elapsed += time.time() - args.start_time
    args.current_state = f"Elapsed performign QC analysis {args.elapsed} seconds"


def html_infolanes(args, infolanes):
    infolanes = infolanes.style.applymap(
        condformat,
        top=args.maxfov,
        bot=args.minfov,
        subset='FOV value')
    infolanes = infolanes.applymap(
        condformat,
        top = args.maxbd,
        bot=args.minbd,
        subset='Binding Density')
    infolanes = infolanes.applymap(
        condformat_LOD,
        subset='limit of detection')
    infolanes = infolanes.applymap(
        condformat,
        top=args.maxlin,
        bot=args.minlin,
        subset='R2')
    infolanes = infolanes.applymap(
        condformat,
        top=args.pbelowbackground,
        bot=0,
        subset='Genes below backg %')
    infolanes = infolanes.applymap(
        condformat,
        top=args.maxscalingfactor,
        bot=args.minscalingfactor,
        subset='scaling factor')

    infolanes.to_html(args.outputfolder / 'info' / 'infolanes.html')

    if args.showbrowserqc:
        webbrowser.open(args.outputfolder / 'info' / 'infolanes.html')

    return infolanes

def runQCfilter(args):
    try:
        runQCfilterpre(args)
        args.current_state = \
            'QC filter applied, ready to perform technical normalization'
        logging.info(args.current_state)
    except Exception as e:
        args.current_state =\
            "Unknown error while QC filtering, check input data and " +\
            f"parameters. Error: {e}"
        print('Error on line {}'.format(sys.exc_info()[-1].tb_lineno), type(e).__name__, e)
        logging.info(args.current_state)

def pipeline1(args):
    args.start_time = time.time()
    if args.pipeline == 'ruvgnorm' and args.deseq2_mor:
        apply_deseq2_mor(args)
    elif args.pipeline == 'scalingfactors':
        technorm(args)
    args.elapsed += time.time() - args.start_time
    args.current_state = f"Pre-normalization succesful. Elapsed: {args.elapsed} seconds"
    logging.info(args.current_state)

def pipeline2(args):
    args.start_time = time.time()
    if args.pipeline == 'ruvgnorm':
        selecting_refgenes(args)
        RUVgnorm2(args)
    elif args.pipeline == 'scalingfactors':
        selecting_refgenes(args)
        contnorm(args)
    args.elapsed += time.time() - args.start_time
    args.current_state = f"Elapsed performing content normalization: {args.elapsed} seconds"
    logging.info(args.current_state)

def RUVgnorm2(args, center=True, round=False, epsilon=1, tolerance=1e-8, isLog=False):
    warnings.filterwarnings("ignore", category=FutureWarning)
    gene_matrix = pd.read_csv(args.outputfolder / 'otherfiles' / 'tnormcounts.csv', index_col=0)

    refgenes = args.refgenessel

    k = args.kvalue
    drop = 0

    if isLog:
        Y = gene_matrix.T
    else:
        Y = gene_matrix.applymap(lambda x: math.log(x + epsilon)).T

    if center:
        Ycenter = Y.apply(lambda col: scale(col, with_mean=True, with_std=False), axis=0)
    else:
        Ycenter = Y

    if drop >= k:
        raise ValueError("'drop' must be less than 'k'.")

    svdWau, svdWad, svdWauv = svd(Ycenter.loc[:, gene_matrix.index.isin(refgenes)], full_matrices=False)
    first = 1 + drop
    k = min(k, max([i for i, val in enumerate(svdWad) if val > tolerance], default=0))
    W = svdWau[:, first-1:k]

    alpha = (W.T @ W) @ W.T @ Y
    correctedY = np.subtract(Y, W @ alpha)


    if not isLog:
        if round:
            correctedY = np.exp(correctedY).subtract(epsilon).round()
            correctedY[correctedY < 0] = 0
        else:
            correctedY = np.exp(correctedY).subtract(epsilon)

    rnormgenes = correctedY.T
    rngg = logarizeoutput(rnormgenes, args)
    rngg.to_csv(args.outputfolder / 'otherfiles' / 'rngg.csv', index=True)
    exportdfgenes(rnormgenes, args)
    pathoutrnormgenes(rnormgenes, args)

def technorm(args):
    dfgenes = pd.read_csv(args.outputfolder / 'otherfiles' / 'dfgenes_qc.csv', index_col=0)
    try:
        if not args.tnormbeforebackgcorr:
            print('Applying background correction...')
            dfgenes = transformlowcounts(dfgenes, args)
            reinfolanes(args)
            print('Performing technical normalization...')
            if args.tecnormeth != 'regression':
                normgenes = normtecnica(dfgenes, args)
            elif args.tecnormeth == 'regression':
                normgenes = regresion(dfgenes, args)
            exportdfgenes(normgenes, args)
        elif args.tnormbeforebackgcorr:
            print('Performing technical normalization...')
            if args.tecnormeth != 'regression':
                normgenes = normtecnica(dfgenes, args)
            elif args.tecnormeth == 'regression':
                normgenes = regresion(dfgenes, args)
            exportdfgenes(normgenes, args)
            print('Applying background correction...')

            reinfolanes(args)
            transformlowcounts(normgenes, args)

        exporttnormgenes(normgenes, args)


    except Exception as e:
        args.current_state = 'Failed technical normalization'
        logging.error(args.current_state)

def selecting_refgenes(args):

    allhkes = getallhkes(args)
    print('Housekeeping genes present in analysis: ', list(allhkes.index))

    selhkes = filter50chkes(allhkes, args)
    if len(selhkes.index) <= 2:
        selhkes = allhkes

        args.current_state = 'All or almost all housekeeping genes are low expressed. Consider re-design experiment. Proceeding with all hkes'
        logging.warning(args.current_state)
    else:
        args.current_state =\
            "Housekeeping genes with more than 50 counts for all lanes: " +\
            f"{list(selhkes.index)}"
        logging.info(args.current_state)

    try:
        refgenes = findrefend(args, selhkes)

        args.current_state =\
            "Refgenes in analysis including housekeepings + best " +\
            f"endogenous selected: {list(refgenes.index)}"
        logging.info(args.current_state)
    except Exception as e:
        logging.warning(
            f"Unable to retrieve candidate ref genes from endogenous, ERROR: {e}")
        refgenes = selhkes

    pathoutrefgenes(refgenes, args)

    if os.path.isfile(args.groupsfile):
        targets = pd.read_csv(args.groupsfile)
        groups = set(targets['GROUP'])
        if len(groups) > 1:
            args.groups = 'yes'

    if args.groups == 'yes':

        args.current_state = '--> Performing kruskal-wallis analysis'
        logging.info(args.current_state)
        ddf = getgroups(args)

        ddfc = list(ddf.items())

        ddfb = list(ddf.values())

        reskrus = calkruskal(*ddfb)


        args.current_state = '--> Performing wilcoxon analysis'
        logging.info(args.current_state)
        reswilcopairs = calwilcopairs(*ddfc)

        flaggedgenes = flagkrus(reskrus)
        flaggedwilcox = flagwilcox(reswilcopairs)
        flaggedwilcox = flagwilcox(reswilcopairs)

        flaggedboth = set(flaggedgenes).intersection(set(flaggedwilcox))

        logging.info(
            f"Flagged genes by kruskal-wallis and/or wilcoxon: {flaggedboth}")

        if args.filtergroupvariation == 'filterkrus':
            refgenes = filterkruskal(flaggedgenes, args)
        elif args.filtergroupvariation == 'filterwilcox':
            refgenes = filterwilcox(flaggedwilcox, args)
        else:
            refgenes = filterkruskal([], args)
        print("Ref genes present in analysis after applying kruskal-wallis " +\
              f"or wilcoxon filtering: {list(refgenes.columns)}")
    elif args.groups == 'no':
        pass

    datarefgenes = refgenes

    eme = measureM(datarefgenes)

    genorm = geNorm(datarefgenes)

    uve = pairwiseV(datarefgenes)

    print('--> Performing geNorm calculations')

    ploteme(eme, args)

    plotavgm(genorm, args)

    plotuve(uve, args)

    args.current_state = '--> Calculating optimal number of reference genes'
    print('--> Calculating optimal number of reference genes')

    if args.chooserefgenes is None:

        names = getnamesrefgenes(uve, genorm, args)

        if args.nrefgenes == None:
            args.current_state = 'Ref. genes selected (auto): ' + str(names)
            args.refgenessel = names
        elif args.nrefgenes != None:
            args.current_state = 'Ref. genes selected (n_manual): ' + str(names)
            args.refgenessel = names

    else:
        names = args.chooserefgenes


        args.current_state = f"Ref. genes selected (manual): {names}"
    args.refgenessel = names


    print('--> Performing feature selection for refgenes evaluation and control.')
    if args.groups == 'yes' or os.path.exists(args.groupsfile):
        if len(targets.columns) > 2:
            targets = targets[['SAMPLE', 'GROUP']].copy()

        metrics = rankfeaturegenes(datarefgenes, targets, args)
        ranking = rankstatsrefgenes(reskrus, reswilcopairs)

        ranking2 = ranking.style.applymap(condformat, top=1, bot=0.05, subset='Kruskal p-value')
        for i in ranking2.columns:
            ranking2 = ranking2.applymap(condformat, top=1, bot=0.05, subset=i)

        ranking2.to_html(args.outputfolder / 'info' / 'ranking_kruskal_wilcox.html')
        ranking.to_csv(args.outputfolder / 'info' / 'ranking_kruskal_wilcox.csv')

    if args.showbrowsercnorm == True:
        webbrowser.open(str(args.outputfolder / 'info' / 'ranking_kruskal_wilcox.html'))

    if args.groups == 'yes':


        metrics2 = metrics.style.applymap(
            condformat,
            top=1.5 / len(groups),
            bot=0.5 / len(groups),
            subset='avg_score')


        metrics2.to_html(
            args.outputfolder / 'reports' / 'metrics_reverse_feature_selection.html')
        metrics.to_csv(
            args.outputfolder / "reports" / "metrics_reverse_feature_selection.csv")

    if args.showbrowsercnorm and args.groups:
        webbrowser.open(
            str(args.outputfolder / 'reports' / "metrics_reverse_feature_selection.html"))

    args.eme = eme

def contnorm(args):

    logging.info('Starting content normalization')

    allgenes = getallgenesdf(args)

    names = args.refgenessel
    bestrefgenes = takerefgenes(names, args)
    eme = args.eme
    eme.set_index('Genes', drop=True, inplace=True)

    if args.contnorm == 'refgenes' or args.contnorm == 'ponderaterefgenes':
        normfactor = getnormfactor(bestrefgenes, eme, args)
    elif args.contnorm == 'all':
        normfactor = getnormfactor(allgenes, eme, args)
    elif args.contnorm == 'topn':
        topngenes = gettopngenesdf(args)
        normfactor = getnormfactor(topngenes, eme, args)

    rnormgenes = refnorm(normfactor, args)
    pathoutrnormgenes(rnormgenes, args)

    if args.adnormalization == 'standarization':
        adnormalization(rnormgenes, args)

    if args.logarizedoutput != 'no':
        rngg = logarizeoutput(rnormgenes, args)
    else:
        rngg = rnormgenes
    rngg.to_csv(str(args.outputfolder / 'otherfiles' / "rngg.csv"), index=True)

def apply_deseq2_mor(args):
    dfgenes = pd.read_csv(
        args.outputfolder / 'otherfiles' / 'dfgenes_qc.csv',
        index_col=0)
    counts, sizefactors = deseq2_norm(dfgenes)
    exporttnormgenes(counts, args)
    exportdfgenes(counts, args)

def plotpcaraw(df, group, args):

    pt = PowerTransformer()
    pt.fit(df)
    normcountst = pt.transform(df)

    pca = PCA(n_components=2)

    pca_f = pca.fit_transform(normcountst.T)
    pca_ve = np.around(pca.explained_variance_ratio_ * 100, 3)
    pca_df = pd.DataFrame(
        data=pca_f,
        index=df.T.index,
        columns=['PC1', 'PC2'],
    )

    conditions = pd.read_csv(args.groupsfile, header=0, index_col=0)

    if 'BATCH' in conditions.columns:
        group = 'BATCH'
    else:
        group = group

    a = pd.merge(pca_df, conditions, left_index=True, right_index=True)
    custom_palette = sns.mpl_palette("turbo", n_colors=len(a[group].unique()))
    sns.kdeplot(data=a, x="PC1", y="PC2", hue=group, fill=True, alpha=0.2, thresh=0.2, levels=3,
                palette=custom_palette, common_norm=False)
    sns.scatterplot(x='PC1', y='PC2', data=a, hue=group, palette=custom_palette)
    plt.xlabel(f"PC1({pca_ve[0]}%)")
    plt.ylabel(f"PC2({pca_ve[1]}%)")
    plt.title('PCA raw counts')
    plt.savefig(args.outputfolder / 'images' / 'pcaraw.png')
    plt.savefig(args.outputfolder / 'images' / 'pcaraw2.png', dpi=80)
    plt.close()

def plotpcanorm(df, group, args):

    pt = PowerTransformer()
    pt.fit(df)
    normcountst = pt.transform(df)

    pca = PCA(n_components=2)

    pca_f = pca.fit_transform(normcountst.T)
    pca_ve = np.around(pca.explained_variance_ratio_ * 100, 3)
    pca_df = pd.DataFrame(
        data=pca_f,
        index=df.T.index,
        columns=['PC1', 'PC2'],
    )

    conditions = pd.read_csv(args.groupsfile, header=0, index_col=0)

    a = pd.merge(pca_df, conditions, left_index=True, right_index=True)
    custom_palette = sns.mpl_palette("turbo", n_colors=len(a[group].unique()))
    sns.kdeplot(data=a, x="PC1", y="PC2", hue=group, fill=True, alpha=0.2, thresh=0.2, levels=3, palette=custom_palette, common_norm=False)
    sns.scatterplot(x='PC1', y='PC2', data=a, hue=group, palette=custom_palette)
    plt.xlabel(f"PC1({pca_ve[0]}%)")
    plt.ylabel(f"PC2({pca_ve[1]}%)")
    plt.title('PCA normalized counts')
    plt.savefig(args.outputfolder / 'images' / 'pcanorm.png')
    plt.savefig(args.outputfolder / 'images' / 'pcanorm2.png', dpi=80)
    plt.close()

def plotevalpcas(args):
    rawfcounts = pd.read_csv(args.outputfolder / 'otherfiles' / 'rawfcounts.csv', index_col=0)
    normcounts = pd.read_csv(args.outputfolder / 'results' / 'rnormcounts.csv', index_col=0)

    print('Plotting raw PCA plot...')
    plotpcaraw(rawfcounts, args.grouppca, args)
    logging.info('Plotted raw PCA plot')

    print('Plotting normalized PCA plot...')
    plotpcanorm(normcounts, args.grouppca, args)
    logging.info('Plotted normalized PCA plot')


def evalnorm(args):
    args.start_time = time.time()
    rawcounts = pd.read_csv(args.outputfolder / 'otherfiles' / "rawcounts2.csv",
                            index_col=0)

    rnormcounts = pd.read_csv(args.outputfolder / 'results' / 'rnormcounts.csv', index_col=0)

    rlegenes = RLEcal(logarizeoutput(rnormcounts, args), args)
    rleraw = RLEcal(logarizeoutput(rawcounts, args), args)

    meaniqr = getmeaniqr(rlegenes)
    meaniqrraw = getmeaniqr(rleraw)

    print('Plotting raw RLE plot...')
    plotevalraw(rawcounts, 'RAW counts', meaniqrraw, args)
    logging.info('Plotted raw plots')

    print('Plotting normalized RLE plot...')
    plotevalnorm(rnormcounts, 'Fully normalized counts', meaniqr, args)
    logging.info('Plotted normalized plots')

    plotevalpcas(args)

    pdfreportnorm(args)

    args.elapsed += time.time() - args.start_time
    args.current_state = '--> Finished. ' +\
        'Elapsed %s seconds ' + str(args.elapsed)
    logging.info(args.current_state)



    return (meaniqrraw, meaniqr)

if __name__ == '__main__':
    args = argParser()
    runQCview(args)
    runQCfilter(args)
    pipeline1(args)
    pipeline2(args)
    evalnorm(args)