import time
import pathlib
import tempfile

class ConfigData:
    def __init__(self):
        self.folder = pathlib.Path(__file__).parent/'examples/d1_COV_GSE183071'
        self.minfov = 0.75
        self.maxfov = 1
        self.minbd = 0.1
        self.maxbd = 1.8
        self.minlin = 0.75
        self.maxlin = 1
        self.minscalingfactor = 0.3
        self.maxscalingfactor = 3
        self.showbrowserrawqc = False
        self.showbrowserqc = False
        self.showbrowsercnorm = False
        self.lowcounts = 'sustract'
        self.modeid = 'filename'
        self.modeview = 'view'
        self.tecnormeth = 'posgeomean'
        self.refendgenes = 'endhkes'
        self.remove = None
        self.background = 'Background'
        self.pbelowbackground = 80
        self.manualbackground = None
        self.chooserefgenes = None
        self.filtergroupvariation = 'filterkrus'
        self.featureselectionneighbors = 4
        self.groups = 'no'
        self.numend = 9
        self.autorename = 'off'
        self.contnorm = 'refgenes'
        self.adnormalization = 'no'
        self.topngenestocontnorm = '100'
        self.mincounthkes = 80
        self.nrefgenes = None
        self.laneremover = 'yes'
        self.groupsinrnormgenes = 'no'
        self.logarizedoutput = '10'
        self.logarizeforeval = '10'
        self.groupsfile = pathlib.Path(__file__).parent/'groups_d1_COV_GSE183071.csv'
        self.start_time = time.time()
        self.current_state = 'Ready to analysis'
        self.current_state = 'Ready to analysis'
        self.badlanes = 'No bad lanes detected'
        self.rawmeaniqr = 'Raw IQR not calculated yet'
        self.normmeaniqr = 'Norm IQR not calculated yet'
        self.firsttransformlowcounts = True
        self.outputfolder = tempfile.gettempdir() + '/guanin_output'
        self.showlastlog = False
        self.refgenessel = ''

    def __str__(self):
        return self.modeid