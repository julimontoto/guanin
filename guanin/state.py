import logging
from pathlib import Path
import tempfile
import time


class ConfigData:
    def __init__(self, *args, **kwargs):
        self.folder = Path(__file__).parent / "examples" / "d1_COV_GSE183071"
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
        self.lowcounts = 'subtract'
        self.modeid = 'filename'
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
        self.groups = False
        self.numend = 9
        self.autorename = False
        self.contnorm = 'refgenes'
        self.adnormalization = 'no'
        self.topngenestocontnorm = '100'
        self.mincounthkes = 80
        self.nrefgenes = None
        self.laneremover = True
        self.groupsinrnormgenes = False
        self.logarizedoutput = '10'
        self.logarizeforeval = '10'
        self.groupsfile = Path(__file__).parent / "examples/groups_d1_COV_GSE183071.csv"
        self.start_time = time.time()
        self.current_state = 'Ready to analysis'
        self.current_state = 'Ready to analysis'
        self.nbadlanes = 'No bad lanes detected'
        self.rawmeaniqr = 'Raw IQR not calculated yet'
        self.normmeaniqr = 'Norm IQR not calculated yet'
        self.tnormbeforebackgcorr = True
        self.outputfolder = kwargs.get(
            "output_folder",
            Path(tempfile.gettempdir()) / "guanin_output")
        self.showlastlog = False
        self.refgenessel = []
        self.indexmethod2 = 0
        self.kvalue = 3
        self.pipeline = 'ruvgnorm'
        self.whatrefgenes = []
        self.eme = None
        self.grouppca = 'GROUP'
        self.deseq2_mor = True
        self.manual_remove = False
        self.badlanes = set()
        self.elapsed = 0.0

    def change_float(self, name, value):
        attr = getattr(self, name, None)
        if attr:
            attr = float(value.replace(",", "."))
        logging.info(f"state.{name} = {attr}")

    def __str__(self):
        return self.modeid
