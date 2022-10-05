import sys
import webbrowser
import time
from PyQt6.QtWidgets import (
    QMainWindow, QApplication, QWidget, QPushButton, QMessageBox, QComboBox, QFileDialog, QLCDNumber, QSlider, QSpinBox,
    QLabel, QToolBar, QStatusBar, QGridLayout, QLineEdit, QTextEdit, QDoubleSpinBox, QHBoxLayout, QVBoxLayout, QFormLayout, QCheckBox
)
from PyQt6.QtGui import QAction, QIcon, QPalette, QColor, QPixmap, QFont
from PyQt6.QtCore import Qt, QSize, QThread

import state
import guanin

class MainWindow(QMainWindow):

    def __init__(self):
        super(MainWindow, self).__init__()
        self.state = state.ConfigData()

        self.setWindowTitle("GUANIN: Nanostring Interactive Normalization")
        self.setWindowIcon(QIcon('image/logoguanin_156x156.png'))
        # self.resize(1080,640)


        self.statusbar = QStatusBar(self)
        self.statusbar.showMessage(self.state.current_state)
        self.setStatusBar(self.statusbar)


        widget = CentralWidget(self)
        self.setCentralWidget(widget)

        toolbar = QToolBar("Toolbar")
        toolbar.setIconSize(QSize(16,16))
        self.addToolBar(toolbar)

        button_action = QAction(QIcon("icons/database_plus.png"), "Load RCC data", self)
        button_action.setStatusTip("Click to select RCC input data")
        button_action.triggered.connect(self.onMyToolBarButtonClick)
        button_action.setCheckable(True)
        toolbar.addAction(button_action)


        menubar = self.menuBar()

        exitAct = QAction(QIcon('minus_exit.png'), '&Exit', self)
        exitAct.setShortcut('Ctrl+Q')
        exitAct.setStatusTip('Exit application')
        exitAct.triggered.connect(QApplication.instance().quit)

        loadRCCsAct = QAction(QIcon("icons/database_plus.png"), 'Load RCCs', self)
        loadRCCsAct.setShortcut('Ctrl+L')
        loadRCCsAct.setStatusTip('Select input RCC files')
        # somehow trigger el click con el abrir files loadRCCsAct.triggered.connect(QApplication.instance().file_open)

        aboutgenvipAct = QAction(QIcon("icons/genvip_24x24.png"), 'About GENVIP', self)
        # aboutgenvipAct.connect(webbrowser.open('http://www.genvip.eu'))

        aboutguaninAct = QAction(QIcon('icons/logoguanin_32x32.png'), 'About GUANIN', self)

        filemenu = menubar.addMenu('File')
        filemenu.addAction(exitAct)
        filemenu.addAction(loadRCCsAct)

        aboutmenu = menubar.addMenu('About')
        aboutmenu.addAction(aboutgenvipAct)
        aboutmenu.addAction(aboutguaninAct)

    def file_open(self):
        print('hay que meter que abra archivos')

    def closeEvent(self, event):

        reply = QMessageBox.question(self, 'Confirm exit',
                    "Are you sure to quit?", QMessageBox.StandardButton.Yes |
                    QMessageBox.StandardButton.No, QMessageBox.StandardButton.No)

        if reply == QMessageBox.StandardButton.Yes:

            event.accept()
        else:

            event.ignore()


    def onMyToolBarButtonClick(self, s):
        print("click", s)

class CentralWidget(QWidget):
    def __init__(self, parent):

        super().__init__()
        self.parent = parent
        self.parent.statusBar().showMessage('Ready to start analysis')
        self.state = state.ConfigData()

        self.initUI()

    def initUI(self):


        layout1 = QHBoxLayout()
        layout2 = QVBoxLayout()
        layout3 = QVBoxLayout()
        layload = QFormLayout()
        layqc = QFormLayout()
        laytnorm = QFormLayout()
        laycnorm = QFormLayout()
        layeval = QFormLayout()
        layexport = QFormLayout()

        layout1.setContentsMargins(10, 5, 5, 5)
        layout1.setSpacing(10)

        loadtitle = QLabel('- Loading data -')
        loadtitle.setAlignment(Qt.AlignmentFlag.AlignCenter)
        headerfont = QFont('Roboto', 11)
        headerfont.setBold(True)
        loadtitle.setFont(headerfont)
        layload.addRow(loadtitle)

        rccfolderbutton = QPushButton('Select folder containing RCC')
        rccfolderbutton.clicked.connect(self.openselectfolder)
        layload.addRow('Select RCC files', rccfolderbutton)

        self.showfoldertextbox = QLabel(str(self.state.folder))
        layload.addRow('Selected folder: ', self.showfoldertextbox)


        csvgroupsbutton = QPushButton('Select csv file')
        csvgroupsbutton.clicked.connect(self.opencsvfile)
        layload.addRow('Select groups file', csvgroupsbutton)

        self.showfiletextbox = QLabel(self.state.groupsfile)
        layload.addRow('Selected groups file: ', self.showfiletextbox)


        sampleidentifier = QComboBox()
        sampleidentifier.addItem('Sample ID')
        sampleidentifier.addItem('Filename')
        sampleidentifier.currentIndexChanged.connect(self.changingmodeid)

        layload.addRow('Sample identifier', sampleidentifier)

        doubleforloading = QHBoxLayout()
        doubleforticksloading = QHBoxLayout()
        doubleforloading.addLayout(doubleforticksloading)
        runloading = QPushButton('Run load RCCs')
        runloading.clicked.connect(self.runloadingrccs)

        runloading.setIcon(QIcon('image/logoguanin_96x96.png'))
        doubleforloading.addWidget(runloading)

        doubleformtic1 = QFormLayout()
        ticpopupinfolanes = QCheckBox()
        ticpopupinfolanes.stateChanged.connect(self.change_showbrowser)
        doubleformtic1.addRow('Pop up infolanes and rawQC report', ticpopupinfolanes)
        doubleforticksloading.addLayout(doubleformtic1)

        layload.addRow(doubleforloading)

        layout2.addLayout(layload)

        backgroundbutton = QComboBox()
        backgroundbutton.addItem('Mean+2std of neg ctrls')
        backgroundbutton.addItem('Max of neg controls')
        backgroundbutton.addItem('Mean of neg controls')
        backgroundbutton.addItem('Alternative background (filtered neg ctrls)')
        backgroundbutton.addItem('Manual background')
        backgroundbutton.currentIndexChanged.connect(self.changingbackgroundbutton)

        qctitle = QLabel('- Quality Control parameters -')
        qctitle.setAlignment(Qt.AlignmentFlag.AlignCenter)
        qctitle.setFont(headerfont)
        layqc.addRow(qctitle)
        layqc.addRow('Choose background', backgroundbutton)

        manualbackground = QSpinBox()
        manualbackground.valueChanged.connect(self.changingmanualbackground)
        layqc.addRow('Set manual background', manualbackground)

        backgroundcorrection = QComboBox()
        backgroundcorrection.addItem('Set as background')
        backgroundcorrection.addItem('Sustract background value')
        backgroundcorrection.addItem('Skip background correction')
        backgroundcorrection.currentIndexChanged.connect(self.changingbackgroundcorrection)

        layqc.addRow('Background correction (low counts)', backgroundcorrection)

        removelanes = QSpinBox()
        removelanes.setValue(80)
        removelanes.setMinimum(0)
        removelanes.setMaximum(100)
        removelanes.valueChanged.connect(self.changingpbelowbackground)
        layqc.addRow('% of low counts for lane remove', removelanes)


        doubleforfov = QHBoxLayout()
        minfovlay = QFormLayout()
        minfovline = QDoubleSpinBox(singleStep = 0.01)
        minfovline.setValue(0.75)
        minfovline.setMinimum(0)
        minfovline.setMaximum(10)
        minfovline.textChanged.connect(self.changeminfov)
        minfovlay.addRow('Min fov', minfovline)
        maxfovlay = QFormLayout()
        maxfovline = QDoubleSpinBox(singleStep = 0.01)
        maxfovline.setValue(1)
        maxfovline.setMinimum(0)
        maxfovline.setMaximum(10)
        maxfovline.textChanged.connect(self.changemaxfov)
        maxfovlay.addRow('Max fov', maxfovline)

        doubleforfov.addLayout(minfovlay)
        doubleforfov.addLayout(maxfovlay)
        layqc.addRow(doubleforfov)


        doubleforbd = QHBoxLayout()
        minbdlay = QFormLayout()
        minbdline = QDoubleSpinBox(singleStep = 0.01)
        minbdline.setValue(0.1)
        minbdline.setMinimum(0)
        minbdline.setMaximum(10)
        minbdline.textChanged.connect(self.changeminbd)
        minbdlay.addRow('Min binding density', minbdline)

        maxbdlay = QFormLayout()
        maxbdline = QDoubleSpinBox(singleStep=0.01)
        maxbdline.setValue(1.8)
        maxbdline.setMinimum(0)
        maxbdline.setMaximum(100)
        maxbdline.textChanged.connect(self.changemaxbd)
        maxbdlay.addRow('Max binding density', maxbdline)

        doubleforbd.addLayout(minbdlay)
        doubleforbd.addLayout(maxbdlay)
        layqc.addRow(doubleforbd)

        doubleforlin = QHBoxLayout()
        minlinlay = QFormLayout()
        minlinline = QDoubleSpinBox(singleStep = 0.01)
        minlinline.setValue(0.75)
        minlinline.setMinimum(0)
        minlinline.setMaximum(1)
        minlinline.textChanged.connect(self.changeminlin)
        minlinlay.addRow('Min linearity', minlinline)

        maxlinlay = QFormLayout()
        maxlinline = QDoubleSpinBox(singleStep = 0.01)
        maxlinline.setValue(1)
        maxlinline.setMinimum(0)
        maxlinline.setMaximum(1)
        maxlinline.textChanged.connect(self.changemaxlin)
        maxlinlay.addRow('Max linearity', maxlinline)

        doubleforlin.addLayout(minlinlay)
        doubleforlin.addLayout(maxlinlay)
        layqc.addRow(doubleforlin)

        doubleforscaf = QHBoxLayout()
        minscaflay = QFormLayout()
        minscafline = QDoubleSpinBox(singleStep = 0.01)
        minscafline.setValue(0.3)
        minscafline.setMinimum(0)
        minscafline.setMaximum(100)
        minscafline.textChanged.connect(self.changeminscaf)
        minscaflay.addRow('Min scaling factor', minscafline)

        maxscaflay = QFormLayout()
        maxscafline = QDoubleSpinBox(singleStep = 0.01)
        maxscafline.setValue(3)
        maxscafline.setMinimum(0)
        maxscafline.setMaximum(100)
        maxscafline.textChanged.connect(self.changemaxscaf)
        maxscaflay.addRow('Max scaling factor', maxscafline)

        doubleforscaf.addLayout(minscaflay)
        doubleforscaf.addLayout(maxscaflay)
        layqc.addRow(doubleforscaf)

        sampleremovercombobox = QComboBox()
        sampleremovercombobox.addItem('Remove auto-QC flagged')
        sampleremovercombobox.addItem('Keep all samples')
        sampleremovercombobox.addItem('Remove manually selected samples')
        sampleremovercombobox.currentIndexChanged.connect(self.changesampleremoving)
        layqc.addRow('Remove bad samples?', sampleremovercombobox)

        manualremoveselection = QLineEdit()
        manualremoveselection.textChanged.connect(self.changemanualremoveinput)
        layqc.addRow('Manual input lanes to remove:', manualremoveselection)

        doubleforqcfiltering = QHBoxLayout()
        doubleforticksqcfiltering = QHBoxLayout()
        doubleforqcfiltering.addLayout(doubleforticksqcfiltering)
        runqcfiltering = QPushButton('Run QC filtering')
        runqcfiltering.setIcon(QIcon('image/logoguanin_96x96.png'))
        runqcfiltering.clicked.connect(self.runqc)
        doubleforqcfiltering.addWidget(runqcfiltering)

        doubleformtic1qc = QFormLayout()
        doubleformtic1qc.addRow('Pop out new infolanes and QC report', QCheckBox())
        doubleforticksqcfiltering.addLayout(doubleformtic1qc)

        layqc.addRow(doubleforqcfiltering)

        layout2.addLayout(layqc)


        tnormtitle = QLabel('- Technical normalization parameters -')
        tnormtitle.setAlignment(Qt.AlignmentFlag.AlignCenter)
        tnormtitle.setFont(headerfont)
        laytnorm.addRow(tnormtitle)

        tnormmethodcombobox = QComboBox()
        tnormmethodcombobox.addItem('Use posgeomean')
        tnormmethodcombobox.addItem('Use summation')
        tnormmethodcombobox.addItem('Use median')
        tnormmethodcombobox.addItem('Use regression')
        tnormmethodcombobox.currentIndexChanged.connect(self.changetnormmethod)
        laytnorm.addRow('Method for technical normalization: ', tnormmethodcombobox)

        doublefortechnorm = QHBoxLayout()
        doublefortickstechnorm = QHBoxLayout()
        doublefortechnorm.addLayout(doublefortickstechnorm)
        runtechnorm = QPushButton('Run technical normalization')
        runtechnorm.setIcon(QIcon('image/logoguanin_96x96.png'))
        runtechnorm.clicked.connect(self.runthetechnorm)
        doublefortechnorm.addWidget(runtechnorm)

        doubleformtic1tn = QFormLayout()
        ticpopupinfolanes2 = QCheckBox()
        ticpopupinfolanes2.stateChanged.connect(self.change_showbrowser)
        doubleformtic1tn.addRow('Pop out infolanes and QC report', ticpopupinfolanes2)

        doublefortickstechnorm.addLayout(doubleformtic1tn)

        laytnorm.addRow(doublefortechnorm)

        layout2.addLayout(laytnorm)

        layout1.addLayout(layout2)


        laycnorm = QFormLayout()
        cnormtitle = QLabel('- Content normalization parameters -')
        cnormtitle.setAlignment(Qt.AlignmentFlag.AlignCenter)
        cnormtitle.setFont(headerfont)
        laycnorm.addRow(cnormtitle)

        filhousekeepingsmincounts = QSpinBox()
        filhousekeepingsmincounts.setValue(50)
        filhousekeepingsmincounts.textChanged.connect(self.change_filhousekeepingmincounts)
        laycnorm.addRow('Filter housekeeping panel genes by min counts', filhousekeepingsmincounts)

        includeerg = QCheckBox()
        includeerg.stateChanged.connect(self.change_includeerg)
        laycnorm.addRow('Include best endogenous as reference candidates', includeerg)

        howmanyergs = QSpinBox()
        howmanyergs.setValue(6)
        howmanyergs.textChanged.connect(self.change_howmanyergs)
        laycnorm.addRow('How many best endogenous genes to include', howmanyergs)

        whatrefgenestouse = QComboBox()
        whatrefgenestouse.addItem('Genorm auto selection (default)')
        whatrefgenestouse.addItem('Top best n genes from genorm ranking')
        whatrefgenestouse.addItem('Top n most expressed endogenous')
        whatrefgenestouse.addItem('All endogenous genes (may be slow)')
        whatrefgenestouse.addItem('Ponderated genorm selection')
        whatrefgenestouse.addItem('Manual selection of genes')
        whatrefgenestouse.currentIndexChanged.connect(self.change_contnormmethod)

        laycnorm.addRow('What reference genes selection to use?', whatrefgenestouse)

        inputngenes = QSpinBox()
        inputngenes.setValue(6)
        inputngenes.setMaximum(50000)
        inputngenes.textChanged.connect(self.change_nrefgenes)
        laycnorm.addRow('If n genes to be set from last option: ', inputngenes)

        inputrefgenesline = QLineEdit()
        inputrefgenesline.textChanged.connect(self.change_inputnamesrefgenes)
        laycnorm.addRow('If manual selection is chosen, input genes:', inputrefgenesline)

        howtofilter = QComboBox()
        howtofilter.addItem('Kruskal-Wallis filtering')
        howtofilter.addItem('Kruskal-Wallis flagging')
        howtofilter.addItem('Wilcoxon filtering')
        howtofilter.addItem('Wilcoxon flagging')
        howtofilter.addItem('No group filtering/flagging')
        howtofilter.currentIndexChanged.connect(self.change_howtofilter)
        laycnorm.addRow('How to filter/flag reference genes', howtofilter)

        showreversefeatureselectionranking = QCheckBox()
        showreversefeatureselectionranking.stateChanged.connect(self.change_showbrowser)
        laycnorm.addRow('Show reverse feature selection ranking?', showreversefeatureselectionranking)

        additionalnormalization = QComboBox()
        additionalnormalization.addItem('No')
        additionalnormalization.addItem('Quantile normalization')
        additionalnormalization.addItem('Standarization')


        additionalnormalization.currentIndexChanged.connect(self.change_adnormalization)

        laycnorm.addRow('Perform additional normalization?', additionalnormalization)


        howtoexport = QComboBox()
        howtoexport.addItem('Normalized count matrix')
        howtoexport.addItem('log2 count matrix')
        howtoexport.addItem('log10 count matrix')
        howtoexport.currentIndexChanged.connect(self.change_exportmethod)
        laycnorm.addRow('Format export results', howtoexport)


        doubleforcnorm = QHBoxLayout()
        doublefortickscnorm = QHBoxLayout()
        doubleforcnorm.addLayout(doublefortickscnorm)
        runcnorm = QPushButton('Run content normalization')
        runcnorm.setIcon(QIcon('image/logoguanin_96x96.png'))
        runcnorm.clicked.connect(self.runcnorm)
        doubleforcnorm.addWidget(runcnorm)

        doubleformtic1cn = QFormLayout()
        ticpopupinfolanes3 = QCheckBox()
        ticpopupinfolanes3.stateChanged.connect(self.change_showbrowser)
        doubleformtic1cn.addRow('Pop out infolanes and QC report', ticpopupinfolanes3)
        doublefortickscnorm.addLayout(doubleformtic1cn)

        laycnorm.addRow(doubleforcnorm)

        layout3.addLayout(laycnorm)

        evaltitle = QLabel('- Normalization evaluation -')
        evaltitle.setAlignment(Qt.AlignmentFlag.AlignCenter)
        evaltitle.setFont(headerfont)
        layeval.addRow(evaltitle)


        doubleforeval = QHBoxLayout()
        doublefortickseval = QHBoxLayout()
        doubleforeval.addLayout(doublefortickseval)
        runevalbutton = QPushButton('Run evaluation')
        runevalbutton.setIcon(QIcon('image/logoguanin_96x96.png'))
        runevalbutton.clicked.connect(self.runeval)
        doubleforeval.addWidget(runevalbutton)

        doubleformtic1ev = QFormLayout()
        ticpopupinfolanes4 = QCheckBox()
        ticpopupinfolanes4.stateChanged.connect(self.change_showbrowser)
        doubleformtic1ev.addRow('Pop out infolanes and final report', ticpopupinfolanes4)

        doublefortickseval.addLayout(doubleformtic1ev)

        # label1 = QLabel(self)
        # pixmap = QPixmap('/home/juli/PycharmProjects/guanin/output/images/rlenormplot2.png')
        # label1.setPixmap(pixmap)
        # label1.setScaledContents(True)

        layeval.addRow(doubleforeval)

        layout3.addLayout(layeval)

        layout1.addLayout(layout3)

        self.setLayout(layout1)



    def openselectfolder(self):
        folder = QFileDialog.getExistingDirectory(self)
        self.state.folder = folder
        self.showfoldertextbox.setText(self.state.folder)

    def opencsvfile(self):
        file = QFileDialog.getOpenFileName(self)
        self.state.groupsfile = file[0]
        self.showfiletextbox.setText(self.state.groupsfile)

    def changingmodeid(self, value):
        if value == 1:
            self.state.modeid = 'filename'
        elif value == 0:
            self.state.modeid = 'sampleID'
        print(self.state.modeid)

    def change_showbrowser(self, checkbox):
        if checkbox == 0:
            self.state.showbrowser = False
        elif checkbox == 2:
            self.state.showbrowser = True
        print(self.state.showbrowser)

    def runloadingrccs(self):
        self.parent.statusBar().showMessage('Loading RCC files...')
        self.parent.statusBar().repaint()
        guanin.runQCview(self.state)
        self.parent.statusBar().showMessage('RCC files loaded, ready to perform QC')


    def changingbackgroundbutton(self, checkbox):
        if checkbox == 0:
            self.state.background = 'Background'
            print(self.state.background)
        elif checkbox == 1:
            self.state.background = 'Background2'
            print(self.state.background)
        elif checkbox == 2:
            self.state.background = 'Background3'
            print(self.state.background)
        elif checkbox == 3:
            self.state.background = 'Backgroundalt'
            print(self.state.background)
        elif checkbox == 4:
            self.state.background = 'Backgroundalt'
            print(self.state.background)

    def changingmanualbackground(self, value):
        self.state.manualbackground = value
        print(self.state.manualbackground)

    def changingbackgroundcorrection(self, checkbox):
        if checkbox == 0:
            self.state.lowcounts = 'asim'
            print(self.state.lowcounts)
        elif checkbox == 1:
            self.state.lowcounts = 'sustract'
            print(self.state.lowcounts)
        elif checkbox == 2:
            self.state.lowcounts = 'skip'
            print(self.state.lowcounts)

    def changingpbelowbackground(self, value):
        self.state.pbelowbackground = value
        print(self.state.pbelowbackground)

    def changeminfov(self, value):
        self.state.minfov = value
        print(self.state.minfov)

    def changemaxfov(self, value):
        self.state.maxfov = value
        print(self.state.maxfov)

    def changeminbd(self, value):
        self.state.minbd = value
        print(self.state.minbd)

    def changemaxbd(self, value):
        self.state.maxbd = value
        print(self.state.maxbd)

    def changeminlin(self, value):
        self.state.minlin = value
        print(self.state.minlin)

    def changemaxlin(self, value):
        self.state.maxlin = value
        print(self.state.maxlin)

    def changeminscaf(self, value):
        self.state.minscaf = value
        print(self.state.minscaf)

    def changemaxscaf(self, value):
        self.state.maxscaf = value
        print(self.state.maxscaf)

    def changesampleremoving(self, checkbox):
        if checkbox == 2:
            self.state.remove = 'variable de luego manual remove'
            print(self.state.remove)
        elif checkbox == 1:
            self.state.laneremover = 'no'
            print(self.state.laneremover + ' lanes removed')
        if checkbox == 0:
            self.state.laneremover = 'yes'
            self.state.remove = None
            print('qc remove')

    def changemanualremoveinput(self, value):
        self.state.remove = value
        print(self.state.remove)

    def runqc(self):
        self.parent.statusBar().showMessage('Aplying filters and QC...')
        self.parent.statusBar().repaint()
        guanin.runQCfilter(self.state)
        self.parent.statusBar().showMessage('QC done, ready to perform technical normalization')


    def changetnormmethod(self, checkbox):
        if checkbox == 0:
            self.state.tecnormeth = 'posgeomean'
            print(self.state.tecnormeth)
        elif checkbox == 1:
            self.state.tecnormeth = 'Sum'
            print(self.state.tecnormeth)
        elif checkbox == 2:
            self.state.tecnormeth = 'Median'
            print(self.state.tecnormeth)
        elif checkbox == 3:
            self.state.tecnormeth = 'regression'
            print(self.state.tecnormeth)

    def runthetechnorm(self):
        self.parent.statusBar().showMessage('Performing technical normalization')
        self.parent.statusBar().repaint()
        guanin.technorm(self.state)
        self.parent.statusBar().showMessage('Technical normalization done, ready to perform content normalization')


    def change_filhousekeepingmincounts(self, value):
        self.state.mincounthkes = value
        print(self.state.mincounthkes)

    def change_includeerg(self, checkbox):
        if checkbox == 0:
            self.state.refendgenes = 'hkes'
            print(self.state.refendgenes)
        elif checkbox == 2:
            self.state.refendgenes = 'endhkes'
            print(self.state.refendgenes)

    def printconf(self):
        print(self.state)

    def change_howmanyergs(self, value):
        self.state.numend = value
        print(self.state.numend)

    def change_contnormmethod(self, checkbox):
        if checkbox == 0:
            self.state.contnorm = 'refgenes'
            print(self.state.contnorm)
        elif checkbox == 1:
            self.state.contnorm = 'refgenes'
            print(self.state.contnorm)
        elif checkbox == 2:
            self.state.contnorm = 'topn'
            print(self.state.contnorm)
        elif checkbox == 3:
            self.state.contnorm = 'all'
            print(self.state.contnorm)
        elif checkbox == 4:
            self.state.contnorm = 'ponderaterefgenes'
            print(self.state.contnorm)
        elif checkbox == 5:
            self.state.contnorm = 'refgenes'
            print(self.state.contnorm)

    def change_nrefgenes(self, value):
        if self.state.contnorm == 'refgenes':
            self.state.nrefgenes = value
            print(self.state.nrefgenes)
        elif self.state.contnorm == 'ponderaterefgenes':
            self.state.nrefgenes = value
            print(self.state.nrefgenes)
        elif self.state.contnorm == 'topn':
            self.state.topngenestocontnorm = value
            print(self.state.topngenestocontnorm)

    def change_inputnamesrefgenes(self, value):
        self.state.chooserefgenes = value
        print(self.state.chooserefgenes)

    def change_howtofilter(self, checkbox):
        if checkbox == 0:
            self.state.filtergroupvariation = 'filterkrus'
        elif checkbox == 1:
            self.state.filtergroupvariation = 'flagkurs'
        elif checkbox == 2:
            self.state.filtergroupvariation = 'filterwilcox'
        elif checkbox == 3:
            self.state.filtergroupvariation = 'flagwilcox'
        elif checkbox == 4:
            self.state.filtergroupvariation = 'nofilter'
        print(self.state.filtergroupvariation)

    def change_adnormalization(self, checkbox):
        if checkbox == 0:
            self.state.adnormalization = 'no'
        elif checkbox == 1:
            self.state.adnormalization = 'quantile'
        elif checkbox == 2:
            self.state.adnormalization = 'standarization'
        print(self.state.adnormalization)

    def runcnorm(self):
        self.parent.statusBar().showMessage('Performing content normalization...')
        self.parent.statusBar().repaint()
        guanin.contnorm(self.state)
        self.parent.statusBar().showMessage('Content normalization done, ready to evaluate normalization')

    def runeval(self):
        self.parent.statusBar().showMessage('Performing evaluation, plotting RLE...')
        self.parent.statusBar().repaint()
        guanin.evalnorm(self.state)
        self.parent.statusBar().showMessage('Evaluation and data export ready, check "output" folder')


    def change_exportmethod(self, checkbox):
        if checkbox == 0:
            self.state.logarizedoutput = 'no'
        elif checkbox == 1:
            self.state.logarizedoutput = '2'
        elif checkbox == 2:
            self.state.logarizedoutput = '10'
        print(self.state.logarizedoutput)


def main():
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    app.exec()

if __name__ == '__main__':
    main()