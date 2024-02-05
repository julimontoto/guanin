import argparse
import logging
import os
import pathlib
import sys
import webbrowser
try:
    from PyQt6.QtWidgets import (
        QMainWindow, QApplication, QWidget, QPushButton, QMessageBox,
        QComboBox, QFileDialog, QSpinBox, QSplashScreen, QLabel, QStatusBar,
        QLineEdit, QDoubleSpinBox, QHBoxLayout, QVBoxLayout, QFormLayout,
        QCheckBox, QPlainTextEdit, QGridLayout, QTabWidget)
    from PyQt6.QtGui import QAction, QIcon, QPixmap, QFont, QGuiApplication
    from PyQt6.QtCore import Qt, QTimer
except ImportError:
    os.environ['LD_LIBRARY_PATH'] = \
        str(pathlib.Path(__file__).parent.parent / 'libraries')

    from PyQt6.QtWidgets import (
        QMainWindow, QApplication, QWidget, QPushButton, QMessageBox,
        QComboBox, QFileDialog, QSpinBox, QSplashScreen, QLabel, QStatusBar,
        QLineEdit, QDoubleSpinBox, QHBoxLayout, QVBoxLayout, QFormLayout,
        QCheckBox, QPlainTextEdit, QGridLayout, QTabWidget)
    from PyQt6.QtGui import QAction, QIcon, QPixmap, QFont, QGuiApplication
    from PyQt6.QtCore import Qt, QTimer


try:
    import guanin.guanin as guanin
    import guanin.state as state
except ModuleNotFoundError:
    import guanin
    import state


class MainWindow(QMainWindow):

    def __init__(self, *args, **kwargs):
        super(MainWindow, self).__init__()
        self.state = state.ConfigData(**kwargs.get("config", {}))
        app_path = pathlib.Path(__file__).parent
        icons_path = app_path / "icons"

        self.setWindowTitle("GUANIN: Nanostring Interactive Normalization")
        self.setWindowIcon(QIcon(
            str(app_path / "image" / "logoguanin_156x156.png")))
        # self.resize(1080,640)

        qtRectangle = self.frameGeometry()
        centerPoint = QGuiApplication.primaryScreen().availableGeometry().center()
        qtRectangle.moveCenter(centerPoint)
        self.move(qtRectangle.topLeft())
        # self.setFixedHeight(820)

        self.statusbar = QStatusBar(self)
        self.statusbar.showMessage(self.state.current_state)
        self.setStatusBar(self.statusbar)

        widget = CentralWidget(self)
        self.setCentralWidget(widget)

        menubar = self.menuBar()

        exitAct = QAction(QIcon(
            str(icons_path / "minus_exit.png")), '&Exit', self)
        exitAct.setShortcut('Ctrl+Q')
        exitAct.setStatusTip('Exit application')
        exitAct.triggered.connect(QApplication.instance().quit)

        aboutgenvipAct = QAction(QIcon(
            str(icons_path / "genvip_24x24.png")), 'About GENVIP', self)
        aboutgenvipAct.triggered.connect(self.popupgenvipweb)

        aboutgenpobTeam = QAction(QIcon(
            str(icons_path / "GenPob_logo.resized.png")),
            'About GENPOB Team', self)
        aboutgenpobTeam.triggered.connect(self.popupgenpobteam)

        aboutguaninAct = QAction(QIcon(
            str(icons_path / "logoguanin_32x32.png")), 'About GUANIN', self)
        aboutguaninAct.triggered.connect(self.popupguaningithub)

        aboutcitationAct = QAction(QIcon(
            str(icons_path / "book-open-bookmark.png")), 'Please cite', self)
        aboutcitationAct.triggered.connect(self.popupguaninpaper)

        aboutlicenseAct = QAction(QIcon(
            str(icons_path / "license-key.png")), 'GPL3 license', self)
        aboutlicenseAct.triggered.connect(self.popuplicenseinfo)

        viewlogAct = QAction(QIcon(
            str(icons_path / "information.png")), 'Analysis information', self)
        viewlogAct.triggered.connect(self.viewlog)

        viewpdfreportAct = QAction(QIcon(
            str(icons_path / "report-paper.png")), 'PDF report', self)
        viewpdfreportAct.triggered.connect(self.viewpdfreport)

        viewsummaryandinfolanesAct = QAction(QIcon(
            str(icons_path / "application-table.png")),
            'Summary and infolanes', self)
        viewsummaryandinfolanesAct.triggered.connect(self.viewsummaryandinfolanes)

        filemenu = menubar.addMenu('Guanin')
        viewmenu = menubar.addMenu('View')
        aboutmenu = menubar.addMenu('About')

        aboutmenu.addAction(aboutguaninAct)
        aboutmenu.addAction(aboutcitationAct)
        aboutmenu.addAction(aboutlicenseAct)
        filemenu.addAction(exitAct)

        viewmenu.addAction(viewlogAct)
        viewmenu.addAction(viewpdfreportAct)
        viewmenu.addAction(viewsummaryandinfolanesAct)

        aboutmenu.addAction(aboutgenvipAct)
        aboutmenu.addAction(aboutgenpobTeam)

    def viewlog(self):
        webbrowser.open(
            str(self.state.outputfolder / "info" / "analysis_description.log"))

    def viewpdfreport(self):
        webbrowser.open(
            str(self.state.outputfolder / "reports" / "QC_inspection.pdf"))
        webbrowser.open(
            str(self.state.outputfolder / "reports" / "norm_report.pdf"))

    def popuplicenseinfo(self):
        webbrowser.open('https://www.gnu.org/licenses/gpl-3.0.html')

    def viewsummaryandinfolanes(self):
        webbrowser.open(str(self.state.outputfolder / "info" / "rawsummary.html"))

    def popupgenpobteam(self):
        webbrowser.open('https://genpob.eu/team/')

    def popupguaninpaper(self):
        webbrowser.open('https://www.google.com/search?q=paper+pending+publication')

    def popupguaningithub(self):
        webbrowser.open('https://github.com/julimontoto/guanin')

    def popupgenvipweb(self):
        webbrowser.open('http://www.genvip.eu')

    def closeEvent(self, event):
        reply = QMessageBox.question(
            self,
            'Confirm exit',
            "Are you sure to quit?",
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
            QMessageBox.StandardButton.No)

        if reply == QMessageBox.StandardButton.Yes:
            logging.info('GUANIN session closed')

            event.accept()
        else:
            event.ignore()


class CentralWidget(QWidget):
    def __init__(self, parent):
        super().__init__()
        self.parent = parent
        self.parent.statusBar().showMessage('Ready to start analysis')
        self.state = parent.state

        self.initUI()

    def initUI(self):
        imgs_path = pathlib.Path(__file__).parent / "image"
        layout = QGridLayout()
        self.setLayout(layout)

        self.tabs = QTabWidget()

        layout.setContentsMargins(10, 5, 5, 5)
        layout.setSpacing(10)

        layload = QFormLayout()
        loadtitle = QLabel('- Loading data -')
        loadtitle.setAlignment(Qt.AlignmentFlag.AlignCenter)
        headerfont = QFont('Helvetica Neue', 11)
        headerfont.setBold(True)
        loadtitle.setFont(headerfont)
        layload.addRow(loadtitle)


        rccfolderbutton = QPushButton('Select folder containing RCC')
        rccfolderbutton.clicked.connect(self.openselectfolder)
        layload.addRow('RCC files location folder', rccfolderbutton)

        self.showfoldertextbox = QLabel(str(self.state.folder))
        layload.addRow('Selected input folder: ', self.showfoldertextbox)


        csvgroupsbutton = QPushButton('Select csv file')
        csvgroupsbutton.clicked.connect(self.opencsvfile)
        layload.addRow('Groups file location', csvgroupsbutton)

        self.showfiletextbox = QLabel(str(self.state.groupsfile))
        layload.addRow('Selected groups file: ', self.showfiletextbox)

        sampleidentifier = QComboBox()
        sampleidentifier.addItem('Filename')
        sampleidentifier.addItem('Sample ID')
        sampleidentifier.currentIndexChanged.connect(self.changingmodeid)

        layload.addRow('Sample identifier', sampleidentifier)

        outputfolderbutton = QPushButton('Select folder to generate output')
        outputfolderbutton.clicked.connect(self.openselectoutputfolder)
        layload.addRow('Output folder location', outputfolderbutton)

        self.showoutputfoldertextbox = QLabel(str(self.state.outputfolder))

        layload.addRow('Selected output folder: ', self.showoutputfoldertextbox)

        doubleforloading = QHBoxLayout()
        doubleforticksloading = QHBoxLayout()
        doubleforloading.addLayout(doubleforticksloading)
        runloading = QPushButton('Run load RCCs')
        runloading.clicked.connect(self.runloadingrccs)

        runloading.setIcon(QIcon(str(imgs_path / "logoguanin_96x96.png")))
        doubleforloading.addWidget(runloading)

        doubleformtic1 = QFormLayout()
        ticpopupinfolanes = QCheckBox()
        ticpopupinfolanes.stateChanged.connect(self.change_showbrowserrawqc)
        doubleformtic1.addRow('Pop up infolanes and rawQC report', ticpopupinfolanes)
        doubleforticksloading.addLayout(doubleformtic1)

        layload.addRow(doubleforloading)

        bigtabload = QWidget()
        bigtabload.setLayout(layload)

        self.tabs.addTab(bigtabload, QIcon(str(imgs_path / "logoguanin_96x96.png")), 'Loading data')

        layqc = QFormLayout()
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
        layqc.addRow('Background calculation method', backgroundbutton)

        manualbackground = QSpinBox()
        manualbackground.valueChanged.connect(self.changingmanualbackground)
        layqc.addRow('Set manual background', manualbackground)

        backgroundcorrection = QComboBox()
        backgroundcorrection.addItem('Subtract background value')
        backgroundcorrection.addItem('Set as background')
        backgroundcorrection.addItem('Skip background correction')
        backgroundcorrection.currentIndexChanged.connect(self.changingbackgroundcorrection)

        layqc.addRow('Background correction method (low counts)', backgroundcorrection)

        removelanes = QSpinBox()
        removelanes.setValue(80)
        removelanes.setMinimum(0)
        removelanes.setMaximum(100)
        removelanes.valueChanged.connect(self.changingpbelowbackground)
        layqc.addRow('% of low counts for lane remove', removelanes)


        doubleforfov = QHBoxLayout()
        minfovlay = QFormLayout()
        minfovline = QDoubleSpinBox(singleStep = 0.01, )
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
        sampleremovercombobox.addItem('Keep all samples, only flag bad samples')
        sampleremovercombobox.addItem('Remove manually selected samples')
        sampleremovercombobox.currentIndexChanged.connect(self.changesampleremoving)
        layqc.addRow('Method to remove low QC samples ', sampleremovercombobox)

        manualremoveselection = QLineEdit()
        manualremoveselection.textChanged.connect(self.changemanualremoveinput)
        layqc.addRow('Manual input lanes to remove:', manualremoveselection)

        doubleforqcfiltering = QHBoxLayout()
        doubleforticksqcfiltering = QHBoxLayout()
        doubleforqcfiltering.addLayout(doubleforticksqcfiltering)
        runqcfiltering = QPushButton('Run QC filtering')
        runqcfiltering.setIcon(QIcon(
            str(imgs_path / "logoguanin_96x96.png")))
        runqcfiltering.clicked.connect(self.runqc)
        runqcfiltering.clicked.connect(self.showflaggedlanes)
        doubleforqcfiltering.addWidget(runqcfiltering)

        doubleformtic1qc = QFormLayout()
        popoutinfolanescheck = QCheckBox()
        popoutinfolanescheck.stateChanged.connect(self.change_showbrowserqc)
        doubleformtic1qc.addRow('Pop out new infolanes and QC report',
                                popoutinfolanescheck)
        doubleforticksqcfiltering.addLayout(doubleformtic1qc)

        layqc.addRow(doubleforqcfiltering)

        self.showingflaggedlanes = QLabel(self.state.nbadlanes)
        layqc.addRow('Flagged lanes: ', self.showingflaggedlanes)

        bigtabqc = QWidget()
        bigtabqc.setLayout(layqc)

        self.tabs.addTab(bigtabqc, QIcon(str(imgs_path / "logoguanin_96x96.png")), 'Quality Control')

        self.laymethod = QFormLayout()
        methodtitle = QLabel('- Normalization method -')
        methodtitle.setAlignment(Qt.AlignmentFlag.AlignCenter)
        methodtitle.setFont(headerfont)
        self.laymethod.addRow(methodtitle)

        self.normethodcombobox = QComboBox()
        self.normethodcombobox.addItem('RUVg (remove unwanted variation through control genes)')
        self.normethodcombobox.addItem('Scaling factors (technical and content normalization)')
        self.normethodcombobox.currentIndexChanged.connect(self.changetnormmethod)
        self.laymethod.addRow('Method for normalization: ', self.normethodcombobox)

        self.list_scaf_options = ['Use posgeomean as scaling factor', 'Use summation as scaling factor', 'Use median as scaling factor', 'Use regression as scaling factor']
        self.list_scaf_tecnormeth = ['posgeomean', 'Sum', 'Median', 'regression']
        self.method2combobox = QComboBox()
        self.method2combobox.addItems(self.list_scaf_options)

        self.method2combobox.currentIndexChanged.connect(self.change_scaling_factor)

        self.inputkvalue = QSpinBox()
        self.inputkvalue.setValue(3)
        self.inputkvalue.valueChanged.connect(self.changekvalue)

        self.laymethod.addRow('k value', self.inputkvalue)

        self.list_text_tick = ['Perform median of ratios pre-normalization (size factors)', 'Perform technical normalization before background correction']

        # self.simpleformethod2combobox = QGridLayout()
        # self.simpleformethod2combobox.addWidget(self.method2combobox)
        #
        self.doublefortechnorm = QHBoxLayout()
        # self.doublefortickstechnorm = QHBoxLayout()
        # self.doublefortechnorm.addLayout(self.doublefortickstechnorm)

        self.doubleformtic1tn = QFormLayout()
        self.ticforaftertransformlowcounts = QCheckBox()
        # self.ticforaftertransformlowcounts.stateChanged.connect(
        #     self.changeaftertransformlowcounts)
        self.ticforaftertransformlowcounts.setChecked(False)

        # self.row_with_tick = [self.list_text_tick[int(self.state.tnormbeforebackgcorr)], self.ticforaftertransformlowcounts]
        self.doubleformtic1tn.addRow(self.list_text_tick[0], self.ticforaftertransformlowcounts)
        self.doubleformtic1tn.setRowVisible(0, True)

        self.ticforaftertransformlowcounts.stateChanged.connect(self.change_ticdeseq2_mor)

        self.doublefortechnorm.addLayout(self.doubleformtic1tn)

        runtechnorm = QPushButton('Continue normalization')
        runtechnorm.setIcon(QIcon(str(imgs_path / "logoguanin_96x96.png")))
        runtechnorm.clicked.connect(self.runthetechnorm)

        self.doublefortechnorm.addWidget(runtechnorm)

        self.laymethod.addRow(self.doublefortechnorm)

        bigtabmethod = QWidget()
        bigtabmethod.setLayout(self.laymethod)

        self.tabs.addTab(bigtabmethod, QIcon(str(imgs_path / "logoguanin_96x96.png")), 'Normalization method')

        self.laycnorm = QFormLayout()
        cnormtitle = QLabel('- Selection of reference genes parameters -')
        cnormtitle.setAlignment(Qt.AlignmentFlag.AlignCenter)
        cnormtitle.setFont(headerfont)
        self.laycnorm.addRow(cnormtitle)

        filhousekeepingsmincounts = QSpinBox()
        filhousekeepingsmincounts.setValue(50)
        filhousekeepingsmincounts.setMaximum(1000)
        filhousekeepingsmincounts.textChanged.connect(
            self.change_filhousekeepingmincounts)
        self.laycnorm.addRow('Min counts for housekeeping genes (all lanes) to not be filtered out of the analysis',
                        filhousekeepingsmincounts)

        includeerg = QCheckBox()
        includeerg.setChecked(True)
        includeerg.stateChanged.connect(self.change_includeerg)
        self.laycnorm.addRow('Include best endogenous as reference candidates',
                        includeerg)

        howmanyergs = QSpinBox()
        howmanyergs.setValue(6)
        howmanyergs.setMaximum(999)
        howmanyergs.textChanged.connect(self.change_howmanyergs)
        self.laycnorm.addRow('Number of top endogenous genes to include', howmanyergs)

        whatrefgenestouse = QComboBox()
        whatrefgenestouse.addItem('Genorm auto selection (default)')
        whatrefgenestouse.addItem('Top best n genes from genorm ranking')
        whatrefgenestouse.addItem('Top n most expressed endogenous')
        whatrefgenestouse.addItem('All endogenous genes')
        whatrefgenestouse.addItem('Ponderated genorm selection')
        whatrefgenestouse.addItem('Manual selection of genes')
        whatrefgenestouse.currentIndexChanged.connect(self.change_contnormmethod)

        self.laycnorm.addRow('Reference genes selection method', whatrefgenestouse)

        inputngenes = QSpinBox()
        inputngenes.setValue(6)
        inputngenes.setMaximum(50000)
        inputngenes.textChanged.connect(self.change_nrefgenes)
        self.laycnorm.addRow('Number of reference genes (if  is selected from last option)', inputngenes)

        inputrefgenesline = QLineEdit()
        inputrefgenesline.textChanged.connect(self.change_inputnamesrefgenes)
        self.laycnorm.addRow('If manual selection is chosen, input genes:',
                        inputrefgenesline)

        howtofilter = QComboBox()
        howtofilter.addItem('Kruskal-Wallis filtering')
        howtofilter.addItem('Wilcoxon filtering')
        howtofilter.addItem('No group filtering')
        howtofilter.currentIndexChanged.connect(self.change_howtofilter)
        self.laycnorm.addRow('How to filter/flag reference genes', howtofilter)

        additionalnormalization = QComboBox()
        additionalnormalization.addItem('No')
        additionalnormalization.addItem('Standarization (0-1)')

        additionalnormalization.currentIndexChanged.connect(self.change_adnormalization)

        self.laycnorm.addRow('Perform additional normalization?', additionalnormalization)

        howtoexport = QComboBox()
        howtoexport.addItem('Normalized count matrix')
        howtoexport.addItem('log2 count matrix')
        howtoexport.addItem('log10 count matrix')
        howtoexport.currentIndexChanged.connect(self.change_exportmethod)
        self.laycnorm.addRow('Format export results', howtoexport)

        doubleforcnorm = QHBoxLayout()
        doublefortickscnorm = QHBoxLayout()
        doubleforcnorm.addLayout(doublefortickscnorm)
        runcnorm = QPushButton('Run normalization')
        runcnorm.setIcon(QIcon(str(imgs_path / "logoguanin_96x96.png")))
        runcnorm.clicked.connect(self.runcnorm)
        doubleforcnorm.addWidget(runcnorm)

        doubleformtic1cn = QFormLayout()
        ticpopupinfolanes3 = QCheckBox()
        ticpopupinfolanes3.stateChanged.connect(self.change_showbrowsercnorm)
        doubleformtic1cn.addRow('Pop out info about ref genes', ticpopupinfolanes3)
        doublefortickscnorm.addLayout(doubleformtic1cn)

        self.laycnorm.addRow(doubleforcnorm)

        bigtabcnorm = QWidget()
        bigtabcnorm.setLayout(self.laycnorm)

        self.tabs.addTab(bigtabcnorm, QIcon(str(imgs_path / "logoguanin_96x96.png")), 'Reference genes selection')

        self.layeval = QFormLayout()
        evaltitle = QLabel('- Normalization evaluation -')
        evaltitle.setAlignment(Qt.AlignmentFlag.AlignCenter)
        evaltitle.setFont(headerfont)
        self.layeval.addRow(evaltitle)

        doubleforeval = QHBoxLayout()
        doublefortickseval = QHBoxLayout()
        doubleforeval.addLayout(doublefortickseval)
        runevalbutton = QPushButton('Run evaluation')
        runevalbutton.setIcon(QIcon(str(imgs_path / "logoguanin_96x96.png")))
        runevalbutton.clicked.connect(self.runeval)

        doubleforeval.addWidget(runevalbutton)

        self.layeval.addRow(doubleforeval)

        doubleforevaltext = QHBoxLayout()
        self.labeltext1 = QLabel('Raw RLE plot')
        self.labeltext1.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.labeltext2 = QLabel('Normalized RLE plot')
        self.labeltext2.setAlignment(Qt.AlignmentFlag.AlignCenter)

        doubleforevaltext.addWidget(self.labeltext1)
        doubleforevaltext.addWidget(self.labeltext2)
        self.layeval.addRow(doubleforevaltext)

        doubleforeval2 = QHBoxLayout()

        self.labelpix1 = QLabel('No raw RLE plot generated yet')
        self.labelpix1.setAlignment(Qt.AlignmentFlag.AlignCenter)

        self.labelpix2 = QLabel('No normalized RLE plot generated yet')
        self.labelpix2.setAlignment(Qt.AlignmentFlag.AlignCenter)

        doubleforeval2.addWidget(self.labelpix1)
        doubleforeval2.addWidget(self.labelpix2)

        self.layeval.addRow(doubleforeval2)

        doubleforevalpcatext = QHBoxLayout()
        self.labelpca1 = QLabel('Raw PCA plot')
        self.labelpca1.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.labelpca2 = QLabel('Normalized PCA plot')
        self.labelpca2.setAlignment(Qt.AlignmentFlag.AlignCenter)

        doubleforevalpcatext.addWidget(self.labelpca1)
        doubleforevalpcatext.addWidget(self.labelpca2)
        self.layeval.addRow(doubleforevalpcatext)

        doubleforevalpca2 = QHBoxLayout()

        self.labelpixpca1 = QLabel('No raw PCA plot generated yet')
        self.labelpixpca1.setAlignment(Qt.AlignmentFlag.AlignCenter)

        self.labelpixpca2 = QLabel('No normalized PCA plot generated yet')
        self.labelpixpca2.setAlignment(Qt.AlignmentFlag.AlignCenter)

        doubleforevalpca2.addWidget(self.labelpixpca1)
        doubleforevalpca2.addWidget(self.labelpixpca2)

        self.layeval.addRow(doubleforevalpca2)

        bigtabeval = QWidget()
        bigtabeval.setLayout(self.layeval)

        self.tabs.addTab(bigtabeval, QIcon(str(imgs_path / "logoguanin_96x96.png")), 'Evaluation')
        layout.addWidget(self.tabs)

        selector_lay = QHBoxLayout()

        # tabswidget = QTabWidget()
        # layout2_w = QWidget()
        # layout2_w.setLayout(layout2)
        # layout3_w = QWidget()
        # layout3_w.setLayout(layout3)
        #
        # tabswidget.addTab(layout2_w, 'Loading data')
        # tabswidget.addTab(layout3_w, 'Normalization')
        #
        # selector_lay.addWidget(tabswidget)

        # mainlayout.addLayout(selector_lay)
        # mainlayout.addLayout(layout1)
        self.setLayout(layout)


    def openselectfolder(self):
        folder = QFileDialog.getExistingDirectory(self)
        self.state.folder = pathlib.Path(folder)
        self.showfoldertextbox.setText(folder)
        logging.debug(f"Folder {folder} selected succesfully")

    def openselectoutputfolder(self):
        folder = QFileDialog.getExistingDirectory(self)
        self.state.outputfolder = pathlib.Path(folder)
        self.showoutputfoldertextbox.setText(folder)

    def opencsvfile(self):
        file = QFileDialog.getOpenFileName(self)
        self.state.groupsfile = file[0]
        if (file[0] == 'No groups defined') | (file[0] == ''):
            self.state.groups = 'no'
            self.showfiletextbox.setText('No groups defined')
        else:
            self.showfiletextbox.setText(str(self.state.groupsfile))
            self.state.groups = 'yes'
        logging.debug(f"state.groupsfile = {self.state.groupsfile}")

    def changingmodeid(self, value):
        if value == 1:
            self.state.modeid = 'sampleID'
        elif value == 0:
            self.state.modeid = 'filename'
        logging.debug(f"state.modeid = {self.state.modeid}")

    def change_showbrowserrawqc(self, checkbox):
        ##DEBUG
        #DUPLICATED TABS WHEN SHOWING INFO IN BROWSER
        self.state.showbrowserrawqc = (checkbox == 2)
        logging.debug(f"state.showbrowserrawqc = {self.state.showbrowserrawqc}")

    def change_showbrowserqc(self, checkbox):
        self.state.showbrowserqc = (checkbox == 2)
        logging.debug(f"state.showbrowserqc = {self.state.showbrowserqc}")

    def change_showbrowsercnorm(self, checkbox):
        self.state.showbrowsercnorm = (checkbox == 2)
        logging.debug(f"state.showbrowsercnorm = {self.state.showbrowsercnorm}")

    def next_tab(self):
        cur_index = self.tabs.currentIndex()
        if cur_index < len(self.tabs) - 1:
            self.tabs.setCurrentIndex(cur_index + 1)

    def runloadingrccs(self):
        self.parent.statusBar().showMessage('Loading RCC files...')
        self.parent.statusBar().repaint()
        logging.debug(f"state.groupsfile = {self.state.groupsfile}")
        logging.debug(f"state.folder = {self.state.folder}")
        self.next_tab()
        guanin.runQCview(self.state)
        self.parent.statusBar().showMessage(self.state.current_state)



    def changingbackgroundbutton(self, checkbox):
        if checkbox == 0:
            self.state.background = 'Background'
        elif checkbox == 1:
            self.state.background = 'Background2'
        elif checkbox == 2:
            self.state.background = 'Background3'
        elif checkbox == 3:
            self.state.background = 'Backgroundalt'
        elif checkbox == 4:
            self.state.background = 'Manual'
        print(f'background: {self.state.background}')
        logging.debug(f"state.background = {self.state.background}")

    def changingmanualbackground(self, value):
        self.state.manualbackground = value
        print(f'manualbackground: {self.state.manualbackground}')
        logging.debug(f"state.manualbackground = {self.state.manualbackground}")

    def changingbackgroundcorrection(self, checkbox):
        if checkbox == 0:
            self.state.lowcounts = 'subtract'
        elif checkbox == 1:
            self.state.lowcounts = 'asim'
        elif checkbox == 2:
            self.state.lowcounts = 'skip'
        print(f'background method: {self.state.lowcounts}')
        logging.debug(f"state.lowcounts = {self.state.lowcounts}")

    def changingpbelowbackground(self, value):
        self.state.pbelowbackground = value
        logging.info(f"state.pbelowbackground = {self.state.pbelowbackground}")

    def changeminfov(self, value):
        self.state.change_float("minfov", value)

    def changemaxfov(self, value):
        self.state.change_float("maxfov", value)

    def changeminbd(self, value):
        ##DEBUG
        #DUPLICATED SELF.STATE.MINBD
        self.state.change_float("minbd", value)
        self.state.minbd = float(value.replace(',', '.'))
        print(self.state.minbd)

    def changemaxbd(self, value):
        self.state.change_float("maxbd", value)
        self.state.maxbd = float(value.replace(',', '.'))
        print(self.state.maxbd)

    def changeminlin(self, value):
        self.state.change_float("minlin", value)

    def changemaxlin(self, value):
        self.state.change_float("maxlin", value)

    def changeminscaf(self, value):
        self.state.change_float("minscaf", value)

    def changemaxscaf(self, value):
        self.state.change_float("maxscaf", value)

    def changesampleremoving(self, checkbox):
        if checkbox == 0:
            self.state.laneremover = True #laneremover removes lanes
            self.state.manual_remove = False
            logging.info(f"lanes removed: {self.state.laneremover}")
        elif checkbox == 1:
            self.state.laneremover = False #no lanes are removed
            self.state.manual_remove = False
            logging.info(f"Lanes removed: {self.state.laneremover}. Check flagged samples in reports")
        elif checkbox == 2:
            self.state.laneremover = True
            self.state.manual_remove = True
            logging.info(f"Lanes to manually remove: {self.state.remove}")

    def changemanualremoveinput(self, value):
        self.state.remove = value
        logging.debug(f"state.remove = {self.state.remove}")

    def runqc(self):
        self.parent.statusBar().showMessage('Aplying filters and QC...')
        self.parent.statusBar().repaint()
        guanin.runQCfilter(self.state)
        self.next_tab()
        self.parent.statusBar().showMessage(
            'QC done, ready to perform technical normalization')

    def changekvalue(self, value):
        self.state.kvalue = value
        logging.info(f'kvalue: {self.state.kvalue}')

    def changetnormmethod(self, checkbox):
        self.laymethod.setRowVisible(2, False)
        self.doubleformtic1tn.setRowVisible(0,False)
        if checkbox == 0:
            self.state.pipeline = 'ruvgnorm'
            self.laymethod.insertRow(2, 'k value', self.inputkvalue)
            self.laymethod.setRowVisible(2, True)
            self.doubleformtic1tn.insertRow(0, self.list_text_tick[checkbox], self.ticforaftertransformlowcounts)
            self.doubleformtic1tn.setRowVisible(0, True)
        elif checkbox == 1:
            self.state.pipeline = 'scalingfactors'
            self.laymethod.insertRow(2, 'Scaling factor', self.method2combobox)
            self.laymethod.setRowVisible(2, True)
            self.doubleformtic1tn.insertRow(0, self.list_text_tick[checkbox], self.ticforaftertransformlowcounts)
            self.doubleformtic1tn.setRowVisible(0, True)
        logging.debug(f"state.tecnormeth = {self.state.tecnormeth}")

    def runthetechnorm(self):
        self.parent.statusBar().showMessage(
            'Performing technical normalization')
        self.parent.statusBar().repaint()
        guanin.pipeline1(self.state)
        self.next_tab()
        self.parent.statusBar().showMessage(
            "First step of normalization done, " +
            "ready to perform content normalization")

    def change_filhousekeepingmincounts(self, value):
        self.state.mincounthkes = int(value.replace(',', '.'))
        print(self.state.mincounthkes)
        logging.debug(f"state.mincounthkes = {self.state.mincounthkes}")

    def change_includeerg(self, checkbox):
        self.state.refendgenes = "hkes" if (checkbox == 0) else "endhkes"
        print(self.state.refendgenes)
        logging.debug(f"state.refendgenes = {self.state.refendgenes}")

    def changeaftertransformlowcounts(self, checkbox):

        self.state.tnormbeforebackgcorr = (checkbox == 2)

        logging.debug(
            f"state.tnormbeforebackgcorr = {self.state.tnormbeforebackgcorr}")

    def change_howmanyergs(self, value):
        self.state.numend = int(value.replace(',', '.'))
        logging.debug(f"state.numend = {self.state.numend}")

    def change_contnormmethod(self, checkbox):
        if checkbox == 0:
            self.state.contnorm = 'refgenes'
        elif checkbox == 1:
            self.state.contnorm = 'refgenes'
        elif checkbox == 2:
            self.state.contnorm = 'topn'
        elif checkbox == 3:
            self.state.contnorm = 'all'
        elif checkbox == 4:
            self.state.contnorm = 'ponderaterefgenes'
        elif checkbox == 5:
            self.state.contnorm = 'refgenes'
        logging.debug(f"state.contnorm = {self.state.contnorm}")

    def change_nrefgenes(self, value):
        if self.state.contnorm == 'refgenes':
            self.state.nrefgenes = int(value)
            logging.debug(f"state.nrefgenes = {self.state.nrefgenes}")
        elif self.state.contnorm == 'ponderaterefgenes':
            self.state.nrefgenes = int(value)
            logging.debug(f"status.nrefgenes = {self.state.nrefgenes}")
        elif self.state.contnorm == 'topn':
            self.state.topngenestocontnorm = int(value)
            logging.debug(
                f"state.topngenestocontnorm = {self.state.topngenestocontnorm}")

    def change_inputnamesrefgenes(self, value):
        self.state.chooserefgenes = value
        logging.debug(f"state.chooserefgenes = {self.state.chooserefgenes}")

    def change_howtofilter(self, checkbox):
        if checkbox == 0:
            self.state.filtergroupvariation = 'filterkrus'
        elif checkbox == 1:
            self.state.filtergroupvariation = 'filterwilcox'
        elif checkbox == 2:
            self.state.filtergroupvariation = 'nofilter'
        logging.debug(
            f"state.filtergroupvariation = {self.state.filtergroupvariation}")

    def change_adnormalization(self, checkbox):
        if checkbox == 0:
            self.state.adnormalization = 'no'
        elif checkbox == 1:
            self.state.adnormalization = 'standarization'
        logging.debug(f"state.adnormalization = {self.state.adnormalization}")

    def runcnorm(self):

        self.parent.statusBar().showMessage(
            'Performing content normalization...')

        self.parent.statusBar().repaint()
        logging.debug(f"state.groupsfile = {self.state.groupsfile}")
        guanin.pipeline2(self.state)
        self.next_tab()
        self.parent.statusBar().showMessage(
            "Content normalization done, ready to evaluate normalization. ")

    def runeval(self):
        self.parent.statusBar().showMessage(
            'Performing evaluation, plotting RLE...')
        self.parent.statusBar().repaint()
        self.next_tab()
        (rawiqr, normiqr) = guanin.evalnorm(self.state)
        self.labeltext1.setText('Raw RLE plot')
        self.labeltext2.setText('Normalized RLE plot')
        self.parent.statusBar().showMessage(
            "Evaluation and data export ready, check 'output' folder")

        pixmap1 = QPixmap(str(self.state.outputfolder
                              / "images" / "rlerawplot2.png"))
        pixmap2 = QPixmap(str(self.state.outputfolder
                              / "images" / "rlenormplot2.png"))
        self.labelpix1.setPixmap(pixmap1)
        self.labelpix2.setPixmap(pixmap2)

        pixmappca1 = QPixmap(str(self.state.outputfolder / "images" / "pcaraw2.png"))
        pixmappca2 = QPixmap(str(self.state.outputfolder / "images" / "pcanorm2.png"))
        self.labelpixpca1.setPixmap(pixmappca1)
        self.labelpixpca2.setPixmap(pixmappca2)

    def change_exportmethod(self, checkbox):
        if checkbox == 0:
            self.state.logarizedoutput = 'no'
        elif checkbox == 1:
            self.state.logarizedoutput = '2'
        elif checkbox == 2:
            self.state.logarizedoutput = '10'
        logging.debug(f"state.logarizedoutput = {self.state.logarizedoutput}")

    def change_showlastlog(self, checkbox):
        self.state.showlastlog = (checkbox == 2)
        logging.debug(f"state.showlastlog = {self.state.showlastlog}")

    def showflaggedlanes(self):
        self.showingflaggedlanes.setText(str(self.state.nbadlanes))

    def showingrawIQR(self):
        self.labeltext1.setText(
            f"Raw RLE plot, IQR: {self.state.rawmeaniqr}")

    def showingnormIQR(self):
        self.labeltext2.setText(
            f"Normalized RLE plot, IQR: {self.state.normmeaniqr}")

    def change_ticdeseq2_mor(self, checkbox):
        self.state.deseq2_mor = (checkbox == 2)
        self.state.tnormbeforebackgcorr = (checkbox == 2)
        print(self.state.tnormbeforebackgcorr)

    def change_scaling_factor(self, index):
        self.state.tecnormeth = self.list_scaf_tecnormeth[index]
        print(self.state.tecnormeth)

class logger(logging.Handler):
    def __init__(self, parent):
        super().__init__(parent)
        self.widget = QPlainTextEdit(parent)
        self.widget.setReadOnly(True)

    def emit(self, record):
        msg = self.format(record)
        self.widget.appendPlainText(msg)


def main(args=None):
    # Set the default output dir to TMP dir + username
    app_dir = guanin.default_output()

    # The the loggers to a File in App Data and Console
    logging_level = logging.INFO
    if args and args.debug:
        logging_level = logging.DEBUG
    logging.basicConfig(
        level=logging_level,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        datefmt="%m/%d/%Y %I:%M:%S %p",
        handlers=[
            logging.FileHandler(
                pathlib.Path(app_dir / "analysis_description.log")),
            logging.StreamHandler()
        ])

    # Launch the app
    app = QApplication(sys.argv)
    pixmap = QPixmap(
        str(pathlib.Path(__file__).parent
            / "image" / "guanin_splashscreen2_HQ.resized.png"))
    splash = QSplashScreen(pixmap)
    logging.info("New GUANIN session started")
    logging.info(f"Output dir set to {app_dir}")

    splash.setWindowFlags(Qt.WindowType.WindowStaysOnTopHint |
                          Qt.WindowType.SplashScreen)
    splash.show()
    QTimer.singleShot(1000, splash.close)
    window = MainWindow(config={"output_folder": app_dir})
    window.resize(860,480)
    window.show()
    app.exec()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Guanin GUI interface.')
    parser.add_argument(
        "--debug",
        dest="debug",
        action="store_true",
        default=False,
        help="Set logging at debug level.")
    args = parser.parse_args()
    main(args)