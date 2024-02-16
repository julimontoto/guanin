import argparse
import pathlib
import time
import tempfile


def argParserCli():
    parser = argparse.ArgumentParser(
        description="Nanostring quality control analysis"
    )
    parser.add_argument(
        "-f",
        "--folder",
        type=str,
        default=pathlib.Path.cwd() / "../examples/d1_COV_GSE183071",
        help="relative folder where RCC set is located. Default: /data",
    )
    parser.add_argument(
        "-minf",
        "--minfov",
        type=float,
        default=0.75,
        help="set manually min fov for QC",
    )
    parser.add_argument(
        "-maxf",
        "--maxfov",
        type=float,
        default=1,
        help="set manually max fov for QC",
    )
    parser.add_argument(
        "-minbd",
        "--minbd",
        type=float,
        default=0.1,
        help="set manually min binding density for QC",
    )
    parser.add_argument(
        "-maxbd",
        "--maxbd",
        type=float,
        default=1.8,
        help="set manually max binding density for QC",
    )
    parser.add_argument(
        "-minlin",
        "--minlin",
        type=float,
        default=0.75,
        help="set manually min linearity for QC",
    )
    parser.add_argument(
        "-maxlin",
        "--maxlin",
        type=float,
        default=1,
        help="set manually max linearity for QC",
    )
    parser.add_argument(
        "-minscaf",
        "--minscalingfactor",
        type=float,
        default=0.3,
        help="set manually min scaling factor for QC",
    )
    parser.add_argument(
        "-maxscaf",
        "--maxscalingfactor",
        type=float,
        default=3,
        help="set manually max scaling factor for QC",
    )
    parser.add_argument(
        "-swbrrq",
        "--showbrowserrawqc",
        type=bool,
        default=False,
        help="pops up infolanes and qc summary",
    )
    parser.add_argument(
        "-swbrq",
        "--showbrowserqc",
        type=bool,
        default=False,
        help="pops up infolanes and qc summary",
    )
    parser.add_argument(
        "-swbrcn",
        "--showbrowsercnorm",
        type=bool,
        default=False,
        help="pops up infolanes and qc summary",
    )
    parser.add_argument(
        "-lc",
        "--lowcounts",
        type=str,
        default="skip",
        choices=["skip", "asim", "subtract"],
        help="what to do with counts below background?",
    )
    parser.add_argument(
        "-mi",
        "--modeid",
        type=str,
        default="filename",
        choices=["sampleID", "filename", "id+filename"],
        help="choose sample identifier. sampleID: optimal if assigned in "
        + "rccs. filenames: easier to be unique. id+filename: care with group "
        + "assignment coherence",
    )
    parser.add_argument(
        "-mv",
        "--modeview",
        type=str,
        default="view",
        choices=["justrun", "view"],
        help="choose if plot graphs or just run calculations",
    )
    parser.add_argument(
        "-tnm",
        "--tecnormeth",
        type=str,
        default="posgeomean",
        choices=["posgeomean", "Sum", "Median", "regression"],
        help="choose method for technical normalization",
    )
    parser.add_argument(
        "-reg",
        "--refendgenes",
        type=str,
        default="endhkes",
        choices=["hkes", "endhkes"],
        help="choose refgenes, housekeeping, or hkes and endogenous",
    )
    parser.add_argument(
        "-re",
        "--remove",
        type=str,
        nargs="+",
        default=None,
        help="lanes to be removed from the analysis",
    )
    parser.add_argument(
        "-bg",
        "--background",
        type=str,
        default="Background",
        choices=[
            "Background",
            "Background2",
            "Background3",
            "Backgroundalt",
            "Manual",
        ],
        help="choose background: b1=meancneg+(2*std), b2=maxcneg, "
        + "b3=meancneg, balt=uses alternative subset of negative controls",
    )
    parser.add_argument(
        "-pbb",
        "--pbelowbackground",
        type=int,
        default=85,
        help="if more than %bb genes are below background, "
        + "sample gets removed from analysis",
    )
    parser.add_argument(
        "-mbg",
        "--manualbackground",
        type=float,
        default=None,
        help="set manually background",
    )
    parser.add_argument(
        "-crg",
        "--chooserefgenes",
        type=list,
        nargs="+",
        default=None,
        help="list of strings like. choose manualy reference genes to use "
        + "over decided-by-program ones",
    )
    parser.add_argument(
        "-fgv",
        "--filtergroupvariation",
        type=str,
        default="filterkrus",
        choices=[
            "filterkrus",
            "filterwilcox",
            "flagkrus",
            "flagwilcox",
            "nofilter",
        ],
        help="¿filter or flag preselected ref genes by significative "
        + "group-driven differences? needs groups to be declared",
    )
    parser.add_argument(
        "-fsn",
        "--featureselectionneighbors",
        type=float,
        default=4,
        help="number of neighbors for feature selection analysis of refgenes. "
        + "recommended 3-6",
    )
    parser.add_argument(
        "-g",
        "--groups",
        type=str,
        default="yes",
        choices=["yes", "no"],
        help="defining groups for kruskal/wilcox/fs analysis?",
    )
    parser.add_argument(
        "-ne",
        "--numend",
        type=int,
        default=6,
        help="number of endogenous tofind by ERgene to include in analysis to "
        + "check viability as refgenes",
    )
    parser.add_argument(
        "-ar",
        "--autorename",
        type=bool,
        default=False,
        help="turn on when sample IDs are not unique, be careful on sample "
        + "identification detail",
    )
    parser.add_argument(
        "-cn",
        "--contnorm",
        type=str,
        default="refgenes",
        choices=["ponderaterefgenes", "refgenes", "all", "topn"],
    )
    parser.add_argument(
        "-an",
        "--adnormalization",
        type=str,
        default="no",
        choices=["no", "standarization"],
        help="perform additional normalization? standarization available",
    )
    parser.add_argument(
        "-tn",
        "--topngenestocontnorm",
        type=int,
        default=100,
        help="set n genes to compute for calculating norm factor from top n "
        + "expressed endogenous genes",
    )
    parser.add_argument(
        "-mch",
        "--mincounthkes",
        type=int,
        default=80,
        help="set n min counts to filter hkes candidate as refgenes",
    )
    parser.add_argument(
        "-nrg",
        "--nrefgenes",
        type=int,
        default=None,
        help="set n refgenes to use, overwriting geNorm calculation",
    )
    parser.add_argument(
        "-lr",
        "--laneremover",
        type=bool,
        default=True,
        choices=[True, False],
        help="option to perform analysis with all lanes if set to no",
    )
    parser.add_argument(
        "-lo",
        "--logarizedoutput",
        type=str,
        default="no",
        choices=["2", "10", "no"],
        help="want normed output to be logarized? in what logbase?",
    )
    parser.add_argument(
        "-le",
        "--logarizeforeval",
        type=str,
        default="10",
        choices=["2", "10", "no"],
        help="logarithm base for RLE calculations",
    )
    parser.add_argument(
        "-gf",
        "--groupsfile",
        type=str,
        default="../examples/groups_d1_COV_GSE183071.csv",
        help="enter file name where groups are defined",
    )
    parser.add_argument("-st", "--start_time", type=float, default=time.time())
    parser.add_argument("-cs", "--current_state", type=str, default="Ready")
    parser.add_argument(
        "-ftl",
        "--tnormbeforebackgcorr",
        type=int,
        default=1,
        help="0= False, 1= True, 2= ruvg",
    )
    parser.add_argument(
        "-of",
        "--outputfolder",
        type=str,
        default=tempfile.gettempdir() + "/guanin_output",
    )
    parser.add_argument("-sll", "--showlastlog", type=bool, default=False)
    parser.add_argument("-rgs", "--refgenessel", type=list, default=[])
    parser.add_argument("-im2", "--indexmethod2", type=int, default=0)
    parser.add_argument("-k", "--kvalue", type=int, default=3)
    parser.add_argument(
        "-pip",
        "--pipeline",
        type=str,
        default="ruvgnorm",
        choices=["ruvgnorm", "scalingfactors"],
    )
    parser.add_argument("-wrg", "--whatrefgenes", type=list, default=[])
    parser.add_argument("-m", "--eme", type=object, default=None)
    parser.add_argument("-gpca", "--grouppca", type=str, default="GROUP")
    parser.add_argument("-dm", "--deseq2_mor", type=bool, default=True)
    parser.add_argument("-mr", "--manual_remove", type=bool, default=False)
    parser.add_argument(
        "-nbl", "--nbadlanes", type=str, default="No badlanes detected"
    )
    parser.add_argument("-bl", "--badlanes", type=set, default=set())
    parser.add_argument("-e", "--elapsed", type=float, default=0.0)

    return parser.parse_args()
