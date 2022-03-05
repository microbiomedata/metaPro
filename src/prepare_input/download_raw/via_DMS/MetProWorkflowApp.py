from Workflow import Workflow
from utility.utils import str2bool

from argparse import RawTextHelpFormatter, ArgumentParser

if __name__ == "__main__":

    parser = ArgumentParser(formatter_class=RawTextHelpFormatter)
    parser.add_argument(
        "-M",
        "--Mode",
        help="\nDifferent Modes to run the workflow?\n"
        "\t  Developer : Automatically generates files at each step!\n"
        "\t  User      : Generates CrossTab/Metric only!\n",
        type=str,
    )
    parser.add_argument(
        "-It",
        "--InputType",
        help="\n"
        "1 : a datapackage ID\n"
        "2 : a list of dataset IDs\n"
        "3 : a list of MSGFjobs Nums\n",
        type=int,
        choices=[1, 2, 3],
    )
    parser.add_argument(
        "-S",
        "--Storage",
        help="Give path where you want to store data & results of pipeline ",
        type=str,
    )
    parser.add_argument(
        "-P", "--ProjectName", nargs="?", help="FICUS project name", type=str
    )
    parser.add_argument(
        "-I",
        "--Input",
        help="\nA valid input\n"
        "InputType: 1, An Integer\n"
        "InputType: 2, A comma-seperated list of Integers\n"
        "InputType: 3, A comma-seperated list of Integers\n",
        type=str,
    )
    parser.add_argument(
        "-C",
        "--CombineDatasets",
        type=str2bool,
        nargs="?",
        const=True,
        default=False,
        help="Combine all dataset's MSGF & MASIC jobs to single file for generating crossTabs.",
    )
    parser.add_argument(
        "-Sa",
        "--SelectAnalysis",
        help="\n"
        "internal : Run internal analysis\n"
        "ficus    : Run ficus analysis\n"
        "both     : Run both analysis\n",
        type=str,
        choices=["internal", "ficus", "both"],
    )
    parser.add_argument(
        "-PLM",
        "--PipeLineMode",
        help="For each dataset: \n"
        "PNNL: Download data[MASIC{SYNOPSIS}, MSGF+{SYN}] from DMS\n"
        "NMDC=: Download data[.raw{DMS}, parameter{DMS}, fasta{JGI}] \n"
        "And run merge\n",
        type=str,
        choices=["PNNL", "NMDC"],
    )
    args = parser.parse_args()
    # start= Workflow.__init__(args.Mode, args.InputType, args.Storage ,args.ProjectName,
    #                  args.Input, args.CombineDatasets, args.SelectAnalysis, args.PipeLineMode)
    start = Workflow(
        args.Mode,
        args.InputType,
        args.Storage,
        args.ProjectName,
        args.Input,
        args.CombineDatasets,
        args.SelectAnalysis,
        args.PipeLineMode,
    )
    print(args.Input)
    start.start_workflow()
