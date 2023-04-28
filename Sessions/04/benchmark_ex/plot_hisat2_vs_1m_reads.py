# Author : Ali Snedden
# Date   : 4/20/23
# License: MIT
#
# Questions
#
# Thoughts :
#
# Note :
"""Generate plot for Ali's talk
"""
import argparse
import pandas as pd
from matplotlib import pyplot as plt

# E.g. python ../benchmark_ex/plot_hisat2_vs_1m_reads.py --path ../benchmark_ex/hisat2_vs_1m_reads.txt --outpath hisat2_align_1m_reads.pdf
def main():
    """

    Args:

    Returns:

    Raises:
    """
    parser = argparse.ArgumentParser(description="Generate plot of aligning synthetic data using hisat2 as a function processor number")
    parser.add_argument('--path', metavar='path', type=str, nargs='*',
                        help='Path to input file')
    parser.add_argument('--outpath', metavar='outpath', type=str, nargs='*',
                        help='Path to output pdf/file (include extension)')
    args = parser.parse_args()
    path = args.path[0]
    outpath = args.outpath[0]
    df = pd.read_csv(path, sep=r"\s+", header=0)
    fig = plt.figure()
    gs = fig.add_gridspec(2,1)
    ax = fig.add_subplot(gs[0,0])
    #
    #ax.set_yticks(ticks=[i for i in range(subM.shape[1])],
    #              labels=electSortL[0:subM.shape[1]]) # Less electrodes for debugging
    #ax.set_title("{}".format(Stage))
    #ax = fig.add_subplot(gs[1,0])
    ax.plot(df['cpu'], df['time'], marker='o', linestyle='--')
    ax.set_xlabel("number of processors")
    ax.set_ylabel("Time (s)")
    ax.set_title("Aligning 1 million single end reads with Hisat2")
    fig.savefig(outpath, dpi=200)
    plt.show()

if __name__ == "__main__":
    main()
