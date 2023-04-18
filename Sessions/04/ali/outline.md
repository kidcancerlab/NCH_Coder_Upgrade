<!--
Compile :
    pandoc -f markdown notes/sample_notes.md - -filter pandoc-crossref -t latex -o notes.pdf

Notes:
    1. http://lierdakil.github.io/pandoc-crossref/

To Do : 
    1. Add \bcancelto from https://tex.stackexchange.com/a/218430/84495
    2. checkmarks https://tex.stackexchange.com/a/132785/84495
-->


<!--
    YAML section
    \usepackage{bbm}
-->
---
title: Outline ``scRGOT $:$ Speeding it up'' 
author: Ali Snedden
date: 2022-05-20
abstract: 

...
---
header-includes:
  - \hypersetup{colorlinks=true,
            urlcolor=blue,
            pdfborderstyle={/S/U/W 1}}
    \usepackage{physics}
    \usepackage{cancel}
    \usepackage{amssymb}
    \usepackage{hyperref}
    \usepackage{listings}
    \usepackage{blochsphere}
    \usepackage{array}
    \usepackage{lstautogobble}
    \usepackage{mathtools}
---
\definecolor{codegray}{gray}{0.9}
\newcommand{\code}[1]{\colorbox{codegray}{\texttt{#1}}}

\maketitle
\tableofcontents




Outline - Slides
==================================
1. Background : 
    a) What is a computer? 
        #. Draw diagram, describe parts (e.g. Disk -> RAM -> Cache -> Register / Processor)
    #) Jargon!
#. Concepts of parallelization
    a) Mostly describe 'embarrassingly parallelizable' problems
        #. For loop, summation or matrix multiplication
        #. Do maybe rope example from SC22
    #) Just b/c it is parallel doesn't make it faster
        #. Sometimes better to request fewer processors but get all work to run concurrently
    #) Shared vs. Distributed memory parallelization
        #. If you see MPI, mvapich, etc, may not be for you
        #. OpenMP, parapply, etc. 
#. How to parallelize :
    a) E.g. use flag (e.g. hisat -p)
        #. Include plot of work
        #. Maybe quiz asking which is the most likely plot. 
    #) Using slurm to loop over samples
        #. Array job
        #. Bash loop
        #. sbatch --wrap="some cmd"
    #) parapply in R
#. Using SLURM
    #) Error handling
        #. When running many jobs, want to know when one out of many encounters problem
#. An actual Single Cell RNASeq example
    
Code - 
==================================
1. Example of parapply in R
#. Example of running hisat2 -p 
#. Example of aligning multiple files at once
#. Example of Single cell sequencing
    a) Data from ?
        #. [Human Cell Atlas](www.humancellatlas.org)
        #. [Tabula Muris](www.tabula-muris.ds.czbiohub.org)


Questions :
==================================
1. Q : What version of R will they be using
#. Q : Should we throttle how many jobs users in training group can use? If so how much?
#. Q : Which 


TO DO :
==================================
1. data set : https://github.com/kidcancerlab/CellTypeAnnRefs
#. Upload git : https://github.com/kidcancerlab/2023CoderUpgrade
