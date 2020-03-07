# Quick start
   ```
   pip install novowrap --user
   python -m novowrap input_file_1 input_file_2 -taxon taxon name
   ```
   Or download [this](link), unzip, and then double click novowrap.exe or
   novowrap.
# Table of Contents
   * [Quick start](#quick start)
   * [Feature](#feature)
      * [Function](#function)
      * [Application](#application)
   * [Prerequisite](#prerequisite)
      * [Hardware](#hardware)
      * [Software](#software)
   * [Installation](#installation)
      * [Portable](#portable)
      * [Install with pip](#Install with pip)
   * [Usage](#usage)
      * [Quick examples](#quick-examples)
      * [Input](#input)
      * [Sequence ID](#sequence-id)
      * [Output](#output)
   * [Options](#options)
      * [Help](#help)
      * [General](#general)
      * [Genbank](#genbank)
      * [Pre-process](#preprocess)
      * [Evaluate](#evaluate)
      * [Primer Design](#primer-design)
   * [Performance](#performance)
# Feature
   * Assembly chloroplast genomes from given NGS data, with minimal parameters
     to set. Also support batch mode.

     Automatic generate uniform conformation with reference (typically, start
     from _trnH-psbA_ and SSC/LSC region have same direction with reference).
   * Merge contigs according to overlapping region. May handle Invert-Repeat
     fragments.
   * Validate assembly results by comparing the synteny and sequence homology
     with given reference (or taxonomy name).
# Prerequisite
## Hardware
     The assembly function will calls NOVOPlasty, which requires 2 GB memory
     for 1 GB uncompressed data.

     The other functions could run in normal computer and have no extra
     requirements for memory, CPU, et al.

     The software requires Internet for the first run to install the lacking
     dependencies (see #software). It could works if offline, but better when
     have connection.
## Software
     For portable version, no dependencies.

     For installing from pip, [Python](https://www.python.org/downloads/) is
     required. Notice that the python version should equal to or newer than
     **3.6**.
# Installation
## Portable
     Download from the [link](url), unpack and run with Internet for the first
     time.
## Install with pip
     1. Install Python from the [link](https://www.python.org/downloads/). 3.6
        or newer is required.
     
     2. Open command line, run
     ```
     pip install novowrap --user
     ```
# Usage
## Command line
     In Linux and MacOS, Python 2 is "python2" and Python 3 is "python3".
     However, in Windows, Python 3 is called "python", too. Please notice the
     difference.

     * Show help information of each module
     ```
     # Windows
     python -m novowrap.assembly -h
     python -m novowrap.validate -h
     python -m novowrap.merge -h
     # Linux and MacOS
     python3 -m novowrap.assembly -h
     python3 -m novowrap.validate -h
     python3 -m novowrap.merge -h
     ```
     * Assembly and validate
     ```
     # Windows, single sample
     python -m novowrap -i [input1] [input2] -taxon [taxonomy]
     # Windows, batch mode for numerous samples
     python -m novowrap -l [list file]
     # Linux and MacOS, single sample
     python3 -m novowrap -i [input1] [input2] -taxon [taxonomy]
     # Linux and MacOS, batch mode for numerous samples
     python3 -m novowrap -l [list file]
     ```
     * Only validate
     ```
     # Windows
     python -m novowrap.validate -input [input file] -taxon [taxonomy]
     # Linux and MacOS
     python3 -m novowrap.validate -input [input file] -taxon [taxonomy]
     ```




## GUI






