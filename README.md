# Quick start
   ```
   pip install novowrap --user
   python -m novowrap
   ```
   Or download [this](link), unzip, and then double click novowrap.exe or
   novowrap.
# Table of Contents
   * [Quick start](#quickstart)
   * [Feature](#feature)
      * [Function](#function)
      * [Application](#application)
   * [Prerequisite](#prerequisite)
      * [Hardware](#hardware)
      * [Software](#software)
   * [Installation](#installation)
      * [Software](#software)
      * [Python module](#python-module)
      * [Internet](#internet)
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







