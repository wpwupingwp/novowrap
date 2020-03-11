# Quick start
   ```shell
   pip install novowrap --user
   # Windows
   python -m novowrap input_file_1 input_file_2 -taxon taxon name
   # Linux and MacOS
   python -m novowrap input_file_1 input_file_2 -taxon taxon name
   ```
   Or download [the package](https://github.com/wpwupingwp/novowrap/releases),
   unzip, and then double click `novowrap.exe` or `novowrap`.
# Table of Contents
   * [Quick start](#quickstart)
   * [Feature](#feature)
      * [Function](#function)
      * [Application](#application)
   * [Prerequisite](#prerequisite)
      * [Hardware](#hardware)
      * [Software](#software)
   * [Installation](#installation)
      * [Portable](#portable)
      * [Install with pip](#Installwithpip)
   * [Usage](#usage)
      * [Command line](#commandline)
      * [Graphical user interface](#graphicaluserintetface)
   * [Options](#options)
      * [Assembly](#assembly)
      * [Validate](#validate)
      * [Merge](#merge)
   * [Performance](#performance)
# Feature
   * Assembly chloroplast genomes from given NGS data, with minimal parameters
     to set. Also support batch mode.

     Automatic generate uniform conformation with reference (typically, start
     from _trnH-psbA_, and, SSC/LSC region have same direction with reference).
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
     * Only merge
     ```
     # Windows
     python -m novowrap.merge -input [input file]
     # Linux and MacOS
     python3 -m novowrap.merge -input [input file]
     ```
## Graphical user interface

     If installed with pip,
     ```
     # Windows
     python -m novowrap
     # Linux and MacOS
     python3 -m novowrap
     ```

     If use portable version, just double click the "novowrap.exe" or
     "novowrap" in the folder.

     Then click the button to choose which module to use. Notice that if one
     of the option was set to the wrong value, the program will refuse to run
     and hint the user to correct the invalid option.
# Options
## Assembly
### General

     These options are for general usage.

     `-h` or `-help`: print help message

     `-input [filenames]`: input filenames, could be single or pair-end,
     support gz and fastq format

     `-list [filenames]`: input list for batch mode. The list should be a csv
     file with three columns,
     ```
     Input 1,Input 2,Taxonomy
     ```
     If only have one input file, just leave the `Input 2` column empty.

     `-ref [filename]`, reference file for assembly and validate, shoud be
     `genbank` format contains *only one* chloroplast genome sequence. Extra
     sequences will be *ignored*. For automatic running, `-taxon` is
     recommended

     `-taxon [taxonomy name]`: taxonomy name of the sample, space is allowed.
     For instance, `-taxon Oryza sativa`, will find reference chloroplast
     genome of _Oryza sativa` from NCBI Refseq database. If not found, will
     find most related species' reference, `Oryza`, `Poaceae`, `Poales` et al.

     If `-ref` and `-taxon` were both not set, will use `_Nicotiana tabacum_`
     to get the reference (which is one of the earliest sequenced chloroplast
     genome

     `-out [folder name]`: output folder, if not set, the program will auto
     create it according to input file's name

### Advanced

     These options are for advanced usage. If not sure, just keep the default
     value.

     `-platform [illumina/ion]`: sequencing platform, the default is
     `illumina`. If use ion-torrent, set `-platform ion`

     `-insert_size [number]`: the insert size of sequencing library, should be
     interger

     `-seed [names]`: gene names as seeds for assembly, separated by
     comma, the default seeds are `rbcL,psaB,psaC,rrn23`

     `-seed_file [filename]`: seed file, will overwrite `-seed` option

     `-split [number]`: split input file, only use `[number]` of them, useful
     for large data while computer memory is limited.  For instance,
     `-split 10000000` will only use 10 million reads

     `-kmer [number]`: kmer size for assembly, should be odd number. Most of
     time it's unnecessary to change it

     `-min [number]`: minimum genome size, default is 100 kB

     `-max [number]`: maximum genome size, default is 200 kB. Only change
     `-min` and `-max` if target genome size is out of the default range. The
     program needn't to know the precise size of the genome

     `-mem [number]`: memory limit, the unit is GB. For instance, `-mem 8`
     will limit the memory usage to 8 GB. Should be integer

     `-perc_identity [number]`: the threshold of minimum percent of identity,
     used for validation with BLAST. The default value is `0.7`. Should be
     float number between 0 and 1.

     `-len_diff [number]`: the threshold of maximum percent of length
     different of query and reference. Used for exliminate invalid assembly
     result. If the sequence length's differnce of assembly and reference
     genome is larger than the value, the assembly result will be discard. The
     default value is `0.2`. Should be float number between 0 and 1.

     `-debug`: print debug information if set

## Validate
### General

     `-h` or `-help`: print help message

     `-input [filename]`: input filename. Only support `fasta` format

     `-ref [filename]`, reference file for assembly and validate, shoud be
     `genbank` or `fasta` format that contains *only one* chloroplast genome
     sequence. Extra sequences will be *ignored*.

     `-taxon [taxonomy name]`: taxonomy name of the the reference's species,
     space is allowed. Recommond to use same genus or family, or higher rank
     if it's well known that the target taxonomy's chloroplast genome is
         conserved.

     `-out [folder name]`: output folder, if not set, the program will auto
     create it according to input file's name
### Advanced

     `-perc_identity [number]`: the threshold of minimum percent of identity,
     used for validation with BLAST. The default value is `0.7`. Should be
     float number between 0 and 1.

     `-len_diff [number]`: the threshold of maximum percent of length
     different of query and reference. Used for exliminate invalid assembly
     result. If the sequence length's differnce of assembly and reference
     genome is larger than the value, the assembly result will be discard. The
     default value is `0.2`. Should be float number between 0 and 1.

     `-debug`: print debug information if set

## Merge

     `-h` or `-help`: print help message

     `-input [filename]`: input filename. Only support `fasta` format

     `-out [folder name]`: output folder, if not set, the program will auto
     create it according to input file's name

# Performance
    The most time-consuming step is assembly. If the chloroplast genome's
    reads in sequencing data is plentiful enough and the computer's memory is
    big enough for the data size, the assembly will be finished in minutes.

    The validation step usually could finished in less than one minute. If
    slower, please check the Internet connection since the program will query
    the NCBI database.

    The merge module could cost seconds or minutes. It depends on input data.
    Complex relationship of contigs requires much more time.

# Reference
# License
# Q&A




