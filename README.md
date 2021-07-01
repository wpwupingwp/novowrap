[![Build Status](https://travis-ci.org/wpwupingwp/novowrap.svg?branch=master)](https://travis-ci.org/wpwupingwp/novowrap)
[![PyPI version](https://badge.fury.io/py/novowrap.svg)](https://badge.fury.io/py/novowrap)
[![Anaconda version](https://anaconda.org/wpwupingwp/novowrap/badges/version.svg)](https://anaconda.org/wpwupingwp/novowrap)
# Quick start
Download [the package](https://github.com/wpwupingwp/novowrap/releases),
unzip, and then double-click `novowrap.exe` or `novowrap`.

__OR__

Make sure you have [Python](https://www.python.org/) (3.7 or higher) or 
[conda](https://docs.conda.io/en/latest/miniconda.html) installed.

Open terminal, run
   ```shell
   # Install, using pip (recommended)
   pip install novowrap --user
   # Or, use conda
   conda install -c wpwupingwp novowrap

   # Initiliaze with Internet
   # Windows
   python -m novowrap init
   # Linux and MacOS
   python3 -m novowrap init

   # Run
   # Windows
   python -m novowrap
   # Linux and MacOS
   python3 -m novowrap
   ```

# Table of Contents
   * [Quick start](#quickstart)
   * [Feature](#feature)
   * [Prerequisite](#prerequisite)
      * [Hardware](#hardware)
      * [Software](#software)
   * [Installation](#installation)
      * [Portable](#portable)
      * [Install with pip](#Installwithpip)
      * [Install with conda](#Installwithconda)
      * [Initialization](#Initialization)
   * [Usage](#usage)
      * [Command line](#commandline)
      * [Graphical user interface](#graphicaluserintetface)
   * [Input](#input)
   * [Output](#output)
   * [Options](#options)
      * [Assembly](#assembly)
      * [Validate](#validate)
      * [Merge](#merge)
   * [Performance](#performance)
# Feature
:heavy_check_mark: Assembly chloroplast genomes from given NGS data, with minimal
parameters to set. Also, it supports batch mode.  

Automatic generate uniform structure with reference (typically, start from 
_trnH-psbA_, and, SSC/LSC region have the same direction with reference).

:heavy_check_mark: Merge contigs according to overlapping region. May handle
Invert-Repeat fragments.

:heavy_check_mark: Validate assembly results by comparing the synteny and sequence
homology with given reference (or taxonomy name).
# Prerequisite
## Hardware
The assembly function will call NOVOPlasty, which requires 2 GB memory for 1
GB uncompressed data.

The other functions could run in normal computers and have no extra
requirements for memory, CPU, et al.

The software requires Internet for the first run to install the missing
dependencies. Then, it could work if offline, but better with the internet.
## Software
For the portable version, nothing needs to be installed manually.

For installing from pip, [Python](https://www.python.org/downloads/) is
required. Notice that the python version should be **3.7** or higher.

:white_check_mark: All third-party dependencies will be automatically
installed with the Internet, including `biopython`, `matplotlib`, `coloredlogs`,
`graphviz` (python packages), and
[perl(for Windows only)](http://strawberryperl.com/download/5.30.1.1/strawberry-perl-5.30.1.1-64bit-portable.zip), 
[NOVOPlasty](https://github.com/ndierckx/NOVOPlasty),
[BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download).

# Installation
## Portable
Download from the [link](https://github.com/wpwupingwp/novowrap/releases),
unpack and run with Internet for the first time.
## Install with pip
1. Install [Python](https://www.python.org/downloads/). *3.7 or newer* is
   required.
     
2. Open the command line, run
```shell
pip install novowrap --user
```
## Install with conda
After installed [conda](https://docs.conda.io/en/latest/miniconda.html), run
```
conda install -c wpwupingwp novowrap
```
## Initialization
During the first running, `novowrap` will check and initialize the running environment. 
Missing dependencies will be automatically installed.
This step requires the Internet connection.
```shell
# Windows
python -m novowrap init
# Linux and MacOS
python3 -m novowrap init
```
# Usage
## Command line
:exclamation: In Linux and macOS, Python2 is `python2` and Python3 is `python3`.  However,
in Windows, Python3 is called `python`, too. Please notice the difference.

 * Show help information of each module
 ```shell
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
 ```shell
 # Windows, single sample
 python -m novowrap -input [input1] [input2] -taxon [taxonomy]
 # Windows, batch mode for numerous samples
 python -m novowrap -list [list file]
 # Linux and MacOS, single sample
 python3 -m novowrap -input [input1] [input2] -taxon [taxonomy]
 # Linux and MacOS, batch mode for numerous samples
 python3 -m novowrap -list [list file]
 ```
 * Only validate
 ```shell
 # Windows
 python -m novowrap.validate -input [input file] -taxon [taxonomy]
 # Linux and MacOS
 python3 -m novowrap.validate -input [input file] -taxon [taxonomy]
 ```
 * Only merge
 ```shell
 # Windows
 python -m novowrap.merge -input [input file]
 # Linux and MacOS
 python3 -m novowrap.merge -input [input file]
 ```
## Graphical user interface
 If installed with pip,
```shell
# Windows
python -m novowrap
# Linux and MacOS
python3 -m novowrap
```

If use the portable version, just double-click the `novowrap.exe` or `novowrap` in
the folder.

Then click the button to choose which module to use. Notice that if one of the
option was set to the wrong value, the program will refuse to run and hint the
user to correct the invalid option.
# Input
The `assembly` module accepts `gz` or `fastq` format as input. If use `input
list`, the list file should be `csv` format. If use `reference` file instead
of automatically get from NCBI, the file format should be `genbank`.

The `merge` module accepts `fasta` format as input.

The `validate` module accepts `fasta` format as input. If use `reference` file
instead of automatically get from NCBI, the file format should be `genbank`
**or** `fasta` as long as it is a complete chloroplast genome.

# Output
`.gb` files: **genbank** format sequence, with annotation of the boundary of
LSC/SSC/IR regions.

`.rotate` files: rotated sequence as **fasta** format; start from `trnH-psbA`, same direction with
reference

`.pdf` files: figure of validation of assembly

`_RC_` files: if filenames contain `_RC_`, it means one of the regions of the
sequence was adjusted according to the reference. The unadjusted sequence
could be found in `Temp` folder.

# Options
## Assembly
### General
These options are for general usage.

`-h` or `-help`: print help message

`-input [filenames]`: input filenames, could be single or pair-end, support gz
and fastq format

`-list [filenames]`: input list for batch mode. The list should be a csv file
with three columns,
```csv
Input 1,Input 2,Taxonomy
```
If only has one input file, just leave the `Input 2` column empty.

Please use *full path* of file names, for instance, `d:\data\sample-1\forward.fastq`
instead of `forward.fastq` or `sample-1\forward.fastq`.

`-ref [filename]`: reference file for assembly and validation, should be
`genbank` format contains *only one* chloroplast genome sequence. Extra
sequences will be *ignored*. For automatic running, `-taxon` is recommended

`-taxon [taxonomy name]`: taxonomy name of the sample, space is allowed.  For
instance, `-taxon Oryza sativa`, will find reference chloroplast genome of
`_Oryza sativa_` from NCBI RefSeq database. If not found, will find most related
species' reference, `Oryza`, `Poaceae`, `Poales` et al.

If `-ref` and `-taxon` were both not set, will use `_Nicotiana tabacum_`
to get the reference (which is one of the earliest sequenced chloroplast
genomes)

`-out [folder name]`: output folder, if not set, the program will auto
create it according to the input file's name

### Advanced
These options are for advanced users. If not sure, just keep the default
value.

`-platform [illumina/ion]`: sequencing platform, the default is `illumina`. If
use ion-torrent, set `-platform ion`

`-insert_size [number]`: the insert size of the sequencing library, should be
integer

`-seed [names]`: gene names as seeds for assembly, separated by commas, the
default seeds are `rbcL,psaB,psaC,rrn23`

`-seed_file [filename]`: seed file, will overwrite `-seed` option

`-split [number]`: split input file, only use `[number]` of them, useful for
large data while computer memory is limited.  For instance, `-split 10000000`
will only use 10 million reads

`-kmer [number]`: kmer size for assembly, should be an odd number. Most time
it's unnecessary to change it

`-min [number]`: minimum genome size, default is 100 kB

`-max [number]`: maximum genome size, default is 200 kB. Only change `-min`
and `-max` if the target genome size is out of the default range. The program
needn't know the precise size of the genome

`-mem [number]`: memory limit, the unit is GB. For instance, `-mem 8` will
limit the memory usage to 8 GB. Should be integer

`-perc_identity [number]`: the threshold of minimum percent of identity, used
for validation with BLAST. The default value is `0.7`. Should be a float number
between 0 and 1.

`-len_diff [number]`: the threshold of maximum percent of length different of
query and reference. Used for eliminating invalid assembly results. If the
sequence length's difference of assembly and reference genome is larger than
the value, the assembly result will be discarded. The default value is `0.2`.
Should be a float number between 0 and 1.

`-debug`: print debug information if set

`-mt`: for mitochondria genomes (experimental function)

`-simple_validate`: for chloroplast genomes without quadripartite structure

## Validate
### General
`-h` or `-help`: print help message

`-input [filename]`: input filename. Only support `fasta` format

`-ref [filename]`, reference file for assembly and validation, should be
`genbank` or `fasta` format that contains *only one* chloroplast genome
sequence. Extra sequences will be *ignored*.

`-taxon [taxonomy name]`: taxonomy name of the reference's species, space
is allowed. Recommend to use same genus or family, or higher rank if it's well
known that the target taxonomy's chloroplast genome is conserved.

`-out [folder name]`: output folder, if not set, the program will automatically 
create it according to input file's name
### Advanced
`-perc_identity [number]`: the threshold of minimum percent of identity, used
for validation with BLAST. The default value is `0.7`. Should be a float number
between 0 and 1.

`-len_diff [number]`: the threshold of maximum percent of length different of
query and reference. Used for eliminating invalid assembly results. If the
sequence length's difference of assembly and reference genome is larger than
the value, the assembly result will be discarded. The default value is `0.2`.
Should be a float number between 0 and 1.

`-debug`: print debug information if set

`-mt`: for mitochondria genomes (experimental function)

`-simple_validate`: for chloroplast genomes without quadripartite structure

## Merge
`-h` or `-help`: print help message

`-input [filename]`: input filename. Only support `fasta` format

`-out [folder name]`: output folder, if not set, the program will auto create
it according to the input file's name

# Performance
The most time-consuming step is assembly. If the chloroplast genome's reads in
sequencing data is plentiful enough, and the computer's memory is big enough
for the data size, the assembly will be finished in minutes.

The validation step usually could finish in less than one minute. If slower,
please check the Internet connection since the program may query the NCBI
database.

The merge module could cost seconds or minutes. It depends on input data.
The complex relationship of contigs requires much more time.

# Citation
Wu, P., Xu, C., Chen, H., Yang, J., Zhang, X. and Zhou, S. (2021), NOVOWrap: An automated solution for plastid genome assembly and structure standardization. Molecular Ecology Resources. https://doi.org/10.1111/1755-0998.13410
# License
The software itself is licensed under
[AGPL-3.0](https://github.com/wpwupingwp/novowrap/blob/master/LICENSE) (**not include third-party
software**).

# Q&A
Please submit your questions in the
[Issue](https://github.com/wpwupingwp/novowrap/issues) page :smiley:

* Q: I can't see the full UI, some part was missing.

  A: Please try to drag the corner of the window to enlarge it. We got reports
  that some users in macOS have this issue.

* Q: I got the error message that the program failed to install
  perl/BLAST/NOVOPlasty.

  A: Uncommonly, users in a specific area will have connection issue for those
  websites. Users have to manually download packages and install (see
  [Software](#software) for the download links).

  For Windows users, please download and unpack files into
  `%HOMEDRIVE%%HOMEPATH%/.novowrap`.

  For Linux  and MacOS users, please download and unpack files into
  `~/.novowrap`.

* Q: I got the error message that I don't have `tkinter `module installed.

  A: If you want to run GUI on Linux or macOS, this error may happen because the
  Python you used did not include tkinter as default package (kind of weird). Run
  ```
  # Debian and Ubuntu
  sudo apt install python3-tk
  # CentOS
  sudo yum install python3-tk
  ```
  may help.
  
  For macOS users or linux users without root privilege, please try to install the 
  newest version of Python or to use conda, see [conda](https://docs.conda.io/en/latest/miniconda.html) 
  and [Python](https://www.python.org/download/mac/tcltk/) for details.
  
* Q: It says my input is invalid, but I'm sure it's OK!

  A: Please check your files' path. The `space` character in the folder name
  or filename may cause this error.
* Q: It says `ImportError: Bio.Alphabet has been removed from Biopython` and 
  the program failed to start.
  
  A: In *2020.9*, Biopython removed Bio.Alphabet module in 
  [v1.78](https://biopython.org/wiki/Alphabet), which may cause this trouble 
  in the old version of novowrap. Please upgrade your `novowrap` to `v0.97` 
  or higher. If you find difficult to upgrade novowrap, please try to use the 
  portable packages.
* Q: I want to assemble mitochondria genomes.

  A: add `-mt` option or click the related checkbutton on the GUI.
  Since mitochondria genomes do not have a stable 
  and uniform structure like chloroplast, wet-lab experiments may be necessary 
  for verification.
* Q: I want to assemble chloroplast genomes without quadripartite structure.

  A: add `-simple_validate` option in command-line or click the related checkbutton 
  on the GUI. Note that without quadripartite structure, the Validate module 
  will skip the adjustment of the structure of the sequences. We have test this function on
  _Cytinus hypocistis_, a parasitic plant having a 19kb plastid genome.
* Q: I am a `conda` user...

  A: Install `novowrap` with `conda install -c wpwupingwp novowrap` and the usage is the same.
  In order to avoid potential conflicts with other packages, it is highly recommended 
  to create a new running environment with `conda` before installation. For example,
  ```
  conda create -n test
  conda activate test
  conda install -c wpwupingwp novowrap
  ```
  
  If fail, try to download conda packages from [this link](https://github.com/wpwupingwp/novowrap/releases),
  and run
  ```shell
  conda install [filename]

  ```
* Q: I want to use novowrap in parallel.

  A: Ideally, each novowrap instance is independent. However, the initialization process 
  could affect other running instances. To avoid this, please running "python3 -m novowrap init"
  one time before running in parallel.
* Q: The assembly failed.

  A: Make sure your sequencing data is okay. Otherwise, try to set the "-min" and "-max" options to 
  restrict the size of the target genome. Sometimes it may fail due to inappropriate size values. The 
  default size is 100-200 kb, which is okay for the common situation.
