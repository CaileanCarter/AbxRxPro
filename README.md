# <b>AbxRxPro:</b> antibiotic resistance profiler

AbxRxPro is a web-page visualisation tool for phenotypic antibiotic resistance using the graphing module <em>Plotly</em>. The plot can be annotated with genotypic data from antibiotic resistance genes identifying programmes like RGI (resistance gene identifier), staramr and amrfinder. The AbxRxPro package comes with a Galaxy workflow for the inputting of genome sequences through these antibiotic resistance genes identifiers. Frequencies for genes across your given isolates are plotted in a seperate tab.

<img src="https://github.com/CaileanCarter/AbxRxPro/pic/main_plot.jpg"
     alt="Main display"
     style="float: left; margin-right: 10px;"/>

<br>

## Features
- Save plots as profiles for faster loading
- Can specify which antibiotics to plot as a commandline argument
- Show genes present in your isolate which confer resistance to the antibiotic class
- Accepts genotypic data sources: staramr, RGI and amrfinder  
- Plots gene frequencies and in which isolates
- Can personalise the colour scheme and have text overlay
- Can export your plots as interactive HTML files which load in your web browser (via Plotly)

<br>

## Requirements
- Python 3.8.2 (may also be compatible with earlier and later versions of Python 3)
- Plotly 4.8.1
- Pandas 1.1.3
- Chrome web browser

<br>

## Input files
### Phenotype file
The main input file to create any antibiotic resistance profile is a 'phenotype file' which contains the information regarding phenotypic resistance. This file should be in an Excel format with antibiotic names as column headers and isolate IDs as row headers. Specify phenotypic antibiotic resistance status using the following identifiers:

Notation | Status
:---: | :---:
R | Resistant
I | Intermediate
S | Susceptible
U | Undetermined <i>optional</i> 

Missing cells will be assigned 'U' (undetermined). 

<img src="https://github.com/CaileanCarter/AbxRxPro/pic/Pheno_template.jpg"
     alt="Phenotypic Excel file"
     style="float: left; margin-right: 10px;"/>



<br> <br>
To create a new profile, you will need to specify the phenotypic file using the `-P` flag:
```
abxrxpro -P path\to\pheno\phenotypic_data.xlsx
```

If your file path contains spaces you will need to put the file path in double quotation marks `"path\to\pheno\phenotypic data.xlsx"`.

<br>

### Genotype files
AbxRxPro accepts certain output files from RGI, staramr and amrfinder. These file types are directly outputted from the Galaxy workflow provided in the AbxRxPro package. The accepted output file for each programme is:

Source | output file
:---: | :---
RGI | summary.txt (tabular)
staramr | summary.tsv (tabular)
amrfinder | result.tabular

<br>

> ### Note
> Genotypic data files are optional, only the phenotype file is required for a profile. But once a profile is saved without genotypic information, it cannot be added to the profile later. 

<br>

Whilst the workflow and antibiotic resistance genes identifying tools themselves allow multiple isolates to be processed simultaneously, it is essential that each isolate/genome is processed singularly. Thus, each output file should be tied to a single isolate ID. This is because AbxRxPro requires the file names to be in the format of isolateID_source.format in order to be detected. The files can be in either tsv, tabular or text format but requires the information to be tab seperated (which is the default for the output files).

<br>

>### Important: 
>The file names for genotypic data must state the isolate ID and the source of genotype data (RGI, staramr & amrfinder). File names are case sensitive.<br>

```
isolateID_RGI.tabular
isolateID_staramr.tabular
isolateID_amrfinder.tabular
```

<br>

## Create a profile
Instead of providing lengthy commandline arguments each time to load the same antibiotic resistance profile, you can save and load profiles.
<br><br>
To save a profile and display it:
```
abxrxpro -n profile_name -P file\path\pheno.xlsx
```
Or save a profile without displaying it:
```
abxrxpro -b profile_name -P file\path\pheno.xlsx
```
<br>

To load a saved profile:
```
abxrxpro -l profile_name
```
To see all your saved profiles as a summary:
```
abxrxpro -p
```

<br>

## Export a profile
The <em>Plotly</em> module allows for graph objects (figures) to be exported as HTML files which can be launched into your web browser. AbxRxPro uses this module to export profiles as HTML files which you can do with the `-e` flag. This can be useful for retrieving an antibiotic resistance profile without interacting with AbxRxPro. It is recommended to export a profile only when saving or loading a profile, otherwise your HTML file will be saved with the name False.
<br><br>
To export a profile when saving (`-b` can be used instead of `-n`):
```
abxrxpro -P file\path\pheno.xlsx -n profile_name -e
```
To export a profile when loading:
```
abxrxpro -l profile_name -e
```

<br>

## Input arguments

Short flag | Long flag | Argument | Description
--- | --- | --- | --- 
`-P` | `--pheno` | file path | Required for building a profile. Layout of Excel file is defined under <i>Input Files</i> section.
`-R` | `--RGI` | folder path | See <i>Input Files</i> for details.
`-S` | `--staramr` | folder path | See <i>Input Files</i> for details.
`-A` | `--amrfinder` | folder path | See <i>Input Files</i> for details.
`-a` | `--antibiotics` | space seperated list | Specify antibiotics to be displayed in plot. Default are those provided in phenotype file.
`-l` | `--load` | profile name | Load a saved antibiotic resistance profile
`-c` | `--colours` | space seperated list | Personalise colour scheme with a list of RGB values in order of <br>resistant, intermediate, susceptible and undetermined. Default values:<br>R : (255, 65, 54), <br>I : (255, 144, 14), <br>S : (44, 160, 101), <br>U : (93, 164, 214). <br>Usage: `abxrxpro.py -c (1,2,3) (4,5,6) (7,8,9) (10,11,12)`
`-d` | `--delete` | profile name | Delete a saved profile
`-p` | `--profiles` |  | List saved profiles
`-e` | `--export` |  | Export a saved profile when loading a saved profile or creating a new one. <br>Usage:<br>`abxrxpro.py -l profile_name -e` <br>or <br>`abxrxpro -P file -b profile_name -e`
&nbsp; | `--hide` |  | Do not show the gene frequency plot when genotype files are given with input. <br>Default is show gene frequency plot. Note: gene frequency plot is not displayed if no genotypic data is available.
 `-b` | `--build` | profile name | Create and save an antibiotic resistance profile without plotting. 
 `-n` | `--new` | profile name | Create, save and plot a new antibiotic resistance profile.
&nbsp; | `--find_log` |  | Identify where log files are being stored so can be sent as part of a bug report.
 


<br>

## Help
To load the help page, run the programme with no arguments or run with the `-h` or `--help` flag:
```
abxrxpro -h
```
If you are still experiencing problems, submit a bug report either in the Github page or email the author with a log file. Use `--find_log` flag will tell you where to find the log files. Including a log file as part of your bug report helps with debugging and understanding where things are going wrong.
