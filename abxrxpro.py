#!/usr/bin/python3

"""
AbxRxPro: Antibiotic Resistance Profiler

Version:                    2.0.1-alpha
Last modified:              04/12/2020
Github:                     https://github.com/CaileanCarter/AbxRxPro
Author:                     Cailean Carter
Email:                      cailean.carter@quadram.ac.uk
Institute affiliation:      Quadram Institute, Norwich Research Park


AbxRxPro is a tool for the visualisation of phenotypic antibiotic resistance as a bubble plot.
It allows for the inclusion of genotypic data from major antibiotic resistance genotype identifying programmes
like RGI (resistance gene identifier), staramr and amrfinder. Gene frequencies are also plotted alongside.
Plots can be saved as profiles which can be loaded directly. Plots can also be exported as an interactive
HTML file. Plots include text overlay and user-defined colour scheme for the plot.

Plotting is done with Plotly and the plot will be displayed in your browser.

For help, run AbxRxPro without any arguments or the -h or --help flag.
Alternatively, check out the AbxRxPro documentation.
"""

import argparse
import glob
import json
import logging
# import re

from itertools import product
from datetime import datetime
from os import path, remove

import plotly.graph_objects as go
import plotly.io as pio
import pandas as pd


relativepath = path.dirname(path.abspath(__file__))

date = datetime.now()

logging.basicConfig(
    filename = path.join(relativepath, "logs", f"{date.strftime('%y-%m-%d')}.log"),
    format = '%(name)s - %(levelname)s - %(message)s', 
    level = logging.INFO
    )


### BASIC FUNCTIONS ###

def delete_profile(profile):

    logging.info(f"""{date.strftime("%d/%m/%Y %H:%M")}   Starting new log. \nUser requested to delete profile: {profile}""")

    # Check user wants to delete this profile
    response = input(f"Are you sure you want to delete {profile}?: (y/n)\n")
    # Handle response
    if response == "n":
        return
    elif response == "y":
        logging.info("User confirmed deletion of profile")
        print(f"Deleting {profile}...")
    else:
        print("Sorry, response was not recognised. Please try again and respond with 'y' or 'n'.")
        logging.warning(f"User did not provide a valid response. User entered: {response}")
        return

    with open(rf"{relativepath}\settings.json", 'r') as settings:
        data = json.load(settings)

    try:
        del(data['profiles'][profile])
        logging.info("Deleted profile from imported settings.")

    except:
        logging.error(f"Profile ({profile}) wasn't found and exception raised.")
        raise ValueError("Profile not found. Please enter a valid profile ID. See profile list for available profiles.")

    with open(rf"{relativepath}\settings.json", 'w') as settings:
        logging.info("Settings with updated profile list has been saved.")
        json.dump(data, settings, indent=4)

    try:
        remove(rf"{relativepath}\profiles\{profile}.json")
    except FileNotFoundError:
        logging.error(f"Profile files weren't found for {profile}. Exception raised.")
        raise FileNotFoundError(f"Could not find profile file associated with {profile}")

    try:
        remove(rf"{relativepath}\profiles\{profile}_gf.json")
    except FileNotFoundError:
        logging.warning(f"No gene frequency file could be found for {profile}. Perhaps no genotype data available.")

    logging.info("""Profile associated files were successfully removed. \nDone. \n \n""")


def list_profiles():
    "List all saved antibiotic resistance profiles."

    logging.info(f"""{date.strftime("%d/%m/%Y %H:%M")}   Starting new log. \nUser requested to see profiles...""")

    with open(rf"{relativepath}\settings.json", 'r') as settings:
        data = json.load(settings)
    profiles = data["profiles"]
    profiletable = pd.DataFrame.from_dict(profiles, orient="index")
    print(profiletable.to_string())
    logging.info("""Profiles were successfully displayed. \nEnd. \n \n""")
    parser.exit()


def _find_log():
    logging.info("User is requesting location of log file...")
    log = path.join(relativepath, "logs", date.strftime('%y-%m-%d'))
    print(f"Location of relevant log file: \n{log}")
    

### HANDLES INPUT ###
class DataHandler:

    "Takes user inputs relating to data files."

    def __init__(self, antibiotics=None, profile=None, build=None, export=False, relativepath=None, show_genefrequency=False):
        
        # user input
        self.build = build
        self.export = export
        self.profile = profile
        self.antibiotics = [antibiotic.capitalize() for antibiotic in antibiotics] if antibiotics else None
        self.show_genefrequency = show_genefrequency 

        # internal variables
        self.Colours = {
            "R" : 'rgb(238, 102, 119)',
            "I" : 'rgb(204, 187, 68)',
            "S" : 'rgb(102, 204, 238)',
            "U" : 'rgb(68, 119, 170)'
            }
        # self.Colours = {
        #     "R" : 'rgb(255, 65, 54)',
        #     "I" : 'rgb(255, 144, 14)',
        #     "S" : 'rgb(44, 160, 101)',
        #     "U" : 'rgb(93, 164, 214)'
        #     } 
        self.relativepath = relativepath if relativepath else path.dirname(path.abspath(__file__))        
        self.ClassSelection = {}
        self.Genotype = False

        # Data for plotting
        self.isolateIDs = []
        self.GeneFrequencies = {}
        self.antibioticcolours = []
        self.GeneCount = []
        self.textoverlay = []
        self.AbxAnnots = []
        self.data = {}


    def __call__(self, pheno=None, RGI=None, staramr=None, amrfinder=None, colours=None, dev_1=False, dev_2=False):

        logging.info(
f"""{date.strftime("%d/%m/%Y %H:%M")}   Starting new log.

User input parameters:

Phenotype file:     {pheno} 
RGI folder:         {RGI} 
staramr folder:     {staramr}
amrfinder folder:   {amrfinder}
Antibiotics:        {self.antibiotics}
Profile name:       {self.profile}
Build?              {self.build}
Export?             {self.export}

""")

        if colours:
            logging.info("Assigning colours")
            self.set_colours(colours)

        if self.build:
            self.profile = self.build

        if not pheno:
            logging.error("No phenotypic data file was given.")
            raise FileNotFoundError("No phenotypic data file given")

        logging.info("Getting phenotypic data.")
        self.get_pheno(pheno)

        if not self.antibiotics:
            logging.error(f"User did not provide any antibiotics to plot. \nAntibiotics variable: {self.antibiotics}")
            raise ValueError("No list of antibiotics provided in phenotype file or using -a flag. \nCheck documentation for how antibiotics are declared.")

        logging.info(f"Isolate IDs identified:  {self.isolateIDs}")

        logging.info("Getting antibiotics/classes from settings.")
        self.load_antibiotic_settings()
        self.antibiotics.sort()
        self.data = {IsolateID : {} for IsolateID in self.isolateIDs}

        if any([RGI, staramr, amrfinder]):
            logging.info("Genotype files detected")
            self.Genotype = True
        else:
            logging.info("No genotype files given")
            self.show_genefrequency = True

        if RGI:
            logging.info("Getting RGI data...")
            self.get_RGI(RGI)

        if staramr:
            logging.info("Getting staramr data...")
            self.get_staramr(staramr)

        if amrfinder:
            logging.info("Getting amrfinder data...")
            self.get_amrfinder(amrfinder)

        if dev_1:
            logging.warning("Exporting self.data; this is a development tool")
            self._output_data()

        if dev_2:
            logging.warning("Exporting self.GeneFrequencies; this is a development tool")
            self._output_genefrequencies()

        logging.info("Creating plot annotations.")
        self.make_annotations()
        

    def load_antibiotic_settings(self):
        with open(rf"{self.relativepath}\settings.json", 'r') as settings:
            data = json.load(settings)
        antibioticClass = data["antibiotics"]
        for antibiotic in self.antibiotics:
            Class = antibioticClass.get(antibiotic)
            if Class:
                self.ClassSelection.update({antibiotic : Class})
        logging.info("Retrieved antibiotic classes.")
             

    def set_colours(self, userScheme):
        NewScheme = {}
        for status, rgb in zip(self.Colours.keys(), userScheme):
            NewScheme[status] = f"rgb{rgb}"
        self.Colours = NewScheme
        print("User defined colour scheme has been set.")


    def get_pheno(self, pheno):

        print("Reading antibiotic resistance phenotypic data...")

        pheno_df = pd.read_excel(pheno, index_col=0)
        pheno_df = pheno_df.sort_index(axis=1)
        pheno_df.columns = map(str.capitalize, pheno_df.columns)

        #Antibiotics on header
        if not self.antibiotics:
            self.antibiotics = list(pheno_df.keys())
        else:
            self.antibiotics.sort()
            if any(filter(lambda antibiotic : antibiotic not in pheno_df.columns, self.antibiotics)):
                raise ValueError(f"""Antibiotic given with -a flag which is not in the phenotypic Excel file column names. \nCheck that {self.antibiotics} are present in the Excel file column headers.""")

        self.isolateIDs = list(pheno_df.index)

        # Replace missing values with 'U' for 'undetermined' for antibiotic resistance provided by user.
        # Then get the colour pallet for given identifier.
        pheno_df = pheno_df.filter(items=self.antibiotics).fillna('U')

        # Get the text overlay which is the user's original input
        text = pheno_df.values.tolist()
        self.textoverlay = [val for sublist in text for val in sublist]

        try:
            pheno_df = pheno_df.applymap(lambda user_colour: self.Colours[user_colour]) 
        except KeyError as key:
            raise KeyError(f"Invalid phenotypic resistance value: {key}. Accepted values are: R I, S, U.")

        Colour_map = pheno_df.values.tolist()
        self.antibioticcolours = [val for sublist in Colour_map for val in sublist]

        logging.info("Phenotypic data loaded.")
        print("Done.")
    

    def assign_gene(self, ID, antibiotic, gene):

        if antibiotic in self.data[ID].keys():
            self.data[ID][antibiotic].append(gene)
        else:
            self.data[ID][antibiotic] = [gene]

        if gene in self.GeneFrequencies.keys():
            self.GeneFrequencies[gene]["isolates"].append(ID)
        else:
            self.GeneFrequencies[gene] = {"isolates" : [ID]}


    def get_staramr(self, filepath):

        print("Reading staramr data...")

        all_files = glob.glob(path.join(filepath, "*_staramr.t*"))
        logging.info("Attempting to concat files for staramr")

        try:
            staramr = pd.concat((pd.read_csv(file, sep="\t", index_col=0).assign(filename = path.basename(file)) for file in all_files))
        except ValueError:
            logging.error("Could not detect any staramr data files.")
            raise FileNotFoundError("No files detected for staramr data")

        try:
            staramr = staramr[~staramr.Genotype.str.contains("None")]                                                                   # Remove empty datasets
        except TypeError:
            logging.error("User provided resfinder file version.") 
            raise ValueError("You have provided the resfinder output of staramr instead of summary.tsv") 

        staramr["Predicted Phenotype"] = staramr["Predicted Phenotype"].str.title().str.split(", ")                                     # split the row into a list
        staramr["Genotype"] = staramr["Genotype"].str.split(", ")
        staramr["filename"] = staramr["filename"].str.replace("_staramr.t.*$", "", regex=True)
        for row in staramr.itertuples():
            for antibiotic, gene in zip(row[2], row[1]):                                                                                # iterate over all the antibiotic - gene pairs
                self.assign_gene(row[3], antibiotic, gene)                                                                              # give to the assign_gene function to sort out along with isolate ID.
        logging.info("staramr data loaded.")

        print("Done.")


    def get_RGI(self, filepath):

        print("Reading RGI data...")

        all_files = glob.glob(path.join(filepath, "*_RGI.t*"))        
        logging.info("Attempting to concat files for RGI")

        try:
            RGI = pd.concat((pd.read_csv(file, sep="\t").assign(filename = path.basename(file)) for file in all_files))     # Concat all the .tsv files and include file name which has the isolate ID
        except ValueError:
            logging.error("Could not detect any RGI files")
            raise FileNotFoundError("No files detected for RGI data")

        RGI = RGI.filter(items=["Best_Hit_ARO", "Drug Class", "filename"])
        RGI = RGI[~RGI["Drug Class"].isna()]    

        RGI["Drug Class"] = RGI["Drug Class"] \
            .str.replace(" antibiotic", "") \
            .str.title() \
            .str.split("; ")                                                                                                # Tidies the antibiotic class list
        RGI["filename"] = RGI["filename"].str.replace("_RGI.t.*$", "", regex=True)

        def filter_antibiotics(Class):
            "Checks to see if antibiotic class in one of those provided by the user."
            antibiotics = [antibiotic + 's' for antibiotic in Class]
            return list(filter(lambda x: x in self.ClassSelection.values(), antibiotics))

        RGI["Drug Class"] = RGI["Drug Class"].apply(lambda Class: filter_antibiotics(Class)) 
        RGI = RGI[RGI["Drug Class"].map(lambda x: len(x)) > 0]                                                              # remove rows with empty lists.
        for row in RGI.itertuples():
            for antibiotic in row[2]:                                                                                       # iterate over all the antibiotics
                self.assign_gene(row[3], antibiotic, row[1])                                                                # give to the assign_gene function to sort out along with isolate ID.
        
        logging.info("RGI data loaded.")
        print("Done.")


    def get_amrfinder(self, filepath):

        print("Reading amrfinder data...")

        all_files = glob.glob(path.join(filepath, "*_amrfinder.t*"))
        logging.info("Attempting to concat amrfinder files")
        try:
            amrfinder = pd.concat((pd.read_csv(file, sep="\t", index_col=0).assign(filename = path.basename(file)) for file in all_files))
        except ValueError:
            logging.error("Could not detect any amrfinder data files.")
            raise FileNotFoundError("No files detected for amrfinder data")

        amrfinder = amrfinder.filter(items=["Gene symbol", "Subclass", "filename"])                                                     # extract relevant info
        amrfinder = amrfinder[~amrfinder["Subclass"].isna()]    

        amrfinder["Subclass"] = amrfinder["Subclass"].str.capitalize()                                                                  # tidy antibiotic names
        amrfinder["filename"] = amrfinder["filename"].str.replace("_amrfinder.t.*$", "", regex=True)
        for row in amrfinder.itertuples():
            self.assign_gene(row[3], row[2], row[1])                                                                                    # give to the assign_gene function to sort out along with isolate ID.

        logging.info("amrfinder data loaded.")
        print("Done.")


    def make_annotations(self):

        print("Writing graph annotations...")

        for isolateID, antibiotic in product(self.isolateIDs, self.antibiotics):
            annotation = ""
            count = 5                          # minimum marker size

            if antibiotic in self.data[isolateID].keys():
                annotation = f"<b>{antibiotic}:</b> <br>" + "<br>".join(set(self.data[isolateID][antibiotic]))
                count += len(self.data[isolateID][antibiotic])

            Class = self.ClassSelection.get(antibiotic)
            if Class in self.data[isolateID].keys():

                if annotation:
                    annotation += "<br>"

                annotation += f"<b>{Class}:</b> <br>" + "<br>".join(set(self.data[isolateID][Class]))
                count += len(self.data[isolateID][Class])

            self.GeneCount.append(count)
            if annotation:
                self.AbxAnnots.append(annotation)
            else:
                self.AbxAnnots.append("<b>No genotype identified</b>")   
        print("Done.")    


    def _output_data(self):
        if path.exists(rf"{self.relativepath}\data.json"):
            raise FileExistsError
        with open(rf"{self.relativepath}\data.json", "w") as outputdata:
            json.dump(self.data, outputdata, indent=4)

    def _output_genefrequencies(self):
        if path.exists(rf"{self.relativepath}\genefrequencies.json"):
            raise FileExistsError
        with open(rf"{self.relativepath}\genefrequencies.json", "w") as output:
            json.dump(self.GeneFrequencies, output, indent=4)


### PLOTS PROFILES ###
class Display(DataHandler):

    "Plots antibiotic resistance files."

    def open_profile(self, profile):
        logging.info(f"""{date.strftime("%d/%m/%Y %H:%M")}   Starting new log.
        
Opening profile: {profile}""")

        print(f"Loading {profile}...")

        fig = pio.read_json(rf"{self.relativepath}\profiles\{profile}.json")                # open main profile

        try:
            gffig = pio.read_json(rf"{self.relativepath}\profiles\{profile}_gf.json")       # open gene frequency plot 
        except FileNotFoundError:
            self.show_genefrequency = True
            logging.warning("No gene frequency file associated with this profile. Possibly no genotypic data available.")

        logging.info("Profile associated files openned.")

        if self.export:
            self.export_HTML(profile, fig, gffig)                                           # export the profile if user requested 

        print(f"{profile} has loaded... \nNow plotting...")

        pio.show(fig)
        logging.info("Main plot displayed")

        if not self.show_genefrequency:
            pio.show(gffig)
            logging.info("Gene frequency plot displayed")
            logging.info("Done. \n \n")
        
        print(f"{profile} is now displaying...")
        
        
    def create_profile(self, profile, fig, gffig):

        "Create a profile from inputted data which can be launched from again instead of inputting data and processing."

        logging.info(f"Creating new profile with name: {profile}")

        with open(rf"{self.relativepath}\settings.json", 'r') as settings:
            profiles = json.load(settings)

        if profile in profiles["profiles"].keys():
            logging.error(f"Creating profile with name {profile} when it already exists.")
            raise ValueError("Profile already exists. Please provide a different name or delete existing profile.")

        time = datetime.now()
                                                                                    # Create summary of profile to be shown in -p
        profiles['profiles'][profile] = {
            "Date created" : time.strftime("%d/%m/%Y %H:%M"),
            "Antibiotic count" : len(self.antibiotics),
            "Number of isolates" : len(self.isolateIDs),
            "Antibiotics" : ', '.join(self.antibiotics[:5]) + '...',
            "Isolates" : ', '.join(self.isolateIDs[:5]) + '...'}

        with open(rf"{self.relativepath}\settings.json", 'w') as settings:
            logging.info("Saving new profile to settings.")
            json.dump(profiles, settings, indent=4)

        pio.write_json(fig, rf"{self.relativepath}\profiles\{profile}.json")
        if self.Genotype:
            pio.write_json(gffig, rf"{self.relativepath}\profiles\{profile}_gf.json")

        logging.info("Profile saved.")
        print("New profile has been saved.")
    

    def export_HTML(self, profile, fig, gffig):

        "Exports profiles to a HTML format which can be launched to a web browser and retain interactivity."
        
        logging.info("Exporting profile.")

        Download_path = path.join(path.expanduser('~'),'Downloads')

        if path.exists(path.join(Download_path, f'{profile}.html')):
            logging.warning("Profile already exported to Downloads folder.")
            raise FileExistsError("Profile already in Downloads folder.")

        Downloads = path.join(Download_path, f'{profile}.html')
        go.Figure.write_html(fig, Downloads)
        
        if self.Genotype:
            Downloads = path.join(Download_path, f'{profile}_gf.html')
            go.Figure.write_html(gffig, Downloads)

        logging.info(f"Profiles were exported to: {Download_path}\\{profile}.html & {profile}_gf.html")

        print(f"""Your antibiotic resistance profile has been exported to your Downloads folder:
{Download_path}\\{profile}.html """ + "& {profile}_gf.html" if self.Genotype else "")


    def plot_frequencies(self):

        "Plots gene frequencies"

        logging.info("Plotting gene frequencies")

        gf = pd.DataFrame.from_dict(self.GeneFrequencies, orient="index")
        gf["isolates"] = gf["isolates"].apply(lambda x : set(x))
        gf["frequency"] = gf["isolates"].apply(lambda x : len(x))
        gf["frequency"] = gf["frequency"].divide(len(self.isolateIDs)).mul(100)
        gf["isolates"] = gf["isolates"].apply(lambda x : "<b>Found in isolates:</b><br>" + "<br>".join(x))

        gffig = go.Figure(data=[go.Bar(
            x=gf.index,
            y=gf["frequency"],
            text= gf["isolates"])])

        gffig.update_layout(
            title = "<b>AbxRxPro:</b> {profile} gene frequencies".format(profile = self.profile if self.profile else ""),
            title_font_size = 24,
            xaxis_title = "Resistance gene",
            yaxis_title = "Gene frequency (%)"
        )

        gffig.update_xaxes(tickangle=45)

        return gffig


    def plot(self):
            
        logging.info("Plotting antibiotic resistance profile.")

        x = []
        y = []
        for x_value, y_value in product(range(1, len(self.isolateIDs) + 1), range(1, len(self.antibiotics) + 1)):
            x.append(x_value), y.append(y_value)

        # Plotting data
        fig = go.Figure(data=[go.Scatter(
            x=x, y=y,
            text = self.AbxAnnots,
            mode = 'markers',
            customdata = self.textoverlay,
            texttemplate = "%{customdata}",
            textposition = "middle center",
            textfont = dict(size=16),
            showlegend = False,
            marker = dict(
                color = self.antibioticcolours,
                size = self.GeneCount,
                sizemode ='area',
                sizeref = 2.*max(self.GeneCount)/(60.**2),
            )
        )])

        fig.update_layout(
            title = "<b>AbxRxPro:</b> {profile} antibiotic resistance profile".format(profile = self.profile if self.profile else ""),
            title_font_size = 24,
            xaxis_title = "Isolate ID",
            xaxis = dict(
                tickmode = 'array',
                tickvals = [x for x in range(1, len(self.isolateIDs) + 1)],
                ticktext = self.isolateIDs
            ),
            yaxis_title = "Antibiotics",
            yaxis = dict(
                tickmode = 'array',
                tickvals = [x for x in range(1, len(self.antibiotics) + 1)],
                ticktext = self.antibiotics
            )
        )

        fig.update_layout(
                    updatemenus = [
                        dict(
                            type = "buttons",
                            x = 1.1,
                            y = 1.1,
                            buttons = list([
                                dict(
                                    args=[{"mode" : "markers"}],
                                    args2=[{"mode" : "markers+text"}],
                                    label = "Text overlay",
                                    method = "update")]
                            ))])

        fig.add_trace(go.Scatter(
            x = [None], 
            y = [None],
            visible = "legendonly",
            name = "Resistant (R)",
            marker = dict(color = self.Colours.get("R")),
            texttemplate = "R",
            textposition = "middle center"
            ))
        fig.add_trace(go.Scatter(
            x = [None], 
            y = [None],
            visible = "legendonly",
            name = "Intermediate (I)",
            marker = dict(color = self.Colours.get("I")),
            texttemplate = "I",
            textposition = "middle center"
            ))
        fig.add_trace(go.Scatter(
            x = [None], 
            y = [None],
            visible = "legendonly",
            name = "Susceptible (S)",
            marker = dict(color = self.Colours.get("S")),
            texttemplate = "S",
            textposition = "middle center"
            ))
        fig.add_trace(go.Scatter(
            x = [None], 
            y = [None],
            visible = "legendonly",
            name = "Undetermined (U)",
            marker = dict(color = self.Colours.get("U")),
            texttemplate = "U",
            textposition = "middle center"
            ))
        fig.update_layout(legend_font_size=16)

        if self.Genotype:
            gffig = self.plot_frequencies()                     # create gene frequency plot
        else:
            gffig = None

        if self.profile:
            self.create_profile(self.profile, fig, gffig)   # plot the given profile name

        if self.export:
            self.export_HTML(self.profile, fig, gffig)      # export to HTML

        if not self.build:
            print("Loading plot...")
            fig.show()                                      # plot the graphs...
            print("Now plotting your antibiotic resistance profile...")

            if not self.show_genefrequency:
                gffig.show()                                # plot the gene frequencies
                print("Gene frequencies now loaded.")

        logging.info("End. \n \n")


#### Argument parsing from commandline ####

parser = argparse.ArgumentParser(
description="Create or load an antibiotic resistance profile.", formatter_class=argparse.RawTextHelpFormatter)

Paths = parser.add_argument_group(title="File/folder inputs:")

Paths.add_argument('-P', '--pheno', metavar="Excel file path", 
help="""Required for building a profile.
File path to Excel file containing isolate names in first column.
If phenotypic resistance is included, header of each column should be antibiotic name.
Indicate phenotype with initialisations:
R : resistant
I : intermediate
S : susceptible
U : undetermined (optional: can be left blank)""")

Paths.add_argument('-R', '--RGI', metavar="Folder path", 
help="""Path to folder with RGI summary files. Ensure each file name contains the 
isolate name in the format of IsolateID_RGI. Accepted formats are: TSV, tabular and text.
Values must be tab seperated.""")

Paths.add_argument('-S', '--staramr', metavar="Folder path", 
help="""Path to folder with staramr files. Ensure each file name contains the isolate name in 
the format of IsolateID_staramr. Accepted formats are: TSV, tabular and text.
Values must be tab seperated.""")

Paths.add_argument('-A', '--amrfinder', metavar="Folder path", 
help="""Path to folder with amrfinder files. Ensure each file name contains the isolate name in
 the format of IsolateID_amrfinder. Accepted formats are: TSV, tabular and text.
Values must be tab seperated.""")

#provide own list of antibiotics
parser.add_argument('-a', '--antibiotics', nargs='+', metavar="name", 
help="""Select which antibiotics from phenotypic antibiotic resistance file are to be displayed (space seperated).
Default is all the antibiotics given in column headers of phenotype file.
Usage: abxrxpro.py -a Ampicillin Tetracycline Trimethoprim""")

# load a profile from a list
parser.add_argument('-l', '--load', metavar="profile name", help="Load antibiotic resistance profile.")

# change colours of plot
parser.add_argument('-c', '--colours', nargs=4, metavar=("(RGB)", "(RGB)", "(RGB)", "(RGB)"), 
help="""Specify colour indications for phenotypic resistance.
Provide space seperated list of RGB colours in order of
resistant, intermediate, susceptible, undetermined.
Default values:
R : (238,102,119),
I : (204,187,68),
S : (102,204,238),
U : (68,119,170)
Usage: abxrxpro.py -c (1,2,3) (4,5,6) (7,8,9) (10,11,12)""") 

parser.add_argument('-d', '--delete', type=str, metavar='profile name', help="Delete a profile.")
parser.add_argument('-p', '--profiles', action="store_true", default=False, help="List antibiotic resistance profiles.")
parser.add_argument('-e', '--export', action="store_true", default=False, help="Export profile to HTML")
parser.add_argument('--hide', action="store_true", default=False, help="Do not show gene frequencies.")
parser.add_argument('-D1', action="store_true", default=False, help=argparse.SUPPRESS)         # Exports the main data container self.data (dev tool)
parser.add_argument('-D2', action="store_true", default=False, help=argparse.SUPPRESS)         # Exports the gene frequency container self.GeneFrequencies (dev tool)
parser.add_argument('--find_log', action="store_true", default=False, help="Move log file to Downloads folder to be submitted as part of a bug report")


group1 = parser.add_mutually_exclusive_group()
# only builds the profile and does not plot
group1.add_argument('-b', '--build', default=False, type=str, metavar="profile name", 
help="""Create and save an antibiotic resistance profile without plotting.
Requires a name for the profile as argument.""")
# saves profile and loads plot
group1.add_argument('-n', '--new', default=False, type=str, metavar="profile name", 
help="""Create, save and plot a new antibiotic resistance profile.
Requires a name for the profile as argument.""")

args = parser.parse_args()


### HANDLE INPUT ###

if args.profiles:
    list_profiles()                                 
    parser.exit()

if args.delete:
    delete_profile(args.delete)                     # Delete profile
    parser.exit()

# if a file was given without phenotypic file as reference:
if not args.pheno and any([args.RGI, args.staramr, args.amrfinder]):
    logging.error(
f"""{date.strftime("%d/%m/%Y %H:%M")}   Starting new log.

User provided a genotype data file without a phenotypic file:

Pheno:      {args.pheno}

RGI:        {args.RGI}
staramr:    {args.staramr}
amrfinder:  {args.amrfinder}

End. \n \n"""
)
    raise FileNotFoundError("""A file was given without phenotypic file for reference of isolate IDs.
Please provide an Excel file containing at least 1 column with isolate names.
Use -P flag to give file path for this file.""")

# if no input is given, the help message is shown
if not any(vars(args).values()):
    print(__doc__)
    parser.print_help()
    parser.exit()

if args.find_log:
    _find_log()
    parser.exit()

# if a profile is given, the plot will be launched
if args.load:
    ShowProfile = Display(export=args.export, relativepath=relativepath, show_genefrequency=args.hide)
    ShowProfile.open_profile(args.load)

else: # if plotting a new profile
    ShowProfile = Display(args.antibiotics, 
        profile=args.new, 
        build=args.build, 
        export=args.export, 
        relativepath=relativepath, 
        show_genefrequency=args.hide)

    ShowProfile(args.pheno, 
        args.RGI, 
        args.staramr, 
        args.amrfinder, 
        colours=args.colours,
        dev_1 = args.D1,
        dev_2 = args.D2)

    ShowProfile.plot()

print("\nThank you for using AbxRxPro.")