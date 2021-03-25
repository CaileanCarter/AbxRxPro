"""
AbxRxPro: Antibiotic Resistance Profiler

Version:                    2.1.1-alpha
Last modified:              25/03/2021
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

import pandas as pd 
import json


class abxcorr:

    def __init__(self, data, pheno):
        self.data = data
        self.phen = pheno.replace({"R" : 1, "S" : 0, "U" : 0, "I" : 0.5})
        self.class_antibiotic = json.load(open("settings.json"))["antibiotics"]


    def __call__(self):
        antibiotic_gene_pair = self.get_antibiotic_gene_pairs()

        for antib, genes in antibiotic_gene_pair.items():
            for gene in genes:

                a = self.get_abx(gene, antib)
                r = self.fetch_results(antib, gene)

                for b in a:
                    if b in self.phen.columns:
                        cor = self.phen[b].corr(r)
                        st = self.strength(cor)
                        if st != "No relationship":
                            print(gene, b, cor, st)


    def get_antibiotic_gene_pairs(self):
        antibiotic_gene_pair = {}

        for value in self.data.values():
            for antibiotic, genes in value.items():

                if antibiotic in antibiotic_gene_pair.keys():
                    antibiotic_gene_pair[antibiotic] += genes
                else:
                    antibiotic_gene_pair[antibiotic] = genes

        for antibiotic, genes in antibiotic_gene_pair.items():
            antibiotic_gene_pair[antibiotic] = set(genes)

        return antibiotic_gene_pair


    def get_abx(self, gene, antibiotic):
        cl = True if antibiotic in self.class_antibiotic.values() else False

        if cl:
            antibiotics = []
            for abx, x in self.class_antibiotic.items():
                if x == antibiotic and abx in self.phen.columns:
                    antibiotics.append(abx)
            return antibiotics

        else:
            return [antibiotic]


    def fetch_results(self, antibiotic, gene):
        result = {sample : 0 for sample in self.data.keys()}
        for sample, values in self.data.items():
            if antibiotic in values.keys() and gene in values[antibiotic]:
                result[sample] += 1
        r = pd.Series(result)
        return r


    def strength(self, x):

        if x >= 0.7:
            return "very strong positive relationship"

        elif 0.69 >= x >= 0.4:
            return "strong positive relationship"
        
        elif 0.39 >= x >= 0.3:
            return "moderate positive relationship"

        elif 0.29 >= x >= 0.2:
            return "weak positive relationship"

        elif 0.19 >= x >= 0.01:
            return "No or negligible relationship"

        elif -0.01 >= x >= -0.19:
            return "No or neglibigle relationship"

        elif -0.2 >= x >= -0.29:
            return "weak negative relationship"

        elif -0.3 >= x >= -0.39:
            return "moderate negative relationship"

        elif -0.4 >= x >= -0.69:
            return "strong negative relationship"

        elif -0.7 >= x:
            return "very strong negative relationship"

        else:
            return "No relationship"


    
