# coding: utf-8


# this script generates pileup data from a variant and sequencing data using the pileup_on_1_pos.sh script relying on samtools
# is is used to merge data across samples and can export tables and plots
# it is centerd on the variant class that contains functions to generate and analyse the pileup data
# it is best used with a parser that reads a file containing multiples variants and generates a variant object for each variant

import pandas as pd
import os
from collections import defaultdict
from pyparsing import *
from typing import Dict
from collections import defaultdict
import pyparsing as pp
from pyparsing import pyparsing_common
import matplotlib.pyplot as plt
import subprocess

class variant:
    def __init__(self, var_id=None, pileup_chr=None, pileup_pos=None, pileup_alt=None, proband_id=None, father_id=None, mother_id=None, pileup_validated=False, bam_list=None, pileup_file=None, script_generate_1_pileup=None, ref_genome=None):
        self.var_id = var_id
        self.pileup_chr = pileup_chr
        self.pileup_pos = pileup_pos
        self.pileup_alt = pileup_alt
        self.proband_id = proband_id
        self.father_id = father_id
        self.mother_id = mother_id
        self.pileup_validated = pileup_validated
        self.pileup_file = pileup_file
        self.bam_list = bam_list
        self.pileup_dict = None
        self.script_generate_1_pileup = script_generate_1_pileup
        self.ref_genome = ref_genome
        self.generate_pileup_file()
        self.read_pileup()
        self.compute_alt_count_and_AB()
        self.complete_info_in_pileup_dicts()
        # self.compute_pileup_validation_status()

    def generate_pileup_file(self):
        # if pileup file is provided, return message
        if self.pileup_file != None:
            print("pileup file already provided")
            return False
        # if bam_list is not provided, return message
        if self.bam_list == None:
            print("bam list not provided")
            return False
        # if pileup file is not provided, generate it using pileup_on_1_pos.sh
        # create output folder if it does not exist
        output_dir="output/pileups/"
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        # create the pileup file
        command = [self.script_generate_1_pileup, "--ref", self.ref_genome, "--chr", self.pileup_chr, "--pos", self.pileup_pos, "--bam_list", self.bam_list, "--output_folder", output_dir]
        subprocess.run(command)
        # produces a file chr-pos.pileup
        filename=output_dir + self.pileup_chr + "-" + self.pileup_pos + ".pileup"
        print("***************** Pileup exported as : " + filename)
        self.pileup_file = filename

    
    def read_pileup(self):
        # if pileup file is not provided, return message   
        if self.pileup_file == None:
            print("pileup file not provided")
            return False
        # lire le fichier pileup
        # returns self.pileup_dict with {chr, pos, ref, pileups (list of dict with sample, count, value, quality)}
        pileup = pd.read_csv(self.pileup_file, sep="\t")
        pileup = pileup.fillna(".")
        # le fichier est au format chr pos ref sampleA_count sampleA_value sampleA_quality ...
        # on le transforme en chr pos ref puis une liste de dictionnaires avec pour clé sample et pour valeur [count, value, quality]
        pileup_dict={}
        pileup_dict["chr"]=pileup["chr"][0]
        pileup_dict["pos"]=pileup["pos"][0]
        pileup_dict["ref"]=pileup["ref"][0]
        pileup_dict["pileups"]=[]

        # vérification que le fichier est valide
        if len(pileup) != 1:
            print("pileup file contains more than one line")
            return False
        chr_in_pileup = str(pileup["chr"][0])
        pos_in_pileup = str(pileup["pos"][0])
        chr_in_variant = str(self.pileup_chr)
        pos_in_variant = str(self.pileup_pos)
        if chr_in_pileup != chr_in_variant:
            print("pileup file chr : " + str(chr_in_pileup) + " does not match variant : " + str(chr_in_variant))
            return False
        if pos_in_pileup != pos_in_variant:
            print("pileup file pos does not match variant")
            return False

        # list all samples according to the count column, rename col to delete _count
        samples = [col for col in pileup.columns if "_count" in col]
        samples = [col.replace("_count", "") for col in samples]
        # pour chaque sample, on ajoute un dictionnaire avec les infos count, value, quality
        for sample in samples:
            sample_dict={}
            # pour le sample, on supprime le suffixe .bam ou .cram
            sample_dict["sample"]=sample.replace(".bam", "").replace(".cram", "")
            sample_dict["count"]=pileup[sample+"_count"][0]
            sample_dict["value"]=pileup[sample+"_value"][0]
            sample_dict["quality"]=pileup[sample+"_quality"][0]
            pileup_dict["pileups"].append(sample_dict)
        self.pileup_dict = pileup_dict
        return pileup_dict

    def rename_samples(self, sample_dict):
        for sample in self.pileup_dict["pileups"]:
            if sample["sample"] in sample_dict:
                sample["sample"] = sample_dict[sample["sample"]]
        return True

    def compute_alt_count_and_AB(self):
        def compute_AB(alt_count, total):
            if total == 0:
                return 0
            else:
                return alt_count / total
        # pour chaque sample, on calcule le alt_count et on ajoute le AB
        for sample in self.pileup_dict["pileups"]:
            # get_alt_count_from_pileup_line(input_str, alt_str)
            sample["alt_count"]=get_alt_count_from_pileup_line(sample["value"], self.pileup_alt)
            sample["AB"]=compute_AB(sample["alt_count"], sample["count"])


    def complete_info_in_pileup_dicts(self):
        # iterate on pileup_dict. if sample = proband_id, then "sample_status" = "case" else "sample_status" = "control"
        # if sample = father_id or sample = mother_id, then "sample_status" = "parent"
        for sample in self.pileup_dict["pileups"]:
            # remove ".sorted.groupedUniq.dedupUT" from sample name
            sample["sample"] = sample["sample"].replace(".sorted.groupedUniq.dedupUT", "")
            # adds sample_status
            if sample["sample"] == self.proband_id:
                sample["sample_status"] = "case"
            elif sample["sample"] == self.father_id:
                sample["sample_status"] = "father"
            elif sample["sample"] == self.mother_id:
                sample["sample_status"] = "mother"
            else:
                sample["sample_status"] = "control"
            # adds a string for the plot : "alt_count",  "/"",  "count"
            sample["plot_string"] = str(sample["alt_count"]) + "/" + str(sample["count"])
            # creates the tupple for the plot : (depth, AB, sample_status, plot_string)
            sample["plot_tupple"] = (sample["count"], sample["AB"], sample["sample_status"], sample["plot_string"])
        return self.pileup_dict


    def export_pileup_as_table(self):
        # generate an output folder if it does not exist
        output_dir = "output/pileup_tables/"
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        output_filename = output_dir + self.var_id + ".tsv"
        # export pileup_dicts as a table
        pileup_df = pd.DataFrame(self.pileup_dict["pileups"])
        pileup_df.to_csv(output_filename, sep="\t", index=False)
        print("***************** Pileup table exported as : " + output_filename)
        return output_filename
    

    def generate_plot_AB_vs_depth(self):
        # generate an output folder if it does not exist
        output_dir = "output/plots/"
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        # lists all tupples and launches the plot
        plot_tupples = [sample["plot_tupple"] for sample in self.pileup_dict["pileups"]]
        output_filename = output_dir + self.var_id + ".png"
        dotplot(plot_tupples, output_filename=output_filename, dpi=300, title=self.var_id, subtitle="Positive control : " + self.proband_id)
        print("***************** Plot exported as : " + output_filename)
        return output_filename
    

    def print_data(self):
        print("***************** Variant data : ")
        print("chr : " + str(self.pileup_chr))
        print("pos : " + str(self.pileup_pos))
        print("alt : " + str(self.pileup_alt))
        print("proband_id : " + str(self.proband_id))
        print("father_id : " + str(self.father_id))
        print("mother_id : " + str(self.mother_id))
        print("***************** Pileup data : ")
        print("chr : " + str(self.pileup_dict["chr"]))
        print("pos : " + str(self.pileup_dict["pos"]))
        print("ref : " + str(self.pileup_dict["ref"]))
        print("pileups : ")
        for sample in self.pileup_dict["pileups"]:
            print(sample)
        print("*****************")
    

# functions **********************************************************************************************************************


def dotplot(x_y_color_text_tuples_list, output_filename="dot.png", dpi=300, title="AB versus depth", subtitle="sub"):
    # tuple = (x, y, color, text)
    x_data = []
    y_data = []
    colors = []
    texts = []
    x_label="Depth (X)"
    y_label="Pileup AB"
    color_map = {"father": "blue", "mother": "red", "case": "green", "control": "grey"}

    for tup in x_y_color_text_tuples_list:
        if tup[0] != "" and tup[1] != "" and tup[2] != "":
            x_data.append(tup[0])
            y_data.append(tup[1])
            color = color_map.get(tup[2], "black")
            colors.append(color)
            if tup[3] != "" and tup[1] != 0:
                texts.append(tup[3])
            else:
                texts.append("")

    fig, ax = plt.subplots()
    ax.scatter(x_data, y_data, c=colors, alpha=0.5)
    for i, txt in enumerate(texts):
        if txt != "":
            ax.annotate(txt, (x_data[i], y_data[i]), textcoords="offset points", xytext=(5,5), ha="left", va="bottom", color="gray")
    if title:
        plt.title(title, fontsize=16)
    if subtitle:
        plt.suptitle(subtitle, fontsize=12)
    ax.set_xlabel(x_label, fontsize=12)
    ax.set_ylabel(y_label, fontsize=12)
    # ax.axhline(y=0.5, color='grey', linestyle='--')
    plt.savefig(output_filename, dpi=dpi)
    plt.show()


def get_alt_count_from_pileup_line(input_str: str, alt_string:str) -> Dict[str, int]:
    input_str = input_str.upper()
    ref_nucleotide = pp.oneOf(". ,")
    alt_nucleotide = pp.oneOf("A T G C a t g c *") #if ignore_case else "A T G C *")
    reference_skip = pp.oneOf("< >")
    ignored = pp.oneOf("$ ^")

    # counter is extracting how many nucleotides are to be read
    indel_counter = pp.oneOf("+ -").suppress() + pyparsing_common.integer()
    # capture the entire text matched with original_text_for after decomposing (+-)N and N*nucleotiden
    indel = pp.originalTextFor(pp.countedArray(intExpr=indel_counter, expr=alt_nucleotide))

    parser = ref_nucleotide | alt_nucleotide | reference_skip | indel
    parser.ignore(ignored)

    occurrences = defaultdict(lambda: 0)

    for m in [r[0][0] for r in parser.scanString(input_str)]:
        occurrences[m.upper()] += 1

    return occurrences[alt_string]






# example of use **********************************************************************************************************************
# in practice the variants parameters should be provided in a file and a parser should be created to read the file and generate the variant objects
# below is a manual example for a single variant

# variant
var_id="1-1737942-A-G"
pileup_chr="1"
pileup_pos="1737942" # the position that will be used to generate the pileup
pileup_alt="G" # the alt motif that will be looked for by the get_alt_count_from_pileup_line function

# sequencing data
proband_id="12-04762.sorted.dedup.real.BQSR" # sample name in the bam_list, if used, will be used to color the plot
father_id="" # sample name in the bam_list, if used, will be used to color the plot
mother_id="" # sample name in the bam_list, if used, will be used to color the plot
bam_list="all_bam.tsv"

# script and reference genome
script_generate_1_pileup="./pileup_on_1_pos.sh"
ref_genome="/storage/store-01/RunsPlateforme/Reference/Homo_sapiens/hg19/human_g1k_v37.fasta"

# create the variant
var = variant(var_id=var_id, pileup_chr=pileup_chr, pileup_pos=pileup_pos, pileup_alt=pileup_alt, proband_id=proband_id, father_id=father_id, mother_id=mother_id, bam_list=bam_list, script_generate_1_pileup=script_generate_1_pileup, ref_genome=ref_genome)

# export pileup table
var.export_pileup_as_table()

# create a plot
var.generate_plot_AB_vs_depth()
