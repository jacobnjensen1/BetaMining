#!/usr/bin/env python3
#
# This script provides an interface that combines the functions in the
# "BetaMining" package to analyze .pdb files of interest.
import pandas as pd

#from biopandas.pdb import PandasPdb
#import datetime
import os
import re
import json
from glob import glob
#from pathlib import Path
import importlib.resources
import multiprocessing as mp

#import beta_mining
#from beta_mining import beta_mining_algorithm
from beta_mining import beta_mining_functions

def analyze_structure(filename, config, json, output_dictionary, thread_id):
    """Main beta_mining algorithm.

    Keyword arguments:
    filename -- the .pdb or .pdb.gz structure file being analyzed
    config -- the configuration YAML file
    json -- the json used to guide the analysis
    output_dictionary -- dictionary of what output to generate
    """
    #Note: this current setup means that any changes to the structure string persist for all targets
    #This would be an issue if there were targets with non-overlapping exclusions and inclusions
    #It might be best to call this for each target separately?

    # create ProDy model, meta information dictionary, BioPandas Dataframe and AA sequence list
    result = beta_mining_functions.create_model_metainfo(filename)
    if result != None:
      model, meta_dictionary, af_object, af_sequence = result
    else:
      return None

    # create the dataframe of dihedral angles, pandas series of secondary structures, and a dictionary with secondary structures as keys and a list of residue numbers as values.
    dihedral_df, secondary_structure_series, secondary_structure_reference = beta_mining_functions.calculation_df(model, meta_dictionary["fragment_offset"], json["secondary_structures"], config["units"])

    # create dataframe of per-residue contacts specified in the features json
    # create a dictionary of boolean Series, keys are "include" or "exclude"
    for target in json["target_region_features"]:
        if target["name"] in config["target_names"]:
            if "contacts" in list(target["include"].keys()) +  list(target["exclude"].keys()):
                contacts_dataframe, mask_dict = beta_mining_functions.contacts_df(af_object, json["target_region_features"], secondary_structure_reference, config["target_names"])
                # replace features from include with the specified symbol first
                # replace features from exclude with the specified symbol last
                for condition in config["conditions"]:
                    for series_symbol in mask_dict[condition]:
                        #print(sum(series_symbol[0] == True))
                        #print(series_symbol[0].to_string())
                        #print(series_symbol[1])
                        #print(mask_dict)
                        #print(secondary_structure_series)
                        secondary_structure_series.mask(series_symbol[0], series_symbol[1], inplace = True)
                        #print(secondary_structure_series)
                # left join the calculated dihedral and twist dataframe and the contacts dataframe
                # create a secondary structure symbols string to search with regex
                dihedral_df = pd.merge(dihedral_df, contacts_dataframe, on = "residue_number", how = "left")
            if "out_of_plane_sheets" in target["exclude"].keys():
                series_symbol = beta_mining_functions.handle_out_of_plane(model, target["exclude"]["out_of_plane_sheets"], filename, config["output_filepath"], thread_id)
                secondary_structure_series.mask(series_symbol[0], series_symbol[1], inplace = True)

    structure_symbols_string = "".join(secondary_structure_series)
    structure_length = len(structure_symbols_string)

    #print(structure_symbols_string)

    # generate dataframe of entire protein
    residue_df = beta_mining_functions.polymer_df(meta_dictionary, af_object)

    residue_df = pd.merge(residue_df, dihedral_df, on = "residue_number", how = "left")

    #print(residue_df.columns)
    #print(residue_df["secondary_structure"])
    #print(dihedral_df)
    #add plddt column from residue_df to dihedral_df
    dihedral_df = dihedral_df.join(residue_df["plddt"])
    #print(dihedral_df.to_string())

    #one per task, won't cause collisions
    residue_df["structure_symbol"] = secondary_structure_series
    if output_dictionary["proteome_aa"] != False: # save complete protein dataframe if indicated in config YAML
        output_filename = config["output_filepath"] + config["results_prefix"] + meta_dictionary["accession"] + "_f" + str(meta_dictionary["fragment"]) + "_" + meta_dictionary["depo_date"].lower() + config["output_filename"]
        if config["per_residue_output"]:
          residue_df.to_csv(output_filename, index = False)

    output_strings = [[],[]] #can hold [csv strings], [fasta strings] but only if flags are set
    # look for regex matches in the secondary structure string for each type of target
    print(f"Sites for {filename}:")
    for target in json["target_region_features"]:
        if target["name"] in config["target_names"]:
            compiled_regex = re.compile(target["regex"])
            #print(structure_symbols_string)
            targets_found = re.finditer(compiled_regex, structure_symbols_string)
            fasta_fields = [
                            target["name"],
                            meta_dictionary["id_code"],
                            meta_dictionary["accession"],
                            meta_dictionary["organism_scientific"],
                            meta_dictionary["organism_taxid"],
                            meta_dictionary["full_title"],
                            meta_dictionary["depo_date"],
                            "fragment " + str(meta_dictionary["fragment"])
                        ]
            fasta_flank = target["fasta_flank"]
            for i, match_obj in enumerate(targets_found):
                print(f"Site {i + 1}, residues {match_obj.span()[0] + 1}-{match_obj.span()[1] + 1}:")
                if beta_mining_functions.attribute_filter(target, match_obj, dihedral_df) == False:
                    print("attribute_filter returned False!")
                    continue
                else:
                    print("attribute_filter returned True!")
                    # these values are the 0 index of residues for slicing
                    res_idx_start = max(0, match_obj.span()[0] - fasta_flank)
                    res_idx_end = min(structure_length - 1, match_obj.span()[1] + fasta_flank)

                    # these values are the 1 index of residues for output
                    res_start = res_idx_start + 1
                    res_end = res_idx_end + 1
                    abs_res_start = meta_dictionary["fragment_offset"] + res_start
                    abs_res_end = meta_dictionary["fragment_offset"] + res_end

                    # sequences of amino acids and structures
                    aa_sequence = "".join(af_sequence[res_idx_start:res_idx_end])
                    ss_sequence = structure_symbols_string[res_idx_start:res_idx_end]

                    if output_dictionary["hits_fasta"] != False:
                        fasta_header = ">" + "|".join(fasta_fields) + "|residues " + \
                          str(res_start) + "-" + str(res_end) + "|abs_residues " + str(abs_res_start) + "-" + str(abs_res_end)
                        #output_dictionary["hits_fasta"].write(fasta_header + "\n" + aa_sequence + "\n")
                        output_strings[1].append(fasta_header + "\n" + aa_sequence + "\n")

                    if output_dictionary["hits_aa"] != False:
                        hit_line = [target["name"],
                            meta_dictionary["database"],
                            meta_dictionary["id_code"],
                            meta_dictionary["accession"],
                            meta_dictionary["full_name"],
                            meta_dictionary["full_title"],
                            meta_dictionary["depo_date"],
                            meta_dictionary["fragment"],
                            meta_dictionary["fragment_offset"],
                            meta_dictionary["organism_scientific"],
                            meta_dictionary["organism_taxid"],
                            aa_sequence,
                            ss_sequence]
                        for i, val in enumerate(hit_line):
                            if type(val) == str and "," in val:
                              hit_line[i] = f'"{val}"'
                        for field in config["additional_attributes"].keys():
                            for funct_name in config["additional_attributes"][field]:
                                funct = getattr(pd.Series, funct_name)
                                
                                attr_value = funct(residue_df[field][residue_df["residue_number"].isin([*range(match_obj.span()[0] + 1, match_obj.span()[1]+ 2)])])
                                hit_line.append(attr_value)
                        hit_line = list(map(str, hit_line))
                        #output_dictionary["hits_aa"].write(",".join(hit_line) + "\n")
                        output_strings[0].append(",".join(hit_line) + "\n")
    return output_strings

def listener(q, output_files):
  if output_files["hits_aa"]:
    hits_aa = open(output_files["hits_aa"], "a")
  if output_files["hits_fasta"]:
    hits_fasta = open(output_files["hits_fasta"], "w")
  while True:
    m = q.get()
    #print(m)
    if m == "done":
      print("done")
      if output_files["hits_aa"]:
        hits_aa.close()
      if output_files["hits_fasta"]:
        hits_fasta.close()
      break
    else:
      if output_files["hits_aa"]:
        for aa_entry in m[0]:
          hits_aa.write(aa_entry)
      if output_files["hits_fasta"]:
        for fasta_entry in m[1]:
          hits_fasta.write(fasta_entry)
  

#needs all the params for analyze_structure and the queue
def runner(file, config_settings, target_parameters, bool_output_files, q):
  #thread_id includes results_prefix so that multiple runs can happen at the same time in the same directory
  thread_id = f"{config_settings['results_prefix']}{mp.current_process()._identity[0]}"
  res = analyze_structure(file, config_settings, target_parameters, bool_output_files, thread_id)
  if len(res[0]) + len(res[1]) > 0:
    #print(res)
    q.put(res)

def main(config_settings):
    """Filehandling algorithm.

    Keyword arguments:
    config_dictionary -- the dictionary built from the config YAML file given to the algorithm.
    """
    ### set global run variables from config YAML ###
    input_path = config_settings["input_filepath"]
    output_path = config_settings["output_filepath"]
    results_prefix = config_settings["results_prefix"]
    output_file = config_settings["output_filename"]
    conditions = config_settings["conditions"]

    # output fields for hits_aa .csv
    field_list = ["target",
        "database",
        "id_code",
        "accession",
        "full_name",
        "full_title",
        "depo_date",
        "fragment",
        "fragment_offset",
        "organism_scientific",
        "organism_taxid",
        "aa_seq",
        "ss_seq"]
    for field in config_settings["additional_attributes"].keys():
        for function in config_settings["additional_attributes"][field]:
            field_list.append(field + "_" + function)

    ## import parameter json ##
    if config_settings["custom_json"] == None:
        with importlib.resources.open_text("beta_mining", "structure_dictionary.json") as param_json:
            target_parameters = json.load(param_json)
    else:
        with open(config_settings["custom_json"], "r") as param_json:
            target_parameters = json.load(param_json)

    ### prepare output files ###
    #hits_fasta
    if "hits_fasta" in config_settings["save_files"]:
        hits_fasta_fname = output_path + results_prefix + output_file + ".fasta"
        #hits_fasta = open(hits_fasta_fname, "w")
    else:
        hits_fasta_fname = False
    #hits_aa
    if "hits_aa" in config_settings["save_files"]:
        hits_aa_fname = output_path + results_prefix + output_file + ".csv"
        with open(hits_aa_fname, "w") as hits_aa:
          hits_aa.write(",".join(field_list) + "\n")
    else:
        hits_aa_fname = False
    #proteome_aa files are generated in analyze_structure function
    if "proteome_aa" in config_settings["save_files"]:
        proteome_aa = True
    else:
        proteome_aa = False

    output_files = {"hits_fasta": hits_fasta_fname, "hits_aa": hits_aa_fname, "proteome_aa": proteome_aa}
    bool_output_files = {key:(True if value else False) for key,value in output_files.items()}

    ### read in files and process different formats ###
    if os.path.isdir(output_path) == False:
        os.mkdir(output_path)
    file_list = glob(input_path + "/**.pdb*", recursive = True)
    #print(file_list)
    file_number = 1
    file_total = len(file_list)
    # we will sort the file list so that protein fragments are together,
    # but they won't be in the perfect order because of alphabetic sorting.
    #for file in sorted(file_list):
    #    print("Mining file " + file + ": "+ str(file_number) + " of " + str(file_total) + "...")
    #    analyze_structure(file, config_settings, target_parameters, output_files)
    #    file_number = file_number + 1
    manager = mp.Manager()
    q = manager.Queue()
    p = mp.Pool(int(config_settings["processes"]) + 1) # This value needs to be > 1, slurm ensures it won't use too many cpus

    watcher = p.apply_async(listener, (q,output_files))

    jobs = []
    for file in sorted(file_list):
      job = p.apply_async(runner, (file, config_settings, target_parameters, bool_output_files, q))
      jobs.append(job)

    for i, job in enumerate(jobs):
      print("Mining file " + str(i + 1) + " of " + str(file_total) + "...")
      job.get()

    q.put("done")
    p.close()
    p.join()

if __name__ == "__main__":
    main(config_settings)
