#!/usr/bin/env python3
#
# This script gives a command-line interface to the `beta_mining` python package
# to analyze and log processing of amino acid PDB structure files.

import os
import sys
#sys.path.insert(0,"../")
#sys.path.insert(0, "../beta_mining/")
import pandas as pd
import pkgutil

import importlib.resources
import yaml
import shutil
import signal
import datetime
from argparse import ArgumentParser

from beta_mining import beta_mining_algorithm
from beta_mining import beta_mining_functions

script_description = (
  "This is a Python 3 package for identifying secondary structure patterns from multimer assemblies of predicted monomer structure files."
)

def handler(signal_number, frame):
    """Shuts down the script via `Ctrl-C` or `Ctrl-Z` entered on the command line."""
    if signal_number == signal.SIGINT:
        print("Termination requested, aborting script.")
        sys.exit(0)
    elif signal_number == signal.SIGTSTP:
        print("Suspend not supported, aborting script.")
        sys.exit(0)

def generate_yaml(target_directory = os.getcwd(), prefix = "generated_"):
    """Creates a yaml config file in the current working directory.

    Keyword arguments:
    target_directory -- directory where the YAML config file will be created
    prefix -- the prefix for the created YAML config file
    """
    #config_default_text = pkgutil.get_data("beta_mining","default_config.yaml")
    #print(target_directory)
    #configfile =  open(target_directory + prefix + "config.yaml", "w")
    #configfile.write(config_default_text.decode("utf-8"))
    #configfile.close()
    with importlib.resources.path("beta_mining", "default_config.yaml") as p:
        shutil.copy(p, target_directory)

def main():
    parser = ArgumentParser("beta_mining", description = script_description)
    parser.add_argument("-c", "--config", help = "The config file (YAML format) to use for a run.", type = str, default = "default_config.yaml")
    #parser.add_argument("-f", "--filepath", help = "If running beta_mining in default mode, indicate path to .pdb files here. A default YAML will be generated and a timestamp-based prefix will be used for all outputs.", type = str, default = "no_input")
    parser.add_argument("-d", "--defaults", help = "Generate a default YAML config file in the current working directory.", action = "store_true")
    args = parser.parse_args()

    now = datetime.datetime.now()
    timestamp = "".join([str(now.year), str(now.month), str(now.day), str(now.hour), str(now.minute), str(now.second), "_"])

    # Make default YAML config and quit.
    if args.defaults:
        path = os.getcwd()
        generate_yaml()
        print("Default config YAML generated and saved to: '{}'. Exiting.".format(path))
        sys.exit(0)

    # If a filepath is specified, generate default YAML and run the full program.
    #elif args.filepath != "no_input":
    #    if args.config != "default_config.yaml":
    #        print("Filepath default mode and YAML are mutually exclusive. Please use one or the other. Exiting.")
    #        sys.exit(0)
    #    filepath = args.filepath

    #    generate_yaml(filepath, timestamp)
    #    print(os.getcwd())
    #    config = open(filepath + timestamp + "config.yaml", "r")
    #    settings = yaml.load(config, Loader = yaml.FullLoader)
    #    config.close()
    #    settings["results_prefix"] = timestamp
    #    settings["input_filepath"] = filepath

    #    yaml_stream = open(filepath + timestamp + "config.yaml", "w")
    #    yaml.dump(settings, yaml_stream)
    #    yaml_stream.close()

    # If a config is specified, load all arguments from it.
    if os.path.isfile(args.config):
        config = open(args.config, "r")
        settings = yaml.load(config, Loader = yaml.FullLoader)
        config.close()
        if settings["results_prefix"] == "default":
            settings["results_prefix"] = timestamp
            config = open(args.config, "a")
            config.write("\ntimestamp: " + timestamp)
            config.close()
    else:
        print("Unable to open config file, exiting.")
        sys.exit(0)

    # Create output directory if it doesn't exist
    if os.path.isdir(settings["output_filepath"]) == False:
        os.mkdir(settings["output_filepath"])

    # Run the main algorithm
    beta_mining_algorithm.main(settings)

if __name__ == "__main__":
    # Register Ctrl-C, Ctrl-Z signals
    signal.signal(signal.SIGINT, handler) # Ctrl-C
    signal.signal(signal.SIGTSTP, handler) # Ctrl-Z
    main()
