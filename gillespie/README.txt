This folder contains scripts to generate Gillespie simulations

Steps to generate Gillespie simulations

1) Generate parameter files as mentioned in Example_run.ipynb jupyter notebook.
2) Run param_gillespie.py to execute the simulation.

usage: param_gillespie.py [-h] [-T MAX_TIME] [-I ITEM] [-L ITEM_LEVEL]
                          json_input_file json_output_dir iterations replicates

positional arguments:
  json_input_file       File containing input parameters in json format
  json_output_dir       directory for output.File name will be input base name + rep number
  iterations            Maximum Iteration. Will run simulation till max iter if optional parameters are not provided
  replicates            Number of replicates each file will be saved with replicate number at the end

optional arguments:
  -h, --help            show this help message and exit
  -T MAX_TIME, --max_time MAX_TIME
                        Time limit at which the simuation is cutoff
  -I ITEM, --item ITEM  Name of reactant of product when specifying nummber at which cutoff happens. Should specify
                        --Item_level with it
  -L ITEM_LEVEL, --item_level ITEM_LEVEL
                        level of reactant of product at which cutoff happens. Should specify --item with it

In this example the commands used was: python param_gillespie.py ./param_jsons/example_parameter.json ./output_jsons/ 1000000 2 -T 10000
This runs the simulation with the starting paramter mentioned in ./param_jsons/example_parameter.json for 10000 seconds or 1000000 iterations, whichever occurs first.
This should output two files crrosponding to each replicate: example_parameter.json_rep_0_out.json and example_parameter.json_rep_1_out.json

3) Follow instructions Example_run.ipynb to read the output files and plotting graph.