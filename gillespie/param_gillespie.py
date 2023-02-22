import json
from Gillespie import Gillespie

class Param_Gillespie:
    def __init__(self,json_file,out_file):
        # load objects from json
        with open(json_file) as inf:
            rate_constant_list,reactants_list,product_list,reactant_state = json.load(inf) 
        # create Gillespie object
        self.g = Gillespie(rate_constant_list,reactants_list,product_list,reactant_state)
        # save outfile name for future use
        self.out_file = out_file
        
    def run_till_steps(self,steps):
        with open(self.out_file,'w') as inf:
            json.dump(self.g.autorun(steps),inf)
            
    def run_till_time(self,time,max_iter=100000):
        with open(self.out_file,'w') as inf:
            json.dump(self.g.run_till_time(time,max_iter=max_iter ),inf)
        
    def run_till_product_level(self,product,level,max_iter=100000):        
        with open(self.out_file,'w') as inf:
            json.dump(self.g.run_till_product_level(product,level,max_iter=max_iter),inf)
            
            

import json
import sys
import os
import argparse

# parsing arguments
parser = argparse.ArgumentParser()
# positional arguments
# 1) json_input_file
parser.add_argument("json_input_file", help="File containing input parameters in json format")
# 2) json_output_file
parser.add_argument("json_output_dir", help="directory for output.File name will be input base name + rep number")
# 3) Iterations 
parser.add_argument("iterations", help="Maximum Iteration. Will run simulation till max iter if optional parameters are not provided")
parser.add_argument("replicates", help="Number of replicates each file will be saved with replicate number at the end")
# optional argument
# 1) max_time
parser.add_argument("-T","--max_time", help="Time limit at which the simuation is cutoff")
# 2) item
parser.add_argument("-I", "--item", help="Name of reactant of product when specifying nummber at which cutoff happens. Should specify --Item_level with it")
# 3) item_level
parser.add_argument("-L", "--item_level", help="level of reactant of product at which cutoff happens. Should specify --item with it")

# parse argumenst  
args = parser.parse_args()

for i in range(int(args.replicates)): # replicates    

    # print(args.json_input_file,args.json_output_file,args.iterations,args.max_time,args.item,args.item_level)
    base_name = os.path.basename(args.json_input_file)
    out_file = args.json_output_dir+'/'+base_name+'_rep_'+str(i)+'_out.json'
    pg = Param_Gillespie(args.json_input_file,out_file)

    max_iter = int(args.iterations)

    # exceptions
    if bool(args.item) ^ bool(args.item_level): # xor to check if both are present or both are absent
        print('Provide both --item_level and --item')
        sys.exit()

    if args.max_time and args.item and args.item_level:
        print('Provide either --item_level and --item, or --max_time')
        sys.exit()

    # run till steps
    if not args.max_time and not args.item:
        pg.run_till_steps(max_iter)

    # run till time
    if args.max_time:
        time = float(args.max_time)    
        pg.run_till_time(time,max_iter=max_iter)

    if args.item and args.item_level:
        item = str(args.item)
        level = int(args.item_level)
        pg.run_till_product_level(item,level,max_iter=max_iter)
