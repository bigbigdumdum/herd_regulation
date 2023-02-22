import numpy as np
# Version 0.5 (compactable with 0.3 scripts)
# Change log form 0.2
# +-----------------+
# 1) Propensity is calcualted as rate now instead of rate/sum rate (solves the problem with time)
# 2) update_reactant_log() function rewritten as reactions of type A -> A + B cause writing log twice
# 3) new fundtion: return_current_time() that returns the last time of the reaction, this allows reactions to
#    to proceed to the same time but differnt iterations.
# 4) Bug in calculating time: self.tau = (1/self.sum_prop)*np.log(1/(1-rand1)) is correct instead of older self.tau =       #    (1/self.sum_prop)*(1/(1-rand1)). This was the source of super long time steps. Thanks to Aditya Marodia for   
#    identifying the bug.(change from 0.4)
# 5) Added a max iteration cutoff for run till time and run_till level, default is 1 bil.(change from 0.4)
#
#
# By Mukundan S 
class Gillespie:
    ## USAGE EXAMPLE #
    #================#
#     # Set up inputs
#     n_iter = 2500
#     Kd = 2
#     kf,kr = 1,Kd
#     reactant_state = {'A':500,'B':500,'AB':0} # starting states in number of molecules 
#     rate_constant_list = [kf,kr]
#     reactants_list = [['A','B'],['AB']] ## for 2A it must be ['A','A'] (Here each index is a reaction with products under the same index in product_list)
#     product_list = [['AB'],['A','B']]
#     # Initiate object
#     reaction = Gillespie(rate_constant_list,reactants_list,product_list,reactant_state)
#     time,reactants_log = reaction.autorun(n_iter)
    #================#
    def __init__(self,rate_constant_list,reactants_list,product_list,reactant_state):# n is the total number of the reaction
        self.time_list = [0]     
        self.reaction_log = [] # list of indices of reactions
        self.rate_constant_list = rate_constant_list # INPUT list of rate constants
        self.reactants_list = reactants_list # INPUT nested list of reactants in string fromat eg [['A','B'],['AB'],['A','Bi'],['ABi']] 
        self.product_list = product_list # INPUT nested list of crrosponding products. eg [['AB'],['A','B'],['ABi'],['A','Bi']]
        self.reactant_state = reactant_state # INPUT a dictionary specifying number of molecules of each reactant and product eg {'A':500,'B':500,'AB':0,'Bi':Bi_n,'ABi':0}
        self.reactant_log = {}
        self.prop_log = []
        for reactants,products in zip(reactants_list,product_list):
            for x in reactants:
                self.reactant_log[x] = [reactant_state[x]]
            for x in products:
                self.reactant_log[x] = [reactant_state[x]]
        self.all_compounds = list(self.reactant_log.keys())
        self.net_change_dict = dict([(x,0) for x in self.all_compounds])
        return
        
    def load_propansities(self,prop_list):
        self.prop_list = prop_list
        self.sum_prop = sum(prop_list)
        return
#     def calc_propencity(self):
#         rate_list = []
#         for i,reactants in enumerate(self.reactants_list):
#             reactant_product = 1
#             for reactant in reactants:
#                 reactant_product = reactant_product*self.reactant_log[reactant][-1]
#             rate = self.rate_constant_list[i]*reactant_product
#             rate_list.append(rate)
#         self.prop_list = [rate/sum(rate_list) for rate in rate_list]
#         self.prop_log.append(self.prop_list)
#         self.sum_prop = sum(self.prop_list)
    
    def calc_propencity(self): # calculates propensity of the reaction based on Gillispie 2013 using k as cj
        rate_list = []
        for i,reactants in enumerate(self.reactants_list):
            reactant_product = 1
            if len(set(reactants)) == 1 and len(reactants) == 2: #special case for association reation eg 2A -> B
                rate = self.rate_constant_list[i]*self.reactant_log[reactants[0]]*(self.reactant_log[reactants[0]]-1)
                rate_list.append(rate)
            else:
                for reactant in reactants:
                    reactant_product = reactant_product*self.reactant_log[reactant][-1]
                rate = self.rate_constant_list[i]*reactant_product
                rate_list.append(rate)
        if len(rate_list) != len(self.reactants_list):
            print('Some propensity not calculated!, Bad simulation')# check point for propensity
#         self.prop_list = [rate/sum(rate_list) for rate in rate_list] # propencity is calulate as ratio of rate to total rate
#         self.prop_log.append(self.prop_list)
#         self.sum_prop = sum(self.prop_list)
        self.prop_list = rate_list # propencity is calulate as rate
        self.prop_log.append(self.prop_list)
        self.sum_prop = sum(self.prop_list)
        return
        
    def calc_tau_and_update_time(self,rand1): # Caclulates Tau from random numbers see gillespie 2013 for more.        
        # self.tau = (1/self.sum_prop)*(1/(1-rand1))        
        self.tau = (1/self.sum_prop)*np.log(1/(1-rand1))    
        self.time_list.append(self.time_list[-1]+self.tau)
        return
        
    def pick_reaction(self,rand2): # Monte Carlo sampling to pick reaction based on propensities and a random number
        sigma_a = 0
        for i,prop in enumerate(self.prop_list):
            sigma_a+=prop
            if sigma_a > rand2*self.sum_prop: # < in gillespie 2013 but figured it was a typo
                self.reaction_log.append(i)
                self.lucky_winner = i
                break
        return
                
#     def update_reactant_log(self): # based on the reaction picked the reactant and product counts are updated
#         lucky_list = [] # list to identify reactants and products that are updated in this iteration
#         for reactant in self.reactants_list[self.lucky_winner]:
#             self.reactant_log[reactant].append(self.reactant_log[reactant][-1]-1) # reactant log is a dictionary that maps reactants to a list(list[-1] = last resction state)
#             lucky_list.append(reactant)
#         for product in self.product_list[self.lucky_winner]:
#             self.reactant_log[product].append(self.reactant_log[product][-1]+1)
#             lucky_list.append(product)
#         # updating unreacted list
#         lucky_list = list(set(lucky_list))  
#         unlucky_list = list(self.reactant_log.keys())
#         for compound in list(set(lucky_list)): # Identifying unreacted compounds by removing compounds part of reaction from U
#             unlucky_list.remove(compound)
#         for unlucky in unlucky_list: #  new value is the same value as the previous iteration as this reaction did not take place.
#             self.reactant_log[unlucky].append(self.reactant_log[unlucky][-1])
            
    def update_reactant_log(self):
        # minus one for reactants
        for reactant in self.reactants_list[self.lucky_winner]:
            self.net_change_dict[reactant] -= 1
        # plus one for reactants
        for product in self.product_list[self.lucky_winner]:
            self.net_change_dict[product] += 1
#         print(self.lucky_winner,self.net_change_dict)
        # updating reactant log
        for compound in self.all_compounds:
            self.reactant_log[compound].append(self.reactant_log[compound][-1]+self.net_change_dict[compound])
            self.net_change_dict[compound] = 0
        return
        
    def return_reaction_progression(self): # returns relevant info back to the user
        return self.time_list,self.reactant_log
    
    def return_current_time(self):
        return self.time_list[-1]
    
    def autorun(self,n_iter): # runs iterations n_iter times automatically
        import random
        for i in range(n_iter):
            # generate random numbers
            rand1 = random.random() 
            rand2 = random.random()
            self.calc_propencity()
            self.calc_tau_and_update_time(rand1)
            self.pick_reaction(rand2)
            self.update_reactant_log()
        return self.return_reaction_progression()
    
    def run_till_time(self,time,max_iter=1000000000):
        import random
        for i in range(max_iter):
        # while self.return_current_time() <= time:
            rand1 = random.random() 
            rand2 = random.random()
            self.calc_propencity()
            self.calc_tau_and_update_time(rand1)
            self.pick_reaction(rand2)
            self.update_reactant_log()
            if self.return_current_time() >= time:
                break
        return self.return_reaction_progression()

    def run_till_product_level(self,product,level,max_iter=1000000000):
        import random
        for i in range(max_iter):
        # while self.reactant_log[product][-1] < level:
            rand1 = random.random() 
            rand2 = random.random()
            self.calc_propencity()
            self.calc_tau_and_update_time(rand1)
            self.pick_reaction(rand2)
            self.update_reactant_log()
            if self.reactant_log[product][-1] == level:
                break
        return self.return_reaction_progression()
    
    def add_reactants(self,reactants,increment):
        for i,reactant in enumerate(reactants):
            self.reactant_log[reactant][-1] += increment[i]
        return    
   
    def debug(self): # prints important lists for debugging purpose
#         print('Time list :',self.time_list)
        print('\nReaction list :',self.reaction_log)
        print('\nReactant progression :\n')
        for reactant in self.reactant_log:
            print(reactant, self.reactant_log[reactant],'\n')
        for prop in self.prop_log:
            print(prop,'\n')            
        return