# -*- coding: utf-8 -*-
"""
Created on Wed May 24 10:56:19 2017

@author: Javier Castillo
"""

number_of_colonies = 53.6
plasmid_volume = 2 # expressed in ul
plasmid_concentration = 0.085 # expressed in ug/ul
dilution_factor = 100 
reaction_volume = 1000 # expressed in ul
plated_volume = 100 # expressed in ul
added_plasmid = (plasmid_volume * plasmid_concentration * (plated_volume / reaction_volume)) / dilution_factor 

transformation_efficiency = number_of_colonies / added_plasmid
print("Transformation Efficiency:", transformation_efficiency, "CFU/ugDNA")