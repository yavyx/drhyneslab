# -*- coding: utf-8 -*-
"""
Created on Wed May 24 13:47:35 2017

@author: Javier
"""

molecular_weight = 58.44 # expressed in g/mol
molarity = 0.6
final_volume = 1000 # expressed in ml
conversion_factor = final_volume / 1000 # expressed in ml
needed_moles = molarity * conversion_factor

amount_needed = molecular_weight * needed_moles # expressed in grams

print("Weigh", amount_needed, "grams.")