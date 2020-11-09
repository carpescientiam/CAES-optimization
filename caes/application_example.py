# -*- coding: utf-8 -*-
"""
Created on Mon Nov  2 17:25:30 2020

@author: Fahim Sadat
"""
from caes.model import*



caes = Diabatic(V_cas=310000, P_cmp=60, P_exp=60)
caes.optimize('prices_example', 'cost_example', 10)

