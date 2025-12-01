import numpy as np
import sympy as sp
import re
from Base.coeff_base import coeff
from Base.Combo_Base import Combo
from itertools import combinations

def check_is_duplicate_parameter(ps):
    rests = []
    for p in ps:
        _,rest = p.split_symbol()
        rests.append(rest)
    if len(rests) == len(set(rests)):
        return True
    return False

def create_combos (num_parameters:int):
    parameters_ = ['a00','a10','a20','a01','a02','a11','b00','b01','b02','b10','b20','b11']
    parameters = ['+'+p for p in parameters_]
    parameters += ['-'+p for p in parameters_]
    parameters = [coeff(p) for p in parameters]
    combos = combinations(parameters,num_parameters)
    combos_ret = []
    for comb in combos:
        if check_is_duplicate_parameter(comb):
            combo = Combo(comb)
            combos_ret.append(combo)
    return combos_ret


def create_combos_from_csv (filename):
    import pandas as pd
    df = pd.read_csv(filename)
    combos_ret = []
    for i in range(len(df)):
        ps = [coeff(df.iloc[i]['p1']),coeff(df.iloc[i]['p2']),coeff(df.iloc[i]['p3'])]
        ident = int(df.iloc[i]['id_counter'])
        combo_to_append = Combo(ps,ident)
        combos_ret.append(combo_to_append)
    return combos_ret
    
def create_combos_from_csv_with_conds(filename):
    import pandas as pd
    df = pd.read_csv(filename,index_col = 'Unnamed: 0')
    combos_ret = []
    print(df.to_string())
    for i in range(len(df)):
        ps = [coeff(df.iloc[i]['p1']),coeff(df.iloc[i]['p2']),coeff(df.iloc[i]['p3'])]
        ident = int(df.iloc[i]['id_counter'])
        combo_to_append = Combo(ps,ident)
        combo_to_append.is_two_variable_dynamics = df.iloc[i]['is_two_variable_dynamics']
        combo_to_append.is_connected = df.iloc[i]['is_connected']
        combo_to_append.is_contained_first_quadrant = df.iloc[i]['is_contained_first_quadrant']
        combo_to_append.is_zero_stable = df.iloc[i]['is_zero_stable']
        combo_to_append.is_nullcline_in_first_quadrant = df.iloc[i]['is_nullcline_in_first_quadrant']
        combo_to_append.is_no_other_symmetric = df.iloc[i]['is_no_other_symmetric']
        combos_ret.append(combo_to_append)

    return combos_ret,df,df.columns.to_list()[-6:]

def check_conditions (item,cond_list):
    #cond_list is a list of functions that returns 0 if condition is not fulfilled and anything else otherwise
    #the function return 1 if the item matches all conditions and 0 otherwise
    for cond in cond_list: 
        if cond(item)==0:
            return 0
    return 1

def check_conditions_list (list,cond_list):
    sum = 0
    for item in list:
        add = check_conditions(item,cond_list)
        sum+=add
    return sum

def filter_conditions_list (list,cond_list):
    ret = []
    for item in list:
        add = check_conditions(item,cond_list)
        if add==1:
            ret.append(item)
    return ret
