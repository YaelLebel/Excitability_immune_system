import pandas as pd
import sympy as sp
import numpy as np
def filter_feedback (df):
    """
    Remove circuits that have feedback.
    A circuit is connected if there is at least one interaction in Px or Rx involving y,
    AND at least one interaction in Py or Ry involving x.
    
    In the DataFrame columns:
    x_columns = ['x_on_x_prod', 'x_on_x_rem', 'y_on_x_prod', 'y_on_x_rem']
    y_columns = ['x_on_y_prod', 'x_on_y_rem', 'y_on_y_prod', 'y_on_y_rem']
    
    Connectedness condition (Paper 1.2):
    "no interaction in either Px or Rx involving y" -> y_on_x_prod == 0 AND y_on_x_rem == 0
    "no interaction in Py or Ry involving x" -> x_on_y_prod == 0 AND x_on_y_rem == 0
    """

    mask_connected = []
    for i in range(len(df)):
        if (df.iloc[i]['y_on_x_prod']==0 and df.iloc[i]['y_on_x_rem']==0) or (df.iloc[i]['x_on_y_prod']==0 and df.iloc[i]['x_on_y_rem']==0):
            mask_connected.append(False)
        else:
            mask_connected.append(True)
    df_ret = df[mask_connected].copy().reset_index(drop = True)
    return df_ret

def filter_symmetric(df):
    """
    Remove circuits that are x-y symmetric to other circuits.
    Each circuit can be represented by its interaction vector.
    We swap x and y definitions and check if the resulting circuit is already in the set.
    
    Mapping for swap (x <-> y):
    x_on_x_* <-> y_on_y_*
    y_on_x_* <-> x_on_y_*
    """
    # Define column mapping for symmetry swap
    swap_map = {
        'x_on_x_prod': 'y_on_y_prod',
        'x_on_x_rem': 'y_on_y_rem',
        'y_on_x_prod': 'x_on_y_prod',
        'y_on_x_rem': 'x_on_y_rem',
        
        'y_on_y_prod': 'x_on_x_prod',
        'y_on_y_rem': 'x_on_x_rem',
        'x_on_y_prod': 'y_on_x_prod',
        'x_on_y_rem': 'y_on_x_rem'
    }
    
    seen_hashes = set()
    indices_to_keep = []
    
    for idx, row in df.iterrows():
        # Create a tuple signature for the current circuit
        # Ensure we use a consistent order of columns for the signature
        cols = sorted(swap_map.keys()) # sort to ensure deterministic order
        
        current_sig = tuple(str(row[c]) for c in cols)
        
        # Create signature for the swapped circuit
        swapped_sig = tuple(str(row[swap_map[c]]) for c in cols)
        
        # If we haven't seen this circuit (or its symmetric partner), keep it
        if current_sig not in seen_hashes:
            indices_to_keep.append(idx)
            seen_hashes.add(current_sig)
            seen_hashes.add(swapped_sig)
            
    return df.loc[indices_to_keep].copy()


def filter_diverging(df):
    """
    "Remove circuits that are not contained and diverge to infinity, namely, circuits in which one variable has only positive interactions and no negative interactions with itself or the other variable."
    
    Positive interactions (promote growth):
    - Production activation (val > 0 on prod)
    - Removal inhibition (val < 0 on rem)
    
    Negative interactions (promote stability):
    - Production inhibition (val < 0 on prod)
    - Removal activation (val > 0 on rem)
    
    Condition for divergence (for a variable V):
    Has(Positive) AND NOT(Has_Negative)
    
    Note: Linear degradation (-rx*x) is intrinsic, but paper implies we filter based on regulatory logic.
    "no negative interactions with itself or the other variable"
    """
    mask_contained = []
    for i in range(len(df)):
        prod_list = df.iloc[i][['x_on_x_prod','y_on_x_prod','x_on_y_prod','y_on_y_prod']].values
        prod_set = set(df.iloc[i][['x_on_x_prod','y_on_x_prod','x_on_y_prod','y_on_y_prod']])
        prod_set_no_zeros = set(prod_list[prod_list!=0])
        pos_subset = {1,2}
        rem_list = df.iloc[i][['x_on_x_rem','y_on_x_rem','x_on_y_rem','y_on_y_rem']].values
        rem_set = set(rem_list)
        rem_set_no_zeros = set(rem_list[rem_list!=0])
        neg_subset = {-1,-2}

        conds = [len(prod_list[prod_list!=0])==3 and prod_set_no_zeros.issubset(pos_subset), 
                len(rem_list[rem_list!=0])==3 and rem_set_no_zeros.issubset(neg_subset),
                prod_set_no_zeros.issubset(pos_subset) and rem_set_no_zeros.issubset(neg_subset) ]
        if np.any(conds):
            mask_contained.append(False)
        else:
            mask_contained.append(True)
    # We keep circuits that are NOT diverging
    return df[mask_contained].copy().reset_index(drop = True)

def filter_N_shaped (df):
    mask_threshold = []
    for i in range(len(df)):
        if df.iloc[i]['x_on_x_prod']!=2:
            mask_threshold.append(False)
        else:
            mask_threshold.append(True)
    return df[mask_threshold].copy().reset_index(drop = True)

def filter_fixed_points (df):
    mask_only_one_fixed_point = []
    for i in range(len(df)):
        list_of_conds = [df.iloc[i]['x_on_x_prod']==2 and df.iloc[i]['y_on_x_prod']==-1 and df.iloc[i]['x_on_y_prod']==-1,
                        df.iloc[i]['x_on_x_prod']==2 and df.iloc[i]['y_on_x_prod']==-1 and df.iloc[i]['x_on_y_prod']==-2,
                        df.iloc[i]['x_on_x_prod']==2 and df.iloc[i]['y_on_x_prod']==-2 and df.iloc[i]['x_on_y_prod']==-2,
                        df.iloc[i]['x_on_x_prod']==2 and df.iloc[i]['y_on_x_prod']==-2 and df.iloc[i]['x_on_y_prod']==-1,
                        df.iloc[i]['x_on_x_prod']==2 and df.iloc[i]['y_on_x_prod']==-1 and df.iloc[i]['x_on_y_rem']==1,
                        df.iloc[i]['x_on_x_prod']==2 and df.iloc[i]['y_on_x_prod']==-1 and df.iloc[i]['x_on_y_rem']==2,
                        df.iloc[i]['x_on_x_prod']==2 and df.iloc[i]['y_on_x_prod']==-2 and df.iloc[i]['x_on_y_rem']==1,
                        df.iloc[i]['x_on_x_prod']==2 and df.iloc[i]['y_on_x_prod']==-2 and df.iloc[i]['x_on_y_rem']==2,
                        df.iloc[i]['x_on_x_prod']==2 and df.iloc[i]['y_on_x_rem']==1 and df.iloc[i]['x_on_y_prod']==-1,
                        df.iloc[i]['x_on_x_prod']==2 and df.iloc[i]['y_on_x_rem']==1 and df.iloc[i]['x_on_y_prod']==-2,
                        df.iloc[i]['x_on_x_prod']==2 and df.iloc[i]['y_on_x_rem']==2 and df.iloc[i]['x_on_y_prod']==-1,
                        df.iloc[i]['x_on_x_prod']==2 and df.iloc[i]['y_on_x_rem']==2 and df.iloc[i]['x_on_y_prod']==-2,
                        df.iloc[i]['x_on_x_prod']==2 and df.iloc[i]['y_on_x_rem']==1 and df.iloc[i]['x_on_y_rem']==1,
                        df.iloc[i]['x_on_x_prod']==2 and df.iloc[i]['y_on_x_rem']==1 and df.iloc[i]['x_on_y_rem']==2,
                        df.iloc[i]['x_on_x_prod']==2 and df.iloc[i]['y_on_x_rem']==2 and df.iloc[i]['x_on_y_rem']==1,
                        df.iloc[i]['x_on_x_prod']==2 and df.iloc[i]['y_on_x_rem']==2 and df.iloc[i]['x_on_y_rem']==2,
                        ]
        if np.any(list_of_conds):
            mask_only_one_fixed_point.append(False)
        else:
            mask_only_one_fixed_point.append(True)
    return df[mask_only_one_fixed_point].copy().reset_index(drop = True)
