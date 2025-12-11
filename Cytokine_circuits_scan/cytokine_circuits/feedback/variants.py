import pandas as pd

def generate_feedback_variants(base_df):
    """
    Generate feedback variants for the circuits in base_df.
    Paper 1.5:
    "Starting from each pareto-optimal... generated two 'neighbor' circuits by introducing an auto-regulatory loop on y."
    "Added auto-regulation was applied to the same term (production or removal) that mediates the regulation of y by x."
        """
    variants = []
    df_a = base_df[['A' in base_df.iloc[i]['model.num'] for i in range(len(base_df))]].copy().reset_index(drop = True)
    df_b = base_df[['B' in base_df.iloc[i]['model.num'] for i in range(len(base_df))]].copy().reset_index(drop = True)

    for i in range(len(df_a)):
        row_p = df_a.iloc[i].copy()
        row_m = df_a.iloc[i].copy()
        if df_a.iloc[i]['x_on_y_prod']!=0:
            row_p['y_on_y_prod'] = 1
            row_m['y_on_y_prod'] = -1
        elif df_a.iloc[i]['x_on_y_rem']!=0:
            row_p['y_on_y_rem'] = -1
            row_m['y_on_y_rem'] = 1
        row_p['model.num'] = row_p['model.num']+'+'
        df_a.loc[len(df_a)]= row_p
        row_m['model.num'] = row_m['model.num']+'-'
        df_a.loc[len(df_a)]= row_m
    
    for i in range(len(df_b)):
        row_p = df_b.iloc[i].copy()
        row_m = df_b.iloc[i].copy()
        if df_b.iloc[i]['x_on_y_prod']!=0:
            row_p['y_on_y_prod'] = 1
            row_m['y_on_y_prod'] = -1
        elif df_b.iloc[i]['x_on_y_rem']!=0:
            row_p['y_on_y_rem'] = -1
            row_m['y_on_y_rem'] = 1
        row_p['model.num'] = row_p['model.num']+'+'
        df_b.loc[len(df_b)]= row_p
        row_m['model.num'] = row_m['model.num']+'-'
        df_b.loc[len(df_b)]= row_m
    variants = pd.concat([df_a,df_b],ignore_index=True)

    return pd.DataFrame(variants)
