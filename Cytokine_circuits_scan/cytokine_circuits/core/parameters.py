import sympy as sp
def get_parameters(row):
    """
    Assign parameter values to a circuit based on its topological definition.
    
    Parameters:
        row (pd.Series): A row from the circuit DataFrame representing one circuit.
        
    Returns:
        dict: A dictionary of parameter values keyed by their symbolic names (strings or SymPy symbols).
              e.g. {a_x: 1.0, b_x: 5.0, ...}
    """
    # Placeholder values as requested.
    # User will insert their computation logic here.
    a_x,b_x,c_x,a_y,b_y,c_y,k_2x,k_3x,k_4x,k_2y,k_3y,k_4y = sp.symbols('a_x b_x c_x a_y b_y c_y k_2x k_3x k_4x k_2y k_3y k_4y')
    params = {
        a_x: 0.1, 
        b_x: 10.0, 
        c_x: 1.0, 
        a_y: 0.1, 
        b_y: 10.0, 
        c_y: 1.0,
        
        # Hill constants
        k_2x: 1.0, k_3x: 1.0, k_4x: 1.0,
        k_2y: 1.0, k_3y: 1.0, k_4y: 1.0
    }
    if 'A' in row['model.num']:
        if row['y_on_x_prod']==-1:
            params[a_x] = 0.01
            params[b_x] = 1.75
            params[c_x] = 0.75
            params[k_2x] = 0.5
            params[k_3x] = 0.
            params[k_4x] = 0.

        if row['y_on_x_prod']==-2:
            params[a_x] = 0.01
            params[b_x] = 1.75
            params[c_x] = 0.75
            params[k_2x] = 0.85
            params[k_3x] = 0.
            params[k_4x] = 0.

        if row['y_on_x_rem']==1:
            params[a_x] = 0.01
            params[b_x] = 1.38
            params[c_x] = 0.75
            params[k_4x] = 0.05
            params[k_3x] = 0.
            params[k_2x] = 0.

        if row['y_on_x_rem']==2:
            params[a_x] = 0.01
            params[b_x] = 1.472
            params[c_x] = 0.75
            params[k_4x] = 0.05
            params[k_3x] = 0.
            params[k_2x] = 0.

        if row['x_on_y_prod']==1:
            if '+' in row['model.num']:
                params[a_y] = 0.01
                params[b_y] = 0.5
                params[c_y] = 0.5/10
                params[k_2y] = 0.5
                params[k_3y] = 0.
                params[k_4y] = 0.
            elif '-' in row['model.num']:
                params[a_y] = 0.01
                params[b_y] = 0.12
                params[c_y] = 0.5/10
                params[k_2y] = 0.5
                params[k_3y] = 0.
                params[k_4y] = 0.
            else:
                params[a_y] = 0.01
                params[b_y] = 1./10
                params[c_y] = 0.5/10
                params[k_2y] = 0.5
                params[k_3y] = 0.
                params[k_4y] = 0.

        if row['x_on_y_prod']==2:
            if '+' in row['model.num']:
                params[a_y] = 0.01
                params[b_y] = 17.3
                params[c_y] = 0.5/10
                params[k_2y] = 0.5
                params[k_3y] = 0.
                params[k_4y] = 0.0
            elif '-' in row['model.num']:
                params[a_y] = 0.01
                params[b_y] = 4.4
                params[c_y] = 0.5/10
                params[k_2y] = 0.5
                params[k_3y] = 0.
                params[k_4y] = 0.0
            else:
                params[a_y] = 0.01
                params[b_y] = 1./10
                params[c_y] = 0.5/10
                params[k_2y] = 0.5
                params[k_3y] = 0.
                params[k_4y] = 0.

        if row['x_on_y_rem']==-1:
            if '-' in row['model.num']:
                params[a_y] = 1.7/10
                params[c_y] = 0.5/10
                params[k_3y] = 0.07
                params[k_2y] = 0.
                params[k_4y] = 0.1
                params[b_y] = 0
            elif '+' in row['model.num']:
                params[a_y] = 1.7/10
                params[c_y] = 0.5/10
                params[k_3y] = 0.0175
                params[k_4y] = 0.1
                params[b_y] = 0.
            else:
                params[a_y] = 1.7/10
                params[c_y] = 0.5/10
                params[k_3y] = 0.06
                params[k_2y] = 0.
                params[k_4y] = 0.
                params[b_y] = 0

        if row['x_on_y_rem']==-2:
            if '+' in row['model.num']:
                params[a_y] = 0.17
                params[c_y] = 0.0009
                params[k_3y] = 0.06
                params[k_2y] = 0.
                params[k_4y] = 0.1
                params[b_y] = 0.
            elif '-' in row['model.num']:
                params[a_y] = 0.17
                params[c_y] = 0.0035
                params[k_3y] = 0.06
                params[k_2y] = 0.
                params[k_4y] = 0.1
                params[b_y] = 0.

            else:
                params[a_y] = 0.17
                params[c_y] = 0.0025
                params[k_3y] = 0.06
                params[k_2y] = 0.
                params[k_4y] = 0.
                params[b_y] = 0.

    if 'B' in row['model.num']:
        if row['y_on_x_prod']==1:
            params[a_x] = 0.05
            params[b_x] = 5.0
            params[c_x] = 0.75
            params[k_2x] = 0.5
            params[k_3x] = -1.
            params[k_4x] = -1.
        if row['y_on_x_prod']==2:
            params[a_x] = 0.05
            params[b_x] = 3.26
            params[c_x] = 0.75
            params[k_2x] = 0.25
            params[k_3x] = -1.
            params[k_4x] = -1.
        if row['y_on_x_rem']==-1:
            params[a_x] = 0.05
            params[b_x] = 2.5
            params[c_x] = 0.75
            params[k_2x] = -1.
            params[k_3x] = -1.
            params[k_4x] = 0.4
        if row['y_on_x_rem']==-2:
            params[a_x] = 0.06
            params[b_x] = 2.1
            params[c_x] = 0.75
            params[k_2x] = -1.
            params[k_3x] = -1.
            params[k_4x] = 0.8


        if row['x_on_y_prod']==-1:
            if '+' in row['model.num']:
                params[a_y] = 0.05/10
                params[b_y] = 0.0174
                params[c_y] = 1./10
                params[k_2y] = 0.05
                params[k_3y] = -1.
                params[k_4y] = -1.
            elif '-' in row['model.num']:
                params[a_y] = 0.05/10
                params[b_y] = 0.007
                params[c_y] = 1./10
                params[k_2y] = 0.05
                params[k_3y] = -1.
                params[k_4y] = -1.
            else:
                params[a_y] = 0.05/10
                params[b_y] = 0.05/10
                params[c_y] = 1./10
                params[k_2y] = 0.05
                params[k_3y] = -1.
                params[k_4y] = -1.
                
        if row['x_on_y_prod']==-2:
            if '+' in row['model.num']:
                params[a_y] = 0.05/10
                params[b_y] = 0.0013
                params[c_y] = 1./10
                params[k_2y] = 0.05
                params[k_3y] = -1.
                params[k_4y] = -1.
            elif '-' in row['model.num']:
                params[a_y] = 0.05/10
                params[b_y] = 0.005/10
                params[c_y] = 1./10
                params[k_2y] = 0.05
                params[k_3y] = -1.
                params[k_4y] = -1.
            else:
                params[a_y] = 0.05/10
                params[b_y] = 0.004/10
                params[c_y] = 1./10
                params[k_2y] = 0.05
                params[k_3y] = -1.
                params[k_4y] = -1.


        if row['x_on_y_rem']==1:
            if '+' in row['model.num']:
                params[a_y] = 0.05/10
                params[b_y] = 0.05/10
                params[c_y] = 0.15
                params[k_2y] = -1.
                params[k_3y] = 0.3
                params[k_4y] = 1.0
            elif '-' in row['model.num']:
                params[a_y] = 0.05/10
                params[b_y] = 0.05/10
                params[c_y] = 0.055
                params[k_2y] = -1.
                params[k_3y] = 0.3
                params[k_4y] = 0.1
            else:
                params[a_y] = 0.05/10
                params[b_y] = 0.05/10
                params[c_y] = 1./10
                params[k_2y] = -1.
                params[k_3y] = 0.3
                params[k_4y] = -1.

        if row['x_on_y_rem']==2:
            if '+' in row['model.num']:
                params[a_y] = 0.0
                params[b_y] = 0.05/10
                params[c_y] = 0.15
                params[k_2y] = -1.
                params[k_3y] = 0.3
                params[k_4y] = 0.1
            elif '-' in row['model.num']:
                params[a_y] = 0.0
                params[b_y] = 0.05/10
                params[c_y] = 0.36
                params[k_2y] = -1.
                params[k_3y] = 0.3
                params[k_4y] = 0.1
            else:
                params[a_y] = 0.05/10
                params[b_y] = 0.05/10
                params[c_y] = 0.3
                params[k_2y] = -1.
                params[k_3y] = 0.3
                params[k_4y] = -1.

    return params
