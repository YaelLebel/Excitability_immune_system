import sympy as sp
from .circuit import Circuit

def row_to_sp (row):
    x,y,a_x,b_x,c_x,a_y,b_y,c_y,k_2x, k_3x, k_4x, k_2y, k_3y, k_4y = sp.symbols('x y a_x b_x c_x a_y b_y c_y k_2x k_3x k_4x k_2y k_3y k_4y')
    if row['x_on_x_prod']>=0:
        n_x = k_x = row['x_on_x_prod']
    else:
        n_x = 0
        k_x = abs(row['x_on_x_prod'])
    if row['y_on_x_prod']>=0:
        m_x = l_x = row['y_on_x_prod']
    else:
        m_x = 0
        l_x = abs(row['y_on_x_prod'])
    if row['x_on_x_rem']>=0:
        p_x = r_x = row['x_on_x_rem']
    else:
        p_x = 0
        r_x = abs(row['x_on_x_rem'])
    if row['y_on_x_rem']>=0:
        q_x = s_x = row['y_on_x_rem']
    else:
        q_x = 0
        s_x = abs(row['y_on_x_rem'])

    xdot = a_x
    
    if l_x != 0:
        xdot = xdot + b_x *(x**n_x/(1+x**k_x))*(y**m_x/(k_2x**l_x+y**l_x))
    else:
        xdot = xdot + b_x *(x**n_x/(1+x**k_x))

    if r_x==0:
        if s_x == 0 :
            xdot = xdot - c_x * x
        else:
            xdot = xdot - c_x * x * (y**q_x/(k_4x**s_x+y**s_x))
    else:
        if s_x == 0:
            xdot = xdot - c_x * x * (x**p_x/(k_3x**r_x+x**r_x)) 
        else:
            xdot = xdot - c_x * x * (x**p_x/(k_3x**r_x+x**r_x)) * (y**q_x/(k_4x**s_x+y**s_x))

    if row['x_on_y_prod']>=0:
        n_y = k_y = row['x_on_y_prod']
    else:
        n_y = 0
        k_y = abs(row['x_on_y_prod'])
    if row['y_on_y_prod']>=0:
        m_y = l_y = row['y_on_y_prod']
    else:
        m_y = 0
        l_y = abs(row['y_on_y_prod'])
    if row['x_on_y_rem']>=0:
        p_y = r_y = row['x_on_y_rem']
    else:
        p_y = 0
        r_y = abs(row['x_on_y_rem'])
    if row['y_on_y_rem']>=0:
        q_y = s_y = row['y_on_y_rem']
    else:
        q_y = 0
        s_y = abs(row['y_on_y_rem'])

    if row['x_on_y_prod']!=0 and row['y_on_y_prod']!=0:
        ydot = a_y + b_y *(x**n_y/(k_2y**k_y+x**k_y))*(y**m_y/(1+y**l_y))
    elif row['x_on_y_prod']==0 and row['y_on_y_prod']!=0:
        ydot = a_y + b_y *(y**m_y/(1+y**l_y))
    elif row['x_on_y_prod']!=0 and row['y_on_y_prod']==0:
        ydot = a_y + b_y *(x**n_y/(k_2y**k_y+x**k_y))
    else:
        ydot = a_y

    if r_y==0:
        if s_y == 0 :
            ydot = ydot - c_y * y
        else:
            ydot = ydot - c_y * y * (y**q_y/(k_4y**s_y+y**s_y))
    else:
        if s_y == 0:
            ydot = ydot - c_y * y * (x**p_y/(k_3y**r_y+x**r_y)) 
        else:
            ydot = ydot - c_y * y * (x**p_y/(k_3y**r_y+x**r_y)) * (y**q_y/(k_4y**s_y+y**s_y))

    return xdot, ydot, [x,y],[a_x,b_x,c_x,a_y,b_y,c_y,k_2x, k_3x, k_4x, k_2y, k_3y, k_4y]

def row_to_circuit(row):
    """
    Wrapper to convert a row directly to a Circuit object.
    """
    xdot, ydot, vars, params = row_to_sp(row)
    
    # Calculate number of interactions
    # 8 interaction columns
    cols = ['x_on_x_prod', 'x_on_x_rem', 'y_on_x_prod', 'y_on_x_rem',
            'x_on_y_prod', 'x_on_y_rem', 'y_on_y_prod', 'y_on_y_rem']
    
    num_interactions = sum(1 for c in cols if str(row.get(c, 0)) != '0')
    
    return Circuit(xdot, ydot, params, num_interactions, x=vars[0], y=vars[1])
