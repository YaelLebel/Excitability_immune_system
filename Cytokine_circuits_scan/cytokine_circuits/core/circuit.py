import sympy as sp

class Circuit:
    """
    Represents a Cytokine-Cytokine circuit with 3 interactions.
    Refactored from Combo_pp.
    """
    _id_counter = 1

    def __init__(self, xdot, ydot, parameters, num_interactions,
                 x=sp.symbols('x'), y=sp.symbols('y'),
                 parameter_placement='def', find_nullclines=False):
        self.x = x
        self.y = y
        self.variables = [x, y]
        self.xdot = xdot
        self.ydot = ydot
        self.parameters = parameters
        self.num_dof = num_interactions
        self.parameter_placement = parameter_placement
        
        self.nullclines_x = []
        self.nullclines_y = []
        self.nullclines = []

        if find_nullclines:
            self.find_nullclines()
            
        self.id = Circuit._id_counter
        Circuit._id_counter += 1

    @classmethod
    def reset_id_counter(cls):
        cls._id_counter = 1

    def find_nullclines(self):
        # Result as a function of y (x = f(y))
        # Note: SymPy solve can return multiple solutions
        # We need to handle them carefully.
        
        # Original code logic:
        # nc_x_x = sp.solve(sp.Eq(self.xdot,0),x) #result as a function of y
        nc_x_y = sp.solve(sp.Eq(self.xdot, 0), self.y) #result as a function of x
        
        # nc_y_x = sp.solve(sp.Eq(self.ydot,0),x)
        nc_y_y = sp.solve(sp.Eq(self.ydot, 0), self.y) #result as a function of x (y = f(x))

        self.nullclines_x = []
        # Storing as (x, solution_expression) pairs
        for n in nc_x_y:
            self.nullclines_x.append((self.x, n))

        self.nullclines_y = []
        for n in nc_y_y:
            self.nullclines_y.append((self.x, n))
            
        self.nullclines = self.nullclines_x + self.nullclines_y
