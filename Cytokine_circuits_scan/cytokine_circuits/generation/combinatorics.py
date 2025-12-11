import pandas as pd
from itertools import product, combinations
import numpy as np

def generate_pp_combinations(x_count, y_count,
                             default_value = '0',
                             x_columns = ['x_on_x_prod', 'x_on_x_rem','y_on_x_prod', 'y_on_x_rem',],
                             y_columns = ['x_on_y_prod', 'x_on_y_rem','y_on_y_prod', 'y_on_y_rem',],
                             interaction_values = ['+', '-']):
    """
    Generate a DataFrame containing all combinations of interactions.

    Parameters:
        x_count (int): Number of nonzero interactions for x.
        y_count (int): Number of nonzero interactions for y.
    
    Returns:
        pd.DataFrame: DataFrame with all interaction combinations.
    """
    total_nonzero = x_count + y_count 

    # All columns
    all_columns = x_columns + y_columns

    # Generate all combinations of `total_nonzero` non-zero interactions
    nonzero_combinations = list(product(interaction_values, repeat=total_nonzero))

    # Generate all positions for non-zero values in the columns
    nonzero_positions = list(combinations(range(len(all_columns)), total_nonzero))

    # Prepare the final list of rows
    rows = []

    for positions in nonzero_positions:
        for values in nonzero_combinations:
            # Initialize a row with all zeros
            row = [default_value] * len(all_columns)

            # Assign non-zero values to the specified positions
            for pos, value in zip(positions, values):
                row[pos] = value

            # Check if the division of non-zero terms satisfies the x_count and y_count
            x_nonzero = sum(1 for i in range(len(x_columns)) if row[i] != default_value)
            y_nonzero = sum(1 for i in range(len(x_columns), len(all_columns)) if row[i] != default_value)

            if x_nonzero == x_count and y_nonzero == y_count:
                rows.append(row)

    # Create a DataFrame
    df = pd.DataFrame(rows, columns=all_columns)

    return df


