import numpy as np

def calculate_response_metrics(t, x, y):
    """
    Calculate Response Strength and Speed.
    
    Paper 1.4:
    "Response strength was defined as the maximal x output."
    "Response speed was defined as the maximal x output divided by the time it took to get to this maximum."
    """
    # Assuming x is the cytokine of interest
    max_x = np.max(x)
    max_idx = np.argmax(x)
    time_to_max = t[max_idx]
    
    response_strength = max_x
    
    dx = x.max()-x[0]
    response_speed = dx/time_to_max


    return response_strength, response_speed
