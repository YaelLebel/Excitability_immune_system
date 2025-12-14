# Classical Models of Immune Dynamics

This directory contains implementations of classical mathematical models of immune system dynamics, as described in the **Methods** section of the paper. These models serve as benchmarks for excitability and immune response behaviors.

## Models Included

### 1. Nowak & Bangham (1996)
**Directory:** `nowak_1996/`
- **Description:** A fundamental model of viral dynamics describing the interaction between uninfected cells, infected cells, free virus, and cytotoxic T-lymphocyte (CTL) response.
- **Key Features:** Basic predator-prey like oscillations and steady states.
- **File:** `nowak.ipynb`

### 2. Frank (2002)
**Directory:** `frank_2002/`
- **Description:** A model focusing on the regulation of the immune response, highlighting the trade-offs between effective clearance and immunopathology.
- **Key Features:** Immune response activation and contraction phases.
- **File:** `frank.ipynb`

### 3. Tan et al. (2012)
**Directory:** `tan_2012/`
- **Description:** A kinetic model of the Interferon-Beta (IFN-$\beta$) response pathway, incorporating delay differential equations (DDEs) to model detailed molecular feedbacks.
- **Key Features:** Uses `ddeint` for time-delay simulations.
- **File:** `tan.ipynb`

### 4. Excitable System (Generic)
**Directory:** `excitable_system/`
- **Description:** A reduced 2-dimensional ODE model capable of demonstrating excitable dynamics. This model captures the core phenomenological features of "all-or-none" immune responses.
- **Equations:**
  $$ \frac{dx}{dt} = x - \frac{x^3}{3} - y + I(t) $$
  $$ \frac{dy}{dt} = \epsilon (x + \beta - \gamma y) $$
- **Parameters:** $\epsilon=0.08, \gamma=0.5, \beta=0.7$.
- **Behavior:** The system exhibits a threshold for activation. Short or weak stimuli ($I(t)$) result in sub-threshold excursions, while stimuli exceeding a critical duration/magnitude trigger a large excursion (firing) before returning to the resting state.
- **File:** `excitable.ipynb`

## Usage
Each model is contained within its own subdirectory and is implemented as a Jupyter Notebook. To run a model:
1. Navigate to the specific model directory.
2. Open the `.ipynb` file using Jupyter Lab or Notebook.
3. Run all cells to reproduce the simulations and plots.

## Dependencies
- `numpy`
- `matplotlib`
- `scipy`
- `ddeint` (specifically for the Tan 2012 model)
