# MidairReconfigurableQuadcopter
Supplementary code for analyzing vehicles as described in the paper entitled "Design and Control of a Midair Reconfigurable  Quadcopter using Unactuated Hinges".

The following python scripts are provided:
- `model_visualization.py`: Displays the positions of each of the relevant points (B, Ai, H, C, P) in 3D for the experimental vehicle in each configuration. This is useful for checking that the geometry of the vehicle is defined correctly.
- `control_input_bounds.py`: Computes the bounds on the control inputs (total thrust + torque) that ensure the arms do not fold/unfold unexpectedly. Values are computed for the experimental vehicle in the unfolded and two-arms-folded configurations. This also displays a graph and some additional computations that show how the set of feasible control inputs is reduced by enforcing these bounds.
- `controller_synthesis.py`: Computes and displays the mapping matrices (and their inverses) between the individual thrust forces and the total thrust + torque commands used to control the vehicle in each configuration. Also computes the infinite-horizon LQR state feedback matrix for each configuration. All values are computed for the experimental vehicle.

These scripts utilize the following files:
- `model.py`: Defines the model (e.g. geometry, mass, inertia, etc.) of the vehicle using sympy. This is a general model that can be used to compute symbolic expressions in terms of the various parameters of the vehicle.
- `experimental_vehicle.py`: Defines the numeric values for each parameter used in `model.py` for the experimental vehicle. Future users will likely want to modify these values to reflect the properties of their own vehicles. 

Dependencies: numpy, sympy, scipy, matplotlib
