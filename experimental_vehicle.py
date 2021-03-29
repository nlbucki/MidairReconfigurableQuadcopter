from model import *

# Define numerical values for the experimental vehicle
experimental_vehicle = {}

##############################################
# Inertial parameters
##############################################
# Single arm mass [kg]
experimental_vehicle[m_A] = 0.067
# Body mass [kg]
experimental_vehicle[m_B] = 0.356
# Optional payload mass [kg]
experimental_vehicle[m_D] = 0.083

# MOI of whole vehicle written at c.m. of whole vehicle C in body frame B when unfolded [kg*m^2]
experimental_vehicle[J_Sigma_unfolded] = sm.Matrix.diag([4.21e-3, 3.79e-3, 7.79e-3])

# MOI of whole vehicle written at c.m. of whole vehicle C in body frame B with two arms folded [kg*m^2]
experimental_vehicle[J_Sigma_two_arms] = sm.Matrix(
    [[4.12e-3, -1.07e-3, 0.0], [-1.07e-3, 3.33e-3, 0.0], [0.0, 0.0, 5.68e-3]]
)

# MOI of whole vehicle written at c.m. of whole vehicle C in body frame B with four arms folded [kg*m^2]
experimental_vehicle[J_Sigma_folded] = sm.Matrix.diag([3.53e-3, 2.37e-3, 3.52e-3])

# MOI of single arm at arm c.m. along each axis [kg*m^2]
experimental_vehicle[J_A_xx] = 0.000033
experimental_vehicle[J_A_yy] = 0.000068
experimental_vehicle[J_A_zz] = 0.000076

##############################################
# Geometric parameters
##############################################
# Angle of hinge off of centerline between unfolded diagonal propellers [rad]
experimental_vehicle[theta] = 0.20717058
# Position of propeller thrust axis w.r.t. arm c.m. [m]
experimental_vehicle[d_PA_A_x] = 0.0144
# Position of hinge w.r.t. arm c.m. [m]
experimental_vehicle[d_HA_A_x] = -0.0756
experimental_vehicle[d_HA_A_z] = -0.0143
# Position of central body c.m. w.r.t. hinge 1 [m]
experimental_vehicle[d_BH1_B_x] = -0.045
experimental_vehicle[d_BH1_B_y] = 0.071
experimental_vehicle[d_BH1_B_z] = -0.0017
# Position of payload c.m. w.r.t. central body c.m. when vehicle is holding payload [m]
experimental_vehicle[d_DB_B_z] = -0.157

##############################################
# Propeller parameters
##############################################
# Thrust to torque constant in forward direction [Nm/N]
experimental_vehicle[k_forward] = 0.0172
# Thrust to torque constant in reverse direction [Nm/N]
experimental_vehicle[k_reverse] = 0.038

# Minimum/maximum thrust per propeller [N]
f_min = -3.4
f_max = 7.8

##############################################
# Mapping matrices
##############################################
M_unfolded_num = M_unfolded.subs(experimental_vehicle).evalf()
M_two_arms_num = M_two_arms.subs(experimental_vehicle).evalf()
M_four_arms_tau_num = M_four_arms.subs(experimental_vehicle)[1:, :].evalf()
M_payload_num = M_payload.subs(experimental_vehicle).evalf()
