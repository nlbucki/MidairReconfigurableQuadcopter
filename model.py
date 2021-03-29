import sympy as sm
from sympy.matrices import rot_axis1, rot_axis2, rot_axis3

# Define skew symmetric matrix form for convenience
def S(vec):
    assert len(vec) == 3
    return sm.Matrix([[0, -vec[2], vec[1]], [vec[2], 0, -vec[0]], [-vec[1], vec[0], 0]])


# Indicator variables (0 or 1) that we set depending on the configuration of the vehicle
# Each indicator variable indicates whether the corresponding arm is folded (1) or unfolded (0)
fold_1, fold_2, fold_3, fold_4 = sm.symbols("fold_1 fold_2 fold_3 fold_4")
has_payload = sm.Symbol("has_payload")

# Multipliers defined for convenience
is_unfolded = (1 - fold_1) * (1 - fold_2) * (1 - fold_3) * (1 - fold_4)
# Assumes we fold arms 2 and 4
is_two_arms = (1 - fold_1) * fold_2 * (1 - fold_3) * fold_4
is_folded = fold_1 * fold_2 * fold_3 * fold_4

# Used with 'subs' at the end
unfolded_config = {fold_1: 0, fold_2: 0, fold_3: 0, fold_4: 0, has_payload: 0}
two_arms_folded_config = {fold_1: 0, fold_2: 1, fold_3: 0, fold_4: 1, has_payload: 0}
four_arms_folded_config = {fold_1: 1, fold_2: 1, fold_3: 1, fold_4: 1, has_payload: 0}
holding_payload_config = {fold_1: 0, fold_2: 1, fold_3: 0, fold_4: 1, has_payload: 1}

# Control inputs (these are mapped to individual thrust forces in a way that depends on the configuration)
f_Sigma = sm.Symbol("f_Sigma")
tau_x, tau_y, tau_z = sm.symbols("tau_x tau_y tau_z")
tau_Sigma = sm.Matrix([tau_x, tau_y, tau_z])

# Inertial parameters
# Mass of one arm, mass of central body, mass of optional payload
m_A, m_B, m_D = sm.symbols("m_A m_B m_D")
m_Sigma = 4 * m_A + m_B + has_payload * m_D

# Inertia of entire vehicle. This changes depending on which configuration we're in.
J_Sigma_unfolded = sm.MatrixSymbol("J_Sigma_unfolded", 3, 3)
J_Sigma_two_arms = sm.MatrixSymbol("J_Sigma_two_arms", 3, 3)
J_Sigma_folded = sm.MatrixSymbol("J_Sigma_folded", 3, 3)

# Inertia of one arm (assume diagonal)
J_A_xx, J_A_yy, J_A_zz = sm.symbols("J_A_xx J_A_yy J_A_zz")
J_A = sm.Matrix.diag([J_A_xx, J_A_yy, J_A_zz])

# Angle offset from the diagonal. theta = 0 for a vehicle without angled arms
theta = sm.Symbol("theta")
# Rotations relating body-fixed frame B to hinge frames H_i (these aren't explicity
# defined in the paper, but we define them here for convenience)
R_HB_1 = rot_axis3(-sm.pi / 4 + theta)
R_HB_2 = rot_axis3(-sm.pi / 4 - sm.pi / 2 - theta)
R_HB_3 = rot_axis3(-sm.pi / 4 - sm.pi + theta)
R_HB_4 = rot_axis3(-sm.pi / 4 - 3 * sm.pi / 2 - theta)
R_HB = [R_HB_1, R_HB_2, R_HB_3, R_HB_4]

# Rotations relating hinge frames H_i to arm frames A_i
R_AH_1 = rot_axis2(fold_1 * sm.pi / 2)
R_AH_2 = rot_axis2(fold_2 * sm.pi / 2)
R_AH_3 = rot_axis2(fold_3 * sm.pi / 2)
R_AH_4 = rot_axis2(fold_4 * sm.pi / 2)
R_AH = [R_AH_1, R_AH_2, R_AH_3, R_AH_4]

# Rotation between arm frames A_i and frame B
R_AB = [R_AH[i] * R_HB[i] for i in range(4)]
R_BA = [R_AB[i].T for i in range(4)]

# Vector from arm center of mass A_i to propeller thrust axis P_i written in arm frame A_i
# (assumes that the propeller thrust axis lies in the X-Z plane of the arm frame for each arm)
d_PA_A_x = sm.Symbol("d_PA_A_x")
d_PA_A = sm.Matrix([d_PA_A_x, 0, 0])
d_PA_B = [R_BA[i] * d_PA_A for i in range(4)]

# Vector from arm center of mass A_i to hinge H_i written in arm frame A_i
# (assumes hinge lies on the X-Z plane of the arm frame for each arm)
d_HA_A_x, d_HA_A_z = sm.symbols("d_HA_A_x d_HA_A_z")
d_HA_A = sm.Matrix([d_HA_A_x, 0, d_HA_A_z])
d_HA_B = [R_BA[i] * d_HA_A for i in range(4)]

# Vector from hinge H_i to central body center of mass B written in the body frame B
# Note that the position of the central body center of mass B is not the same as the position of center
# of mass of the vehicle C (C is the center of mass of the central body + the four arms, but B is just
# the center of mass of the central body).
# We define the vectors for arms 2-4 based on the vector for arm 1
d_BH1_B_x, d_BH1_B_y, d_BH1_B_z = sm.symbols("d_BH1_B_x d_BH1_B_y d_BH1_B_z")
d_BH_B_1 = sm.Matrix([d_BH1_B_x, d_BH1_B_y, d_BH1_B_z])
d_BH_B_2 = sm.Matrix([-d_BH1_B_x, d_BH1_B_y, d_BH1_B_z])
d_BH_B_3 = sm.Matrix([-d_BH1_B_x, -d_BH1_B_y, d_BH1_B_z])
d_BH_B_4 = sm.Matrix([d_BH1_B_x, -d_BH1_B_y, d_BH1_B_z])
d_BH_B = [d_BH_B_1, d_BH_B_2, d_BH_B_3, d_BH_B_4]

# Vector from central body c.m. B to payload c.m. D written in body frame B/C
d_DB_B_z = sm.Symbol("d_DB_B_z")
d_DB_B = sm.Matrix([0, 0, d_DB_B_z])

# Position of each arm c.m. A_i w.r.t central body c.m. B written in B frame
d_AB_B = [-d_HA_B[i] - d_BH_B[i] for i in range(4)]

# Position of vehicle c.m. C w.r.t. central body c.m. B written in B frame
d_CB_B = (
    m_A * (d_AB_B[0] + d_AB_B[1] + d_AB_B[2] + d_AB_B[3]) + has_payload * m_D * d_DB_B
) / m_Sigma

# Position of each arm c.m A_i w.r.t. vehicle c.m. written in B frame
d_CA_B = [d_CB_B - d_AB_B[i] for i in range(4)]

# Position of propellers P_i w.r.t. vehicle c.m. C written in body frame B
d_PC_B = [d_PA_B[i] + d_AB_B[i] - d_CB_B for i in range(4)]

# Propeller constants relating thrust to torque (torque = k * thrust)
# Assume k is the same for each propeller
k_forward = sm.Symbol("k_forward")
k_reverse = sm.Symbol("k_reverse")

# Use k_forward when arm is unfolded and k_reverse when arm is folded
# Change the sign of k depending on the handedness of the propeller
k_p = [0, 0, 0, 0]
k_p[0] = -((1 - fold_1) * k_forward + fold_1 * k_reverse)
k_p[1] = (1 - fold_2) * k_forward + fold_2 * k_reverse
k_p[2] = -((1 - fold_3) * k_forward + fold_3 * k_reverse)
k_p[3] = (1 - fold_4) * k_forward + fold_4 * k_reverse

# Compute the thrust conversion matrix for different configurations
z = sm.Matrix([0, 0, 1])
M_f_sigma = sm.Matrix.zeros(1, 4)
M_tau = sm.Matrix.zeros(3, 4)
for i in range(4):
    thrust_dir_i = R_BA[i] * z  # Rotate thrust direction into body frame
    M_f_sigma[0, i] = z.dot(thrust_dir_i)
    M_tau[:, i] = d_PC_B[i].cross(thrust_dir_i) + k_p[i] * thrust_dir_i
M = sm.Matrix.vstack(M_f_sigma, M_tau)

M_unfolded = M.subs(unfolded_config)
M_two_arms = M.subs(two_arms_folded_config)
M_four_arms = M.subs(four_arms_folded_config)
M_payload = M.subs(holding_payload_config)
