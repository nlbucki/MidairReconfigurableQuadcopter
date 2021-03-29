import matplotlib.pyplot as plt
import numpy as np
import scipy.spatial
from sympy import pprint

from experimental_vehicle import *

# Compute the thrust bounds that keep the vehicle in the unfolded configuration

# Attitude of vehicle (this should get cancelled out and should not appear in any final expressions)
r, p, y = sm.symbols("r p y")
R_BE = rot_axis1(r) * rot_axis2(p) * rot_axis3(y)

# Gravity (this should also get cancelled out and should not appear in any final expressions)
g = sm.Matrix([0, 0, sm.Symbol("g")])

# Write individual thrust forces in terms of the total thrust and desired torque
f_p_unfolded = M_unfolded.subs(experimental_vehicle).LUsolve(
    sm.Matrix([f_Sigma, tau_x, tau_y, tau_z])
)
f_p_two_arms = M_two_arms.subs(experimental_vehicle).LUsolve(
    sm.Matrix([f_Sigma, tau_x, tau_y, tau_z])
)
f_p = is_unfolded * f_p_unfolded + is_two_arms * f_p_two_arms

# Acceleration of vehicle c.m. written in body frame B/C
ddot_d_CE_B = R_BE * g + 1 / m_Sigma * z * f_Sigma

# Angular velocity of vehicle written in body frame B/C
# Here we assume the Coriolis term is negligible
J_Sigma_inv = (
    is_unfolded * J_Sigma_unfolded.subs(experimental_vehicle).inv()
    + is_two_arms * J_Sigma_two_arms.subs(experimental_vehicle).inv()
)
dot_omega_CE_B = J_Sigma_inv * tau_Sigma

tau_r = []
for i in range(4):
    # Reaction force acting at hinge H_i written in A_i frame
    f_r_A = (
        -m_A
        * (R_AB[i] * (ddot_d_CE_B - R_BE * g) + R_AB[i] * S(d_CA_B[i]) * dot_omega_CE_B)
        + z * f_p[i]
    )

    # Reaction torque acting at hinge H_i written in A_i frame
    tau_r.append(
        -J_A * R_AB[i] * dot_omega_CE_B
        + S(d_PA_A) * z * f_p[i]
        + z * k_p[i] * f_p[i]
        - S(d_HA_A) * f_r_A
    )

########################################################
# Print bounds for unfolded configuration
########################################################
# Extract y-component of reaction torque and sub in numerical values for unfolded configuration
print("Unfolded configuration bounds:")
unfolded_bound_matrix = sm.Matrix.zeros(4, 4)
for i in range(4):
    print("Arm {} bound (expression must be greater than zero):".format(i + 1))
    tau_r_y_num = (
        tau_r[i].subs(unfolded_config).subs(experimental_vehicle)[1, 0].evalf()
    )
    pprint(-tau_r_y_num)

    for j, term in enumerate([f_Sigma, tau_x, tau_y, tau_z]):
        unfolded_bound_matrix[i, j] = sm.Poly((-tau_r_y_num), term).coeffs()[0]

########################################################
# Print bounds for two-arms-folded configuration
########################################################
# Extract y-component of reaction torque and sub in numerical values for two arms folded configuration
print("\n----------------------------------------------------\n")
print("Two-arms-folded configuration bounds:")
two_arms_bound_matrix = sm.Matrix.zeros(4, 4)
for i in range(4):
    print("Arm {} bound (expression must be greater than zero):".format(i + 1))
    tau_r_y_num = (
        tau_r[i].subs(two_arms_folded_config).subs(experimental_vehicle)[1, 0].evalf()
    )

    if i == 0 or i == 2:
        bound = -tau_r_y_num
    else:
        bound = tau_r_y_num
    pprint(bound)

    for j, term in enumerate([f_Sigma, tau_x, tau_y, tau_z]):
        two_arms_bound_matrix[i, j] = sm.Poly(bound, term).coeffs()[0]

########################################################
# Examine vehicle agility reduction as a result of arm folding/unfolding bounds
########################################################
M_unfolded_inv = np.array(M_unfolded_num.inv()).astype(np.float64)
M_unfolded_f_sigma_tau_z_inv = np.column_stack(
    (M_unfolded_inv[:, 0], M_unfolded_inv[:, 3])
)
unfolded_bound_matrix_f_sigma_tau_z = np.column_stack(
    (unfolded_bound_matrix[:, 0], unfolded_bound_matrix[:, 3])
)
interior_point = np.array([1.0, 0.0])

plt.figure()

# Bounds on individual thrust forces (i.e. regular quadcopter constraints)
A = np.vstack((np.eye(4), -np.eye(4))).dot(M_unfolded_f_sigma_tau_z_inv)
b = np.vstack((-f_max * np.ones((4, 1)), f_min * np.ones((4, 1))))
hs = scipy.spatial.HalfspaceIntersection(np.hstack((A, b)), interior_point)
order = [0, 1, 3, 2]
intersections = [hs.intersections[i] for i in order]
plt.fill(*zip(*intersections))

# Plot region if f_min = 0
b_f_min_zero = np.vstack((-f_max * np.ones((4, 1)), np.zeros((4, 1))))
hs_f_min_zero = scipy.spatial.HalfspaceIntersection(
    np.hstack((A, b_f_min_zero)), interior_point
)
order = [0, 1, 3, 2]
intersections_f_min_zero = [hs_f_min_zero.intersections[i] for i in order]
plt.fill(*zip(*intersections_f_min_zero))

# Add bounds that prevent arms from changing configuration
A_restricted = np.vstack((A, -unfolded_bound_matrix_f_sigma_tau_z))
b_restricted = np.vstack((b, np.zeros((4, 1))))
hs_restricted = scipy.spatial.HalfspaceIntersection(
    np.hstack((A_restricted, b_restricted)), interior_point
)

intersections_restricted = [hs_restricted.intersections[i] for i in order]
plt.fill(*zip(*intersections_restricted), color="g")

plt.xlabel("f_Sigma")
plt.ylabel("tau_z")

########################################################
# Compute reduction of various control inputs at hover
########################################################
print("\n----------------------------------------------------\n")


def get_max_input(index, nominal_values, interior_point):
    A = np.vstack((np.eye(4), -np.eye(4))).dot(M_unfolded_inv[:, index])
    b = np.vstack((-f_max * np.ones((4, 1)), np.zeros((4, 1))))
    b = b + np.vstack((M_unfolded_inv, -M_unfolded_inv)).dot(nominal_values)
    max_no_folding = 9999
    for i in range(A.shape[0]):
        t = b[i] / A[i]
        if t < max_no_folding and t > 0:
            max_no_folding = t

    A_restricted = np.vstack((A.reshape(8, 1), -unfolded_bound_matrix[:, index]))
    b_restricted = np.vstack((b, -unfolded_bound_matrix * nominal_values))
    max_folding = 9999
    for i in range(A_restricted.shape[0]):
        t = b_restricted[i] / A_restricted[i]
        if t < max_folding and t > 0:
            max_folding = t

    return max_no_folding, max_folding


# Yaw torque reduction
mg = 9.81 * m_Sigma.subs(experimental_vehicle).subs(unfolded_config)
hover = np.array([[mg], [0], [0], [0]])
max_tau_z_no_folding, max_tau_z_folding = get_max_input(3, hover, 0.0)
print(
    "Max tau_z without folding bounds (f_min = 0) = {} Nm".format(max_tau_z_no_folding)
)
print("Max tau_z with folding bounds (f_min = 0) = {} Nm".format(max_tau_z_folding))
print("Percent reduction = {}\n".format(1 - max_tau_z_folding / max_tau_z_no_folding))

# Roll torque reduction
max_tau_x_no_folding, max_tau_x_folding = get_max_input(1, hover, 0.0)
print(
    "Max tau_x without folding bounds (f_min = 0) = {} Nm".format(max_tau_x_no_folding)
)
print("Max tau_x with folding bounds (f_min = 0) = {} Nm".format(max_tau_x_folding))
print("Percent reduction = {}\n".format(1 - max_tau_x_folding / max_tau_x_no_folding))

# Pitch torque reduction
max_tau_y_no_folding, max_tau_y_folding = get_max_input(1, hover, 0.0)
print(
    "Max tau_y without folding bounds (f_min = 0) = {} Nm".format(max_tau_y_no_folding)
)
print("Max tau_y with folding bounds (f_min = 0) = {} Nm".format(max_tau_y_folding))
print("Percent reduction = {}\n".format(1 - max_tau_y_folding / max_tau_y_no_folding))

# Thrust reduction
max_f_sigma_no_folding, max_f_sigma_folding = get_max_input(1, np.zeros((4, 1)), 1.0)
print(
    "Max f_sigma without folding bounds (f_min = 0) = {} Nm".format(
        max_f_sigma_no_folding
    )
)
print("Max f_sigma with folding bounds (f_min = 0) = {} Nm".format(max_f_sigma_folding))
print("Percent reduction = {}".format(1 - max_f_sigma_folding / max_f_sigma_no_folding))

plt.show()