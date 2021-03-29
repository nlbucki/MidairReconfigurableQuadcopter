from sympy import pprint
import scipy.linalg
import numpy as np

from experimental_vehicle import *

##########################################################
# Individual thrusts to total thrust + torque mappings
##########################################################
print("M_unfolded = ")
pprint(M_unfolded_num)

print("\nM_two_arms = ")
pprint(M_two_arms_num)

print("\nM_holding_payload = ")
pprint(M_payload_num)

# We can only control the torque with all four arms folded
print("\nM_four_arms_tau = ")
pprint(M_four_arms_tau_num)

##########################################################
# Total thrust + torque to individual thrusts mappings
##########################################################
print("\n----------------------------------------------------\n")
print("inv(M_unfolded) = ")
pprint(M_unfolded_num.inv())

print("\ninv(M_two_arms) = ")
pprint(M_two_arms_num.inv())

print("\ninv(M_payload) = ")
pprint(M_payload_num.inv())

# Here we use the pseudoinverse because we have four thrusts and only three torque components
# The pseudoinverse gives the minimum norm thrust vector for a given desired torque
print("\ninv(M_four_arms_tau) = ")
pprint(M_four_arms_tau_num.pinv())

##########################################################
# Define LQR cost matrices for attitude controller
##########################################################
I = np.diag((1.0, 1.0, 1.0))
A = np.vstack((np.hstack((np.zeros((3, 3)), I)), np.zeros((3, 6))))

# Build the cost matricies
maxRollPitchError = np.deg2rad(15.0)
maxYawError = np.deg2rad(30.0)
maxRollPitchAngVelError = 1000
maxYawAngVelError = 2.0 * maxRollPitchAngVelError
maxForwardThrustDelta = 1.0  # Newtons
maxReverseThrustDelta = maxForwardThrustDelta / 1.5  # Newtons

q_xy_att = 1 / maxRollPitchError ** 2
q_z_att = 1 / maxYawError ** 2
q_xy_ang = 1 / maxRollPitchAngVelError ** 2
q_z_ang = 1 / maxYawAngVelError ** 2
Q = np.diag((q_xy_att, q_xy_att, q_z_att, q_xy_ang, q_xy_ang, q_z_ang))

r_forward = 1 / maxForwardThrustDelta ** 2
r_reverse = 1 / maxReverseThrustDelta ** 2

##########################################################
# LQR controller for unfolded config
##########################################################
M_tau_unfolded_num = np.array(M_unfolded_num[1:, :]).astype(np.float64)
R_f = np.diag((r_forward, r_forward, r_forward, r_forward))
R_unfolded = np.linalg.pinv(M_tau_unfolded_num).T.dot(
    R_f.dot(np.linalg.pinv(M_tau_unfolded_num))
)

B_unfolded = np.vstack(
    (np.zeros((3, 3)), experimental_vehicle[J_Sigma_unfolded].inv())
).astype(np.float64)
X_unfolded = scipy.linalg.solve_continuous_are(A, B_unfolded, Q, R_unfolded)
K_unfolded = np.dot(np.linalg.inv(R_unfolded), (np.dot(B_unfolded.T, X_unfolded)))

print("\n----------------------------------------------------\n")
print("Feedback matrix K for unfolded configuration:")
pprint(sm.Matrix(K_unfolded))

##########################################################
# LQR controller for two-arms-folded config
##########################################################
M_tau_two_arms_num = np.array(M_two_arms_num[1:, :]).astype(np.float64)
R_f = np.diag((r_forward, r_reverse, r_forward, r_reverse))
R_two_arms = np.linalg.pinv(M_tau_two_arms_num).T.dot(
    R_f.dot(np.linalg.pinv(M_tau_two_arms_num))
)

B_two_arms = np.vstack(
    (np.zeros((3, 3)), experimental_vehicle[J_Sigma_two_arms].inv())
).astype(np.float64)
X_two_arms = scipy.linalg.solve_continuous_are(A, B_two_arms, Q, R_two_arms)
K_two_arms = np.dot(np.linalg.inv(R_two_arms), (np.dot(B_two_arms.T, X_two_arms)))

print("\n----------------------------------------------------\n")
print("Feedback matrix K for two-arms-folded configuration:")
pprint(sm.Matrix(K_two_arms))

##########################################################
# LQR controller for two-arms-folded with payload config
##########################################################
M_tau_payload_num = np.array(M_payload_num[1:,:]).astype(np.float64)
R_f = np.diag((r_forward, r_reverse, r_forward, r_reverse))
R_payload = np.linalg.pinv(M_tau_payload_num).T.dot(R_f.dot(np.linalg.pinv(M_tau_payload_num)))

# An approximation of the total vehicle MOI when holding the box
J_Sigma_two_arms_box = experimental_vehicle[J_Sigma_two_arms] + sm.Matrix.diag(
    [0.001, 0.001, 0.0]
)

B_two_arms_box = np.vstack((np.zeros((3, 3)), J_Sigma_two_arms_box.inv())).astype(np.float64)
X_payload = scipy.linalg.solve_continuous_are(A, B_two_arms_box, Q, R_payload)
K_payload = np.dot(np.linalg.inv(R_payload), (np.dot(B_two_arms_box.T, X_payload)))
print("\n----------------------------------------------------\n")
print("Feedback matrix K for two-arms-folded configuration with payload:")
pprint(sm.Matrix(K_payload))

##########################################################
# LQR controller for four-arms-folded config
##########################################################
# Here we ignore yaw because we only care about keeping the z-axis aligned with the gap
maxYawError = 1000.0
maxYawAngVelError = 1000.0

q_z_att = 1 / maxYawError ** 2
q_z_ang = 1 / maxYawAngVelError ** 2
Q = np.diag((q_xy_att, q_xy_att, q_z_att, q_xy_ang, q_xy_ang, q_z_ang))
R_f = np.diag((r_reverse, r_reverse, r_reverse, r_reverse))

M_tau_four_arms_num = np.array(M_four_arms_tau_num).astype(np.float64)
R_four_arms = np.linalg.pinv(M_tau_four_arms_num).T.dot(
    R_f.dot(np.linalg.pinv(M_tau_four_arms_num))
)

B_four_arms = np.vstack((np.zeros((3, 3)), experimental_vehicle[J_Sigma_folded].inv())).astype(np.float64)
X_four_arms = scipy.linalg.solve_continuous_are(A, B_four_arms, Q, R_four_arms)
K_four_arms = np.dot(np.linalg.inv(R_four_arms), (np.dot(B_four_arms.T, X_four_arms)))

print("\n----------------------------------------------------\n")
print("Feedback matrix K for four-arms-folded configuration:")
pprint(sm.Matrix(K_four_arms))