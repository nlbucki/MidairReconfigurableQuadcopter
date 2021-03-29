import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from experimental_vehicle import *

# Draw experimental vehicle to make sure that all dimensions look right
scale = 0.25
configurations = ['Unfolded', 'Two arms folded', 'Four arms folded', 'Holding Payload']
for config in configurations:
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlim(-scale, scale)
    ax.set_ylim(-scale, scale)
    ax.set_zlim(-scale, scale)
    ax.set_title(config)

    # B (center of mass of central body only)
    ax.scatter(0, 0, 0, color='b',label='B')

    # C
    d_CB_B_num = d_CB_B.subs(experimental_vehicle)
    if config == 'Unfolded':
        d_CB_B_num = d_CB_B_num.subs(unfolded_config)
    elif config == 'Two arms folded':
        d_CB_B_num = d_CB_B_num.subs(two_arms_folded_config)
    elif config == 'Four arms folded':
        d_CB_B_num = d_CB_B_num.subs(four_arms_folded_config)
    elif config == 'Holding Payload':
        d_CB_B_num = d_CB_B_num.subs(holding_payload_config)
    ax.scatter(d_CB_B_num[0], d_CB_B_num[1], d_CB_B_num[2], color='r', label='C')
    
    # Frames H_i, A_i, and P_i
    d_HB_B_num = [-d_BH_B[i].subs(experimental_vehicle) for i in range(4)]
    d_AB_B_num = [d_AB_B[i].subs(experimental_vehicle) for i in range(4)]
    d_PB_B_num = [d_PC_B[i].subs(experimental_vehicle) + d_CB_B_num for i in range(4)]
    quads_to_plot = [d_HB_B_num, d_AB_B_num, d_PB_B_num]
    colors = ['m', 'c', 'k']
    labels = ['H', 'A', 'P']
    for points, point_color, point_label in zip(quads_to_plot, colors, labels):
        for i in range(4):
            if config == 'Unfolded':
                points[i] = points[i].subs(unfolded_config)
            elif config == 'Two arms folded':
                points[i] = points[i].subs(two_arms_folded_config)
            elif config == 'Four arms folded':
                points[i] = points[i].subs(four_arms_folded_config)
            elif config == 'Holding Payload':
                points[i] = points[i].subs(holding_payload_config)
            pt = ax.scatter(points[i][0], points[i][1], points[i][2], color=point_color)
        pt.set_label(point_label)
    ax.legend()
        
plt.show()