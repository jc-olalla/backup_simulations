"""Parametric design of a LEI kite airfoil"""
import numpy as np
from math import *
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.optimize import fsolve

# shift + alt: type multiple lines

# file directory
# change directory to your pc setup
root_path = r"/home/kaspermasure/Linux_thesis/LEI_profile_aero"

fig_file_path = Path(root_path) / 'results' / 'para_model.png'
fig_curvature_file_path = Path(root_path) / 'results' / 'para_model_curvature.png'


# Cubic polynomial
def interpolation3(P1, P2, t1, t2, n=100):
    # start and end point coordinates
    x1 = P1[0]
    x2 = P2[0]
    y1 = P1[1]
    y2 = P2[1]

    # slope at end points
    s1, s2 = np.tan(t1), np.tan(t2)

    # AX = Y
    A = np.array([[x1 ** 3, x1 ** 2, x1, 1],
                  [x2 ** 3, x2 ** 2, x2, 1],
                  [3 * x1 ** 2, 2 * x1, 1, 0],
                  [3 * x2 ** 2, 2 * x2, 1, 0]])

    Y = np.array([y1, y2, s1, s2])
    X = np.linalg.solve(A, Y)
    x_list = np.linspace(x1, x2, n)
    y_list = X[0] * x_list ** 3 + X[1] * x_list ** 2 + X[2] * x_list + X[3]

    # 2D array of the x and y coordinate of each point
    points = np.transpose(np.vstack((x_list, y_list)))
    return points


# Cubic Bezier curve
def cubic_bezier(P0, P1, P2, P3, t):
    # Bezier curve x and y coordinates
    x_bezier = (1-t)**3 * P0[0] + 3*(1-t)**2 * t * P1[0] + 3*(1-t) * t**2 * P2[0] + t**3 * P3[0]
    y_bezier = (1-t)**3 * P0[1] + 3*(1-t)**2 * t * P1[1] + 3*(1-t) * t**2 * P2[1] + t**3 * P3[1]

    # Bezier derivatives in x and y to retrieve the slope
    dx_bezier = 3 * (1-t)**2 * (P1[0]-P0[0]) + 6*(1-t)*t*(P2[0]-P1[0]) + 3*t**2*(P3[0] - P2[0])
    dy_bezier = 3 * (1-t)**2 * (P1[1]-P0[1]) + 6*(1-t)*t*(P2[1]-P1[1]) + 3*t**2*(P3[1] - P2[1])
    slope = np.array(dy_bezier/dx_bezier)

    # 2D array of the x and y coordinate of each point
    points = np.transpose(np.vstack((x_bezier, y_bezier)))
    return points, slope


# Seam angle calculator
def LE_seam_angle(tube_size, c_x, c_y):
    for angle in range(0, 120):
        seam_a = np.radians(angle)
        radius = tube_size/ 2
        s = (np.pi / 2) - seam_a
        P1 = [radius * (1 - np.cos(seam_a)), radius * np.sin(seam_a)]
        P2 = [c_x, c_y]
        poly = interpolation3(P1, P2, s, 0)
        maximum = round(max(poly[:-1, 1]), 4)   # For error prevention disregard last point

        # If highest point is below c_y add a safety margin of 3 deg
        if maximum <= c_y:
            print(angle+3)
            return np.radians(angle-10)

        if maximum >= c_y and angle == 120:
            return print("Camber lower than LE thickness")


# LEI kite profile coordinates and control points
def LEI_airfoil(seam_a, tube_size, c_x, c_y, LE_config, LE_tension, e, TE_angle, TE_cam_tension, TE_tension, LE_fillet):
    ### LE ###
    seam_a = LE_seam_angle(tube_size, c_x, c_y)          #np.radians(46)                            # LE_seam_angle(tube_size, c_x, c_y)
    radius = tube_size / 2                                                  # Radius LE tube
    P1_s = np.tan((np.pi / 2) - seam_a)                                     # Tangency at seam location
    P1 = np.array([radius * (1 - np.cos(seam_a)), radius * np.sin(seam_a)]) # Seam location
    P2 = np.array([c_x, c_y])                                               # Max camber position
    P11_max = np.array([P1[0] + (c_y-P1[1])/P1_s, P2[1]])                   # P11 at c_y height


    # Compute control points P11 and P12
    # config 1: only P11 can move, P12 at P11_max (c_y)
    # config 2: only P12 can move, P11 at P11_max (c_y)
    # config 3: P11 and P12
    if LE_config == 1:
        P11 = np.array([(1-LE_tension) * P1[0] + LE_tension * P11_max[0], (1-LE_tension) * P1[1] + LE_tension * P11_max[1]])
        P12 = P11_max

    if LE_config == 2:
        P11 = P11_max
        P12 = np.array([(P2[0] - P11_max[0]) * (1-LE_tension) + P11_max[0], P2[1]])

    if LE_config == 3:
        P11 = np.array([(1-LE_tension) * P1[0] + LE_tension * P11_max[0], (1-LE_tension) * P1[1] + LE_tension * P11_max[1]])
        P12 = np.array([(P2[0] - P11_max[0]) * (1-LE_tension) + P11_max[0], P2[1]])

    # LE bezier curve
    t_LE = np.linspace(0, 1, 80)
    LE_points = cubic_bezier(P1, P11, P12, P2, t_LE)[0]                     # 2d array of points [[x1,y1],[xn, yn]]

    # s = (np.pi / 2) - seam_a
    # LE_points = interpolation3(P1, P2, s, 0, n=80)

    # Calculate slope and curvature of upper surface
    LE_dyu_dx = np.gradient(LE_points[:, 1], LE_points[:, 0])  # First derivative (slope)
    LE_d2yu_dx2 = np.gradient(LE_dyu_dx[1:-1], LE_points[1:-1, 0])  # Second derivative (curvature)

    # print("P1-P11:", np.linalg.norm(P1 - P11))
    # print("P11-P12:", np.linalg.norm(P11 - P12))
    # print("P12-P2:", np.linalg.norm(P12 - P2))


    ### TE ###
    P2 = np.array([c_x, c_y])                                               # max camber location
    P3 = np.array([1, 0])                                                   # TE top side
    P21 = np.array([c_x + TE_cam_tension * (1-c_x), c_y])                   # TE top first control point
    D_cam_TE = sqrt((1-c_x)**2 + c_y**2)                                    # distance from max camber to TE
    reflex_angle = np.radians(TE_angle) + atan(P21[1]/(1-P21[0]))           # reflex angle
    P22 = np.array([P3[0] - D_cam_TE * TE_tension * cos(reflex_angle),
                    D_cam_TE * TE_tension * sin(reflex_angle)])             # TE top second control point, based on tension para related to the TE distance

    # TE bezier curve
    t_TE = np.linspace(0, 1, 100)
    TE_points = cubic_bezier(P2, P21, P22, P3, t_TE)[0]                     # 2d array of points [[x1,y1],[xn, yn]]


    ### TE lower side control points ###
    P5 = np.array([c_x, c_y-e])                                             # Max camber point on lower side
    P51 = np.array([c_x + TE_cam_tension * (1-c_x), c_y-e])                 # TE lower skin first control point, closest to max camber
    P4 = np.array([1 - e * sin(reflex_angle), 0 - e * cos(reflex_angle)])   # TE lower point
    P52 = np.array([P4[0] - D_cam_TE * TE_tension * cos(reflex_angle),
                    D_cam_TE * TE_tension * sin(reflex_angle) + P4[1]])     # TE lower skin 2nd control point


    ### Round TE ###
    round_TE = []                                                           # List of round TE points
    round_TE_mid = np.array([0.5 * (P4[0] + 1), 0.5 * P4[1]])               # Middle point of round TE

    # Discretize the round TE
    for i in np.linspace(0, np.pi, 30):
        round_TE_point = [round_TE_mid[0] + e/2 * sin(reflex_angle + i), round_TE_mid[1] + e/2 * cos(reflex_angle + i)]
        round_TE.append(round_TE_point)
    round_TE_points = np.array(round_TE)                                    # 2d array of points [[x1,y1],[xn, yn]]


    ### LE and TE lower side and slope at P63 ###
    LE_lower_points = LE_points.copy()                                      # Make a copy of LE_points
    LE_lower_points[:, 1] -= e                                              # Subtract e from the column Y
    LE_lower_points = LE_lower_points[:-1]                                  # Remove the last row

    t_TE_lower = np.linspace(0, 1, 100)
    TE_lower_points_init = cubic_bezier(P5, P51, P52, P4, t_TE_lower)[0]    # Initial TE lower points from P5 to P4

    # Find the index of the closest value
    lower_surface_init = np.vstack((LE_lower_points, TE_lower_points_init)) # Initial lower surface points from P4 till LE
    index_lower = np.abs(lower_surface_init[:, 0] - LE_fillet).argmin()     # Defines the intersection point P63
    TE_lower_points = lower_surface_init[index_lower:]                      # Defining TE array based on fillet percentage

    LE_fillet_slopes = cubic_bezier(P1, P11, P12, P2, t_LE)[1]              # Array of fillet slopes
    TE_lower_slopes = cubic_bezier(P5, P51, P52, P4, t_TE_lower)[1]         # Array of TE slopes
    P63_slope = np.hstack((LE_fillet_slopes, TE_lower_slopes))[index_lower] # returns the slope at P63 based on the fillet percentage


    ### LE fillet ###
    fillet_a = np.radians((40))     #np.radians((60 - 150 * LE_fillet))                # Fillet angle at intersection with LE tube based on fillet parameter
    radius = tube_size / 2

    P6_s = -np.tan((np.pi/2) - fillet_a)                                            # Slope at the fillet seam
    P6 = np.array([radius * (1 + np.cos(fillet_a)), radius * np.sin(fillet_a)])     # Fillet seam position
    P63 = np.array([TE_lower_points[0, 0], TE_lower_points[0, 1]])                  # Intersection between LE fillet and TE lower side

    D_fillet1 = 0.02 + 0.1 * LE_fillet
    D_fillet2 = 0.01 + 0.2 * LE_fillet

    # D_fillet1 = 0.02 + 0.15 * LE_fillet     # Distance of the 1st control point(closest to LE) based on the fillet %
    # D_fillet2 = 0.01 + 0.15 * LE_fillet     # Distance of the 2nd control point based on the fillet %

    # LE fillet first control point(closest to LE)
    if fillet_a < 0:
        P61 = np.array([P6[0] + D_fillet1/sqrt(1+P6_s**2), P6[1] - P6_s*(-D_fillet1/sqrt(1+P6_s**2))])
    else:
        P61 = np.array([P6[0] - D_fillet1/sqrt(1+P6_s**2), P6[1] + P6_s*(-D_fillet1/sqrt(1+P6_s**2))])

    P62 = np.array([P63[0] - D_fillet2, P63_slope * - D_fillet2 + P63[1]])      # LE fillet 2nd control point

    # LE fillet bezier curve
    t_LE_u = np.linspace(0, 1, 50)
    fillet_points = cubic_bezier(P6, P61, P62, P63, t_LE_u)[0]                  # 2d array of points [[x1,y1],[xn, yn]]


    ### LE tube ###
    circle_n_points = 80                                                        # Number of LE tube points  # int((np.rad2deg(seam_a) + np.rad2deg(fillet_a) + 180)/3)
    theta = np.linspace(np.pi - seam_a, fillet_a + np.pi*2, circle_n_points)    # Array of LE tube angles
    Origin_LE_tube = [radius, 0]                                                # Origin of the LE tube


    x_cr = Origin_LE_tube[0] + radius * np.cos(theta)                           # x location of the LE tube points
    y_cr = Origin_LE_tube[1] + radius * np.sin(theta)                           # y location of the LE tube points
    LE_tube_points = np.column_stack((x_cr, y_cr))                              # LE tube points from seam to fillet seam, 2d array of points [[x1,y1],[xn, yn]]

    both_array = np.vstack((LE_tube_points[::-1][:-1], LE_points))

    # Calculate slope and curvature of upper surface
    circ_dyu_dx = np.gradient(both_array[:, 1], both_array[:, 0])  # First derivative (slope)
    circ_d2yu_dx2 = np.gradient(circ_dyu_dx, both_array[:, 0])  # Second derivative (curvature)
    return LE_tube_points, P1, P11, P12, LE_points, TE_points, P2, P21, P22, P3, round_TE_points, P4, P5, P51, P52, TE_lower_points, P6, P61, P62, P63, fillet_points, Origin_LE_tube, round_TE_mid, seam_a, LE_dyu_dx, LE_d2yu_dx2, circ_dyu_dx, circ_d2yu_dx2, both_array


"""Boundary Layer and progression functions"""
def wall_height(Re):
    # Given constants
    y_plus = 1.0     # y+ value [-]
    rho = 1.225      # Air density [kg/m^3]

    c = 1.0          # Chord length [m]
    U_inf = Re*1.5e-5     # Free-stream velocity [m/s]

    # Dynamic viscosity [kg/m/s]
    mu = (rho * U_inf * c) / Re

    # Flat plate skin friction coefficient [-]
    c_f = 0.027 * Re**(-1/7)

    # Wall shear stress [kg/m/s^2]
    tau_w = 0.5 * c_f * rho * U_inf**2

    # Friction velocity [m/s]
    u_tau = sqrt(tau_w / rho)

    # First cell layer height [m]
    yw = (y_plus * mu) / (u_tau * rho) * 0.4
    return round(yw, 7)

def progression_f(yw, L, np):
    def F(x):
        return -yw + L * (x - 1) / (x**(np - 1) - 1)

    # Initial guess for fsolve
    x0 = 1.0001

    # Solve for x using fsolve
    x_solution = fsolve(F, x0)
    return x_solution



"""plotting"""
def plot_airfoil(fig_file_path, profile_name, LE_tube_points, P1, P11, P12, LE_points, TE_points, P2, P21, P22, P3, round_TE_points, P4, P5, P51, P52, TE_lower_points, P6, P61, P62, P63, fillet_points, seam_a):
    fig = plt.figure(figsize=(16, 6.5))
    ax = fig.add_subplot(1, 1, 1)
    ax.yaxis.grid(color='gray')
    ax.xaxis.grid(color='gray')
    ax.set_axisbelow(True)

    line_t = 2
    control_s = 50
    line_type = '-'
    control = True
    cfd_fillet = True

    # LE tube
    if cfd_fillet:
        plt.plot(LE_tube_points[:, 0], LE_tube_points[:, 1], '-o', linewidth=line_t, markersize=0.1, label='circ')

    # Plot LE full circle
    eta = np.linspace(0, 2*np.pi, 100)    # array of LE tube angles
    radius = LE_tube_points[:, 0].max() / 2
    Origin_Circle = [radius, 0]     # origin of the LE tube
    x_cr = Origin_Circle[0] + radius * np.cos(eta)    # x location of the LE tube points
    y_cr = Origin_Circle[1] + radius * np.sin(eta)    # y location of the LE tube points
    circ_full = np.column_stack((x_cr, y_cr))         # 2d array of points [[x1,y1],[xn, yn]]
    plt.plot(circ_full[:, 0], circ_full[:, 1], '--', linewidth=line_t, markersize=0.1, color='#3776ab', label="Circular tube")

    # Plot front spline
    plt.plot(LE_points[:, 0], LE_points[:, 1], line_type, color='#ff7f0e', linewidth=line_t, markersize=0.5, label="Front spline")
    if control:
        control_points = np.array([P1, P11, P12, P2])
        plt.plot(control_points[:, 0], control_points[:, 1], '--', linewidth=line_t, color='gray')
        plt.scatter(control_points[:, 0], control_points[:, 1], color='#ff7f0e', s=control_s, label="Control front")


    # Plot Rear spline
    plt.plot(TE_points[:, 0], TE_points[:, 1], line_type, color='#2CA02C', linewidth=line_t, markersize=0.5, label="Rear spline")
    if control:
        control_points = np.array([P2, P21, P22, P3])
        plt.plot(control_points[:, 0], control_points[:, 1], '--', linewidth=line_t, color='gray')
        plt.scatter(control_points[:, 0], control_points[:, 1], color='#2CA02C', s=control_s, label="Control rear")

    # Plot LE fillet
    if cfd_fillet:
        plt.plot(fillet_points[:, 0], fillet_points[:, 1], line_type, color='#D62728', linewidth=line_t, markersize=0.5, label="LE fillet")
        if control:
            control_points = np.array([P6, P61, P62, P63])
            plt.plot(control_points[:, 0], control_points[:, 1], '--', linewidth=line_t, color='gray')
            plt.scatter(control_points[:, 0], control_points[:, 1], color='#D62728', s=control_s, label="Control LE fillet")

    # Plot rear lower spline
    if cfd_fillet:
        plt.plot(TE_lower_points[:, 0], TE_lower_points[:, 1], line_type, color='teal', linewidth=line_t, markersize=0.5, label="TE lower")
        if control:
            control_points = np.array([P5, P51, P52, P4])
            plt.plot(control_points[:, 0], control_points[:, 1], '--', linewidth=line_t, color='gray')
            plt.scatter(control_points[:, 0], control_points[:, 1], color='teal', s=control_s, label="Control TE lower")

    # Plot round TE
    if cfd_fillet:
        plt.plot(round_TE_points[:, 0], round_TE_points[:, 1], '-', markersize=0.5, color='k', label="Round TE")

    plt.scatter(Origin_Circle[0], Origin_Circle[1], marker='*', color='b', s=control_s, label="LE tube centre")
    plt.scatter(TE_points[-1, 0], TE_points[-1, 1], marker='*', color='r', s=control_s, label="TE position")
    plt.scatter(LE_points[0, 0], LE_points[0, 1], marker='*', color='g', s=control_s, label="Tube-canopy intersection")
    plt.scatter(LE_points[-1, 0], LE_points[-1, 1], marker='*', color='k', s=control_s, label="Max. camber position")

    # Reflex angle for legend
    plt.scatter([0.5], [-2], marker='$r$', color='k', s=control_s-10, label="Reflex angle")

    plt.xlabel('x/c [-]')
    plt.ylabel('y/c [-]')
    plt.xticks(np.arange(0, 1.05, 0.1))
    plt.yticks(np.arange(-0.25, 0.12, 0.05))
    plt.axis('equal')
    plt.xlim([0, 1.05])          # plt.xlim([0, 1.1])
    plt.ylim([-0.25, 0.12])      # plt.ylim([-0.1, 0.20])

    plt.legend(loc=8)                # Add legend
    plt.title(profile_name, fontsize=10, pad=5)
    plt.savefig(fig_file_path, format='png', dpi=400)  # Saving figure to results
    plt.close()