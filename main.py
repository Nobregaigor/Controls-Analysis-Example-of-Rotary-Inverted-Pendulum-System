import numpy as np
import numpy.linalg as la
import scipy.linalg as spla
import matplotlib.pyplot as plt

from pprint import pprint as pp

from classes.controls import controls
from classes.ctr_sys import ctr_sys

# Defining system:

if __name__ == '__main__':

    print('\n' + '#'*65 + '\n')
    print('Rotary Inverted Pendulum System Controls Analysis')
    print('Developed by Igor Nobrega')
    print('\n' + '#'*65 + '\n')

    #######################################################################
    # Defining Constants
    #######################################################################

    # Values for pendulum
    A = [[0,0,1,0],[0,0,0,1],[0,149.2751,-0.0104,0],[0,261.6091,-0.0103,0]]
    B = [[0],[0],[49.7275],[49.1494]]
    C = [[1, 1, 0, 0]]

    # Initial position
    x0 = np.array([[0.4],[0.2],[0.],[0.]])

    x0_C = np.array([[0.],[0.],[0.],[0.]])

    # Time interval
    t_ = [0, 10]

    # Values for State Feedback pole placement
    poles = [-1.75,-5,-15,-16]

    # Values for State Feedback LQR
    Q = [[5,0,0,0],[0,50,0,0],[0,0,1,0],[0,0,0,1]]
    R = [[10]]

    # Command following point
    c_point = 1.5

    x1_c_point = x0_C[0] - c_point
    print('x1-c_point = ' + str(x1_c_point))

    # Number of points used to plot responses
    n_points = 600

    #######################################################################
    # Creating pendulum
    #######################################################################

    print('A1-A4:')
    print('The Rotary Inverted Pendulum System:\n')
    print('.'*65 + '\n')
    pendulum = ctr_sys(A, B, C)
    print('_.'*45 + '\n')


    ## A5
    ## Erronous pendulum
    # Aerr = A
    # Aerr[0] = [1,0,0,0]
    #
    # print('.'*65 + '\n')
    # print('\nA5:')
    # print('Erronous Pendulum:\n')
    # print('.'*65 + '\n')
    # pendulum_err = ctr_sys(A, B, C)
    # print('_.'*45 + '\n')

    #######################################################################
    # Applying Control
    #######################################################################

    # Pole placement
    pendulum_pp = controls.sf_pole_placement(pendulum,poles)

    # LQR Optimal solution
    pendulum_LQR = controls.sf_optimal_LQR(pendulum,Q,R)

    # State Feedback , Feedback Feedforward Control
    sffbfw = {'QK': Q, 'RK': R}
    pendulum_sf_fbfw = controls.fbfw(pendulum,sffbfw,type='SF')

    # Observed Based , Feedback Feedforward Control
    poles_OB_FBFW = {'PPK': [-3, -5, -6, -8], 'PPL': [-2.66, -5, -6, -8]}
    pendulum_ob_fbfw = controls.fbfw(pendulum,poles_OB_FBFW,type='OB',k_meth='PP',ob_meth='PP')

    # State Feedback , PI Control
    poles_SF_PI = {'PPK': [-2.66, -5, -6, -8, -15]}
    pendulum_SF_PI = controls.pi(pendulum, poles_SF_PI, type='SF',k_meth='PP')

    # Output Based , PI Control
    poles_OB_PI = {'PPK': [-2.66, -5, -6, -8, -15], 'PPL': [-2.66, -5, -6, -8]}
    pendulum_OB_PI = controls.pi(pendulum, poles_OB_PI, type='OB',k_meth='PP',ob_meth='PP')


    #######################################################################
    # Ploting responses
    #######################################################################

    # Ploting the response of the open loop system
    pendulum.plot_response(x0,t_,open=True,title='Open Loop response')

    # Ploting the response of the stabilized system using Pole placement
    pendulum_pp.plot_response(x0,t_,title='Pole Placement Stabilization')

    # Ploting the response of the stabilized system using LQR
    pendulum_LQR.plot_response(x0,t_,title='LQR Stabilization')

    # Ploting the response of the command followed system using State Feedback Feedback Feedforward control
    pendulum_sf_fbfw.plot_response(x0_C,t_,c_point,title='SF FBFW Control')

    # Ploting the response of the command followed system using Output Based Feedback Feedforward control
    pendulum_ob_fbfw.plot_response(x0_C,t_,c_point,initial_ctrs_params={'x_ob': x0_C + 0.001},title='OB FBFW Control')

    # Ploting the response of the command followed system using State Feedback PI control
    pendulum_SF_PI.plot_response(x0_C,t_,c_point,initial_ctrs_params={'z': 0},title='SF PI Control')

    # Ploting the response of the command followed system using Output Based PI control
    pendulum_OB_PI.plot_response(x0_C,t_,c_point,initial_ctrs_params={'z': 0, 'x_ob': x0_C + 0.1},title='OB PI Control')

    plt.show(block=True)















#############################

# TEST:
    # A = [[0, 1],[0, 0]]
    # B = [[0],[1]]
    # C = [[1, 0]]
    #
    # car = ctr_sys(A, B, C)
    #
    # pI_poles = [-1.75,-5,-15,-16, -8]
    # Q_PI = np.identity(3)
    # R_PI = [[1]]
    #
    # x0_C = np.array([[0],[0]])
    # x1_c_point = x0_C[0] - 1.5
    #
    # # pendulumPI = controls.sf_pi(pendulum,pI_poles,'PP')
    # carPI = controls.sf_pi(car,{'Q': Q_PI,'R':R_PI})
    #
    # carPI.plot_response(x0_C,t_,1.5,title='PI Control',dyn=True,res=n_points)
    # #
    # plt.show(block=True)


################################
