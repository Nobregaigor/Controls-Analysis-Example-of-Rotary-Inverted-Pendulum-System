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
    n_points = 400

    #######################################################################
    # Creating pendulum
    #######################################################################

    print('A1-A4:')
    print('The Rotary Inverted Pendulum System:\n')
    print('.'*65 + '\n')
    pendulum = ctr_sys(A, B, C)
    print('_.'*45 + '\n')


    #A5
    #Erronous pendulum
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

    # # Pole placement
    # pendulum_pp = controls.sf_pole_placement(pendulum,poles)
    #
    # # LQR Optimal solution
    pendulum_LQR = controls.sf_optimal_LQR(pendulum,Q,R)
    #
    # # Feedback Feedforward Control (needs a stabilized system)
    pendulum_fbfw = controls.fbfw(pendulum_LQR)

    #######################################################################
    # Ploting responses
    #######################################################################

    # ## Ploting the response of the open loop system
    # pendulum.plot_response(x0,t_,title='Open Loop response',res=n_points)
    #
    # ## Ploting the response of the stabilized system using Pole placement
    # pendulum_pp.plot_response(x0,t_,title='Pole Placement Stabilization',res=n_points)
    #
    # ## Ploting the response of the stabilized system using LQR
    # pendulum_LQR.plot_response(x0,t_,title='LQR Stabilization',res=n_points)
    #
    # ## Ploting the response of the command followed system using Feedback Feedforward control
    pendulum_fbfw.plot_response(x0_C,t_,c_point,title='Feedback Feedforward Control',res=n_points)
    #
    # plt.show(block=True)


    Q_PI = [[0.5,0,0,0,0],[0,0.1,0,0,0],[0,0,14,0,0],[0,0,0,5,0],[0,0,0,0,15]]
    R_PI = [[25]]

    pendulumPI = controls.sf_pi(pendulum,{'Q': Q_PI,'R': R_PI})

    x0_C = np.array([[0.],[0.],[0.],[0.],[0.]])
    pendulumPI.plot_response(x0_C,t_,1.5,title='PI Control',dyn=True,res=n_points)

    #
    A = [[0, 1],[0, 0]]
    B = [[0],[1]]
    C = [[1, 0]]

    car = ctr_sys(A, B, C)

    pI_poles = [-1.75,-5,-15,-16, -8]
    Q_PI = np.identity(3)
    R_PI = [[1]]

    x0_C = np.array([[0],[0],[0]])
    x1_c_point = x0_C[0] - 1.5

    # pendulumPI = controls.sf_pi(pendulum,pI_poles,'PP')
    carPI = controls.sf_pi(car,{'Q': Q_PI,'R':R_PI})

    carPI.plot_response(x0_C,t_,1,title='PI Control',dyn=True,res=n_points)
    #
    plt.show(block=True)











    # phi = controls.apply_method3(pendulum.A)

    # A = [[1, -3], [4, 2]]
    # B = [[1], [1]]
    # C = [1, 0]
    #
    # ss = ctr_sys(A,B,C)
    #
    # ss1 = controls.sf_pole_placement(ss,[-1,-2])
    # ss1.initialize()
    #
    # x0 = np.array([[3],[5]])
    # ss1.plot_response(x0,[0,10],100)
