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
    C = [[1, 0, 0, 0]]

    # Initial position
    x0 = np.array([[0.4],[0.2],[0.],[0.]])

    x0_C = np.array([[0.],[0.],[0.],[0.]])

    x0_D = np.array([[0.02],[0.04],[-0.02],[-0.04]])

    # Time interval
    t_ = [0, 10]

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
    #
    # print('A1-A4:')
    # print('The Rotary Inverted Pendulum System:\n')
    # print('.'*65 + '\n')
    # pendulum = ctr_sys(A, B, C)
    # print('_.'*45 + '\n')


    # A5
    # Erronous pendulum
    # Aerr = A
    # Aerr[0] = [1,0,0,0]
    #
    # print('.'*65 + '\n')
    # print('\nA5:')
    # print('Erronous Pendulum:\n')
    # print('.'*65 + '\n')
    # pendulum_err = ctr_sys(Aerr, B, C)
    # print('_.'*45 + '\n')

    #######################################################################
    # Applying Control
    #######################################################################

    # # Pole placement
    # poles = [-2.25,-5,-15,-16]
    # pendulum_pp = controls.sf_pole_placement(pendulum,poles)
    # # pendulum_pp.check_gain_margin(plot=True, x=x0,time=t_,title='PP Gain Margin')
    #
    # # LQR Optimal solution
    # pendulum_LQR = controls.sf_optimal_LQR(pendulum,Q,R)
    # # pendulum_LQR.check_gain_margin(plot=True, x=x0,time=t_,title='LQR Gain Margin')
    #
    # # # State Feedback , Feedback Feedforward Control
    # sffbfw = {'QK': Q, 'RK': R}
    # pendulum_sf_fbfw = controls.fbfw(pendulum,sffbfw,type='SF')
    #
    # # Observed Based , Feedback Feedforward Control
    # poles_OB_FBFW = {'PPK': [-3, -5, -6, -8], 'PPL': [-2.66, -5, -6, -8]}
    # pendulum_ob_fbfw = controls.fbfw(pendulum,poles_OB_FBFW,type='OB',k_meth='PP',ob_meth='PP')
    #
    # # State Feedback , PI Control
    # poles_SF_PI = {'PPK': [-2.66, -5, -6, -8, -15]}
    # pendulum_SF_PI = controls.pi(pendulum, poles_SF_PI, type='SF',k_meth='PP')
    # # pendulum_SF_PI.check_gain_margin(plot=True, x=x0_C,time=t_,title='SF PI Gain Margin')
    #
    # # Output Based , PI Control
    # poles_OB_PI = {'PPK': [-1, -5, -6, -8, -15], 'PPL': [-4, -11, -13, -24]}
    # pendulum_OB_PI = controls.pi(pendulum, poles_OB_PI, type='OB',k_meth='PP',ob_meth='PP')
    # pendulum_OB_PI.init_ctrs_storage({'z': 0, 'x_ob': x0 })
    # pendulum_OB_PI.check_gain_margin()


    #######################################################################
    # Ploting responses
    #######################################################################

    # # Ploting the response of the open loop system
    # pendulum.plot_response(x0,t_,open=True,title='Open Loop response')
    #
    # # Ploting the response of the stabilized system using Pole placement
    # pendulum_pp.plot_response(x0,t_,title='Pole Placement Stabilization')
    #
    # # Ploting the response of the stabilized system using LQR
    # pendulum_LQR.plot_response(x0,t_,title='LQR Stabilization')
    #
    # # Ploting the response of the command followed system using State Feedback Feedback Feedforward control
    # pendulum_sf_fbfw.plot_response(x0_C,t_,c_point,title='SF FBFW Control')
    #
    # # Ploting the response of the command followed system using Output Based Feedback Feedforward control
    # pendulum_ob_fbfw.plot_response(x0_D,t_,c_point,initial_ctrs_params={'x_ob': x0_C + 0.001},title='OB FBFW Control')
    #
    # # Ploting the response of the command followed system using State Feedback PI control
    # pendulum_SF_PI.plot_response(x0_C,t_,c_point,initial_ctrs_params={'z': 0},title='SF PI Control')
    #
    # # Ploting the response of the command followed system using Output Based PI control
    # X_res, U_res = pendulum_OB_PI.plot_response(x0_D,t_,c_point,initial_ctrs_params={'z': 0, 'x_ob': x0_D + 0.15},title='OB PI Control')
    # X_ob_res = pendulum_OB_PI.ctrs.res['x_ob']
    #
    # time = t_
    # res = 800
    # plt.ion()
    # # Create a figure
    # fig = plt.figure()
    # plot = fig.add_subplot(211)
    # plot.grid(color='#9c94af', linestyle='--', linewidth=0.5)
    #
    # U_plot = fig.add_subplot(212)
    # U_plot.grid(color='#9c94af', linestyle='--', linewidth=0.5)
    #
    # time_ = np.linspace(time[0],time[1],res)
    # range_n = range(pendulum_OB_PI.n)
    # X_labels = []
    # for i in range_n:
    #     plot.plot(time_,X_res[i,:])
    #     l = 'x'+str(i)
    #     X_labels.append(l)
    #
    #
    # print(X_ob_res[0])
    # for i in range_n:
    #     x__ = []
    #     for j in range(res):
    #         x__.append(X_ob_res[j][i])
    #     plot.plot(time_,x__,'--')
    #     l = 'x_ob'+str(i)
    #     X_labels.append(l)
    #
    # range_m = range(pendulum_OB_PI.m)
    # U_labels = []
    # for i in range_m:
    #     U_plot.plot(time_,U_res[i,:])
    #     l = 'u'+str(i)
    #     U_labels.append(l)
    #
    #
    # fig.suptitle('OB PI Control',fontsize=16)
    # plt.subplots_adjust(wspace = 0.2,hspace = 0.5)
    #
    # plot.legend(X_labels)
    # plot.set_xlabel('time [s]')
    # plot.set_ylabel('x matrix')
    # plot.set_title('X response')
    #
    # U_plot.legend(U_labels)
    # U_plot.set_xlabel('time [s]')
    # U_plot.set_ylabel('U matrix')
    # U_plot.set_title('U response')
    # #
    # plt.show(block=True)




    A = [[-7,1],[0,-2]]
    B = [[1],[2]]
    C = [[1, 0]]

    my_sys = ctr_sys(A, B, C)

    sffbfw = {'PPK': [-2,-4]}
    my_sys = controls.fbfw(my_sys,sffbfw,type='SF',k_meth='PP')
    print('K1: ' + str(my_sys.K1))
    print('K2: ' + str(my_sys.K2))












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
