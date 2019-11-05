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

    # creating pendulum as a control system:
    A = [[0,0,1,0],[0,0,0,1],[0,149.2751,-0.0104,0],[0,261.6091,-0.0103,0]]
    B = [[0],[0],[49.7275],[49.1494]]
    C = [1, 0, 0, 0]

    x0 = np.array([[3],[5],[2],[1]])
    p1 = [-1,-2,-3,-4]
    p2 = [-11,-7,-2,-5]

    #######################################################################
    # Creating objects
    #######################################################################

    pendulum = ctr_sys(A, B, C)

    pendulum_1 = controls.sf_pole_placement(pendulum,p1)
    pendulum_2 = controls.sf_pole_placement(pendulum,p2)

    #######################################################################
    # Ploting responses
    #######################################################################

    pendulum.plot_response(x0,[0,10],100)
    pendulum_1.plot_response(x0,[0,10],100)
    pendulum_2.plot_response(x0,[0,10],100)


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
