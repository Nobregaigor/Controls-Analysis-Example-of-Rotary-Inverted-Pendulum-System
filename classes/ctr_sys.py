import numpy as np
from pprint import pprint as pp
import matplotlib.pyplot as plt

from classes.controls import controls


class ctr_sys():
    def __init__(self,A,B,C,D=0,type='LTI'):
        self.A = np.array(A)
        self.B = np.array(B)
        self.C = np.array(C)
        self.D = np.array(D)

        self.type = type

        self.n = len(self.A)
        self.m = len(self.B)
        self.r = len(self.C)

        self.stbly = None
        self.ctrb = None
        self.obsv = None

        self.initialize()


    def initialize(self):
        self.stbly = controls.stbly(self.A)
        self.ctrb = controls.ctrb(self)
        self.obsv = controls.obsv(self)
        self.print_system()

    def print_system(self):
        pp('System:')
        print('A:\n' + str(self.A))
        print('B:\n' + str(self.B))
        print('C:\n' + str(self.C))
        print('\n' + '-'*65 + '\n')

        # Checking stability
        pp('Stability: ' + self.stbly.type)
        pp('Reason: ' + self.stbly.reason)
        print('\n' + '-'*65 + '\n')

        # Check Controlability
        pp('Controlability: ' + self.ctrb.type)
        print('Controlability gramian:\n' + str(self.ctrb.gram))
        print('\n' + '-'*65 + '\n')

        # Check Observability
        pp('Observability: ' + self.obsv.type)
        print('Observability gramian:\n' + str(self.obsv.gram))
        print('\n' + '-'*65 + '\n')

    def plot_response(self,x0,time,res=100):
        dt = (time[1] - time[0]) / res

        time_ = np.linspace(time[0],time[1],res)
        res = np.zeros([self.n,res])

        x = x0
        r_t = range(100)
        for i in r_t:
            x = x + dt*(np.matmul(self.A,x))
            res[:,i] = x[:,0]

        # Create a plot
        plt.ion()
        # Create a figure
        fig = plt.figure()
        plot = fig.add_subplot(111)
        plot.grid(color='#9c94af', linestyle='--', linewidth=0.5)
        range_n = range(self.n)
        for i in range_n:
            plot.plot(time_,res[i,:])
