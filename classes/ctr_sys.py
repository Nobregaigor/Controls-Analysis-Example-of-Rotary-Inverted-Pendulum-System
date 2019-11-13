import numpy as np
from pprint import pprint as pp
import matplotlib.pyplot as plt

from classes.controls import controls

class ctrs():
    def __init__(self):
        self.applied = []
        self.K1 = None
        self.K2 = None


class ctr_sys():
    def __init__(self,A,B,C,D=0,type='LTI'):
        self.A = np.array(A)
        self.B = np.array(B)
        self.C = np.array(C)
        self.D = np.array(D)

        self.type = type

        self.n = len(self.A)
        self.m = self.B.shape[1]
        self.r = len(self.C)

        self.stbly = None
        self.ctrb = None
        self.obsv = None

        self.ctrs = ctrs()

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

    def plot_response(self,x0,time,c_point=0,title='Untitled',res=100):
        dt = (time[1] - time[0]) / res
        r_t = range(res)

        time_ = np.linspace(time[0],time[1],res)
        X_res = np.zeros([self.n,res])
        U_res = np.zeros([self.m,res])
        x = x0
        
        if len(self.ctrs.applied) != 0:
            if self.ctrs.K1.any() != None and self.ctrs.K2 == None:
                K1 = self.ctrs.K1
                K2 = np.array([0])
            elif self.ctrs.K2.any() != None:
                K1 = self.ctrs.K1
                K2 = self.ctrs.K2
            else:
                K1 = np.zeros([1,self.n])
                K2 = np.array([0])
        else:
            K1 = K1 = np.zeros([1,self.n])
            K2 = np.array([0])

        for i in r_t:
            x = x + dt*(np.matmul(self.A,x) + self.B*c_point)
            X_res[:,i] = x[:,0]
            u = np.matmul(-K1,x) + K2 * c_point
            U_res[:,i] = u[:,0]

        # Create a plot
        plt.ion()
        # Create a figure
        fig = plt.figure()
        plot = fig.add_subplot(211)
        plot.grid(color='#9c94af', linestyle='--', linewidth=0.5)

        U_plot = fig.add_subplot(212)
        U_plot.grid(color='#9c94af', linestyle='--', linewidth=0.5)

        range_n = range(self.n)
        X_labels = []
        for i in range_n:
            plot.plot(time_,X_res[i,:])
            l = 'x'+str(i)
            X_labels.append(l)

        range_m = range(self.m)
        U_labels = []
        for i in range_m:
            U_plot.plot(time_,U_res[i,:])
            l = 'u'+str(i)
            U_labels.append(l)


        fig.suptitle(title,fontsize=16)
        plt.subplots_adjust(wspace = 0.2,hspace = 0.5)

        plot.legend(X_labels)
        plot.set_xlabel('time [s]')
        plot.set_ylabel('x matrix')
        plot.set_title('X response')

        U_plot.legend(U_labels)
        U_plot.set_xlabel('time [s]')
        U_plot.set_ylabel('U matrix')
        U_plot.set_title('U response')
