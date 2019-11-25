import numpy as np
from pprint import pprint as pp
import matplotlib.pyplot as plt

from classes.controls import controls

class ctrs():
    def __init__(self):
        self.applied = []
        self.last = None
        self.K1 = None
        self.K2 = None
        self.L = None

class ctr_sys():
    def __init__(self,A,B,C,D=0,type='LTI'):
        self.A = np.array(A)
        self.B = np.array(B)
        self.C = np.array(C)
        self.D = np.array(D)

        self.Abar = np.array(A)
        self.Bbar = np.array(B)

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

    def update_shapes(self):
        self.n = len(self.A)
        self.m = self.B.shape[1]
        self.r = len(self.C)

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

    def check_K_values(self):
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

        return K1, K2

    def plot_response(self,
                      x,    # Initial x --> x0
                      time, # range of time for plotting
                      c_point=0, # convergence point
                      z=0, # Initial z for PI controller
                      title='Untitled',
                      res=100):

        K1, K2 = self.check_K_values()

        dt = (time[1] - time[0]) / res

        X_res = np.zeros([self.n,res])
        U_res = np.zeros([self.m,res])
        Z_res = np.zeros([1,res]) if (self.ctrs.last == 'PI') else None

        x_ob = x + 0.001
        print('x_ob:')
        print(x_ob)

        r_t = range(res)

        for i in r_t:

            if self.ctrs.last == 'PI':
                val = -z
                z = z + dt*(np.subtract(np.matmul(self.C,x), c_point)[0])
                Z_res[:,i] = z
            elif self.ctrs.last == 'OB_FBBK':
                a = np.matmul(self.A,x_ob) + np.matmul(self.B,np.matmul(-K1,x_ob))
                # print('a:')
                # print(a)
                b = self.ctrs.L * np.subtract(np.matmul(self.C,x), np.matmul(self.C,x_ob))
                # print('b:')
                # print(b.transpose())
                x_ob = a + b.transpose()
                val = x - x_ob
                # print('val:')
                # print(val)
            else:
                val = c_point

            x = x + dt*(np.matmul(self.Abar,x) + self.Bbar*val)
            u = np.matmul(-K1,x) + K2 * val

            X_res[:,i] = x[:,0]
            # U_res[:,i] = u[:,0]

        # Create a plot
        plt.ion()
        # Create a figure
        fig = plt.figure()
        plot = fig.add_subplot(211)
        plot.grid(color='#9c94af', linestyle='--', linewidth=0.5)

        U_plot = fig.add_subplot(212)
        U_plot.grid(color='#9c94af', linestyle='--', linewidth=0.5)

        time_ = np.linspace(time[0],time[1],res)
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
