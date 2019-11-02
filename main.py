import numpy as np
import numpy.linalg as la

from pprint import pprint as pp

# Defining system:

class pendulum:
    def __init__(self):
        self.A = np.array([[0,0,0,0],[0,0,0,1],[0,149.2751,-0.0104,0],[0,261.6091,-0.0103,0]])
        self.B = np.array([[0],[0],[49.7275],[49.1494]])
        self.C = np.array([1, 0, 0, 0])

class controls:
    def checkStability(mat):

        eigVals = la.eig(mat)[0]
        pp(eigVals)

        zeros = 0
        positives = 0
        repeated = False
        eAllowed = 0.001

        idx = 0
        n = len(eigVals)
        for eig in eigVals:
            idx += 1

            rest = abs(eigVals[idx:])
            for val in rest:
                abs_eig = abs(eig)
                err = abs(abs_eig - val)/abs_eig
                if err < eAllowed:
                    repeated = True

            if eig == 0:
                zeros += 1
            elif eig > 0:
                positives += 1
            else:
                pass

        if repeated == False:
            if zeros == positives == 0:
                stability = 'gas'
            elif zeros > 0 and positives == 0:
                stability = 'stable'
            else:
                stability = 'unstable'
        else:
            if zeros == 0 and positives == 0:
                stability = 'gs'
            else:
                stability = 'unstable'

        return stability


if __name__ == '__main__':

    pendulum = pendulum()
    pp(pendulum)

    stability = controls.checkStability(pendulum.A)
    pp(stability)
