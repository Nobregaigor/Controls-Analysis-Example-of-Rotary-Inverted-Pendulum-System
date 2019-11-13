import numpy as np
import numpy.linalg as la
import scipy.linalg as spla
import scipy.signal as spsg

from pprint import pprint as pp
import copy

class controls:

    #____________________________________________________________
    # Stability
    #____________________________________________________________

    def stbly(mat):
        """
            Stability function.
            Calculates system stability based on eigenvalues of system's matrix A.
            Requires:
                Numpy matrix (system A matrix)
            Returns:
             Object with two paramenters:
                .type: (string) 'gas, gs, stable, unstable,' --> value of the system stability
                .reason: (string) --> explanation for why the system has the given stability
        """
        eigVals = la.eig(mat)[0]
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

                err = abs(val - abs_eig)/val

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
                reason = '(1) Non-repeating eigvalues. (2) All eigvalues are negative.'
            elif zeros > 0 and positives == 0:
                stability = 'stable'
                reason = '(1) Non-repeating eigvalues. (2) One or more eigvalue are zero and the remaining are negative.'
            else:
                stability = 'unstable'
                reason = '(1) Non-repeating eigvalues. (2) One or more eigvalue is positive.'
        else:
            if zeros == 0 and positives == 0:
                stability = 'gs'
                reason = '(1) Repeating eigvalues. (2) All eigvalues are negative.'
            else:
                stability = 'unstable'
                reason = '(1) Repeating eigvalues. (2) One ore more eigvalues are positive.'

        class stbly_obj:
            def __init__(self,type,reason):
                self.type = type
                self.reason = reason

        return stbly_obj(stability, reason)


    #____________________________________________________________
    # Gramian
    #____________________________________________________________

    def gramian(system,type):
        """
            Gramian function.
            Calculates system gramian for Observability or Controlability
            Requires:
                System object
                string ('ctrb' or 'obsv')
            Returns:
                Numpy array: value of the requested gramian
        """
        range_n = range(system.n)
        gram = np.zeros([system.n,system.n])

        if type == 'ctrb':
            for i in range_n:
                val = np.matmul(la.matrix_power(system.A,i), system.B)
                for j in range_n:
                    gram[j][i] = val[j]
            return gram

        elif type == 'obsv':
            for i in range_n:
                val = np.matmul(system.C,la.matrix_power(system.A,i))
                for j in range_n:
                    gram[i][j] = val[j]
            return gram

    #____________________________________________________________
    # Controlability
    #____________________________________________________________

    def ctrb(system):
        """
            Controlability function.
            Calculates system gramian for Controlability
            Requires:
                System object
            Returns:
             Object with four paramenters:
                .type: (string) -> value for system Controlability
                .rank: (int) -> value for the rank of the gramian
                .val: (numpy array) -> if sys is CC returns gram,
                    otherwise returns the ctrb decomposition
                .gram: (numpy array) -> value for the ctrb gramian
        """
        class ctrb_obj:
            def __init__(self,type,rank,val,gram):
                self.type = type
                self.rank = rank
                self.val = val
                self.gram = gram

        gram = controls.gramian(system,'ctrb')
        rank = la.matrix_rank(gram)

        if rank == system.n:
            return ctrb_obj('CC',rank,gram,gram)
        else:
            val = controls.apply_decomposition(system,gram)

            A22 = val[0][system.n-1][system.n-1]
            type = 'NCC. Can make closed loop system AS.'

            if np.isscalar(A22):
                if A22.real >= 0:
                    type = 'NC. Needs more actuators.'
            else:
                eigVals = la.eig(A22)
                for val in eigVals:
                    if val.real >= 0:
                        type = 'NC. Needs more actuators.'

            return ctrb_obj(type,rank,val,gram)

    #____________________________________________________________
    # Observability
    #____________________________________________________________
    def obsv(system):
        """
            Observability function.
            Calculates system gramian for Observability
            Requires:
                System object
            Returns:
             Object with four paramenters:
                .type: (string) -> value for system Observability
                .rank: (int) -> value for the rank of the gramian
                .val: (numpy array) -> if sys is CC returns gram,
                    otherwise returns the obsv decomposition
                .gram: (numpy array) -> value for the obsv gramian
        """
        class obsv_obj:
            def __init__(self,type,rank,val,gram):
                self.type = type
                self.rank = rank
                self.val = val
                self.gram = gram

        gram = controls.gramian(system,'obsv')
        rank = la.matrix_rank(gram)
        if rank == system.n:
            return obsv_obj('CO',rank,gram,gram)
        else:
            val = controls.apply_decomposition(system,gram)
            return obsv_obj('NCO',rank,val,gram)

    #____________________________________________________________
    # Other functions
    #____________________________________________________________

    def sim_transf(system,T):
        """
            Similarity transformation function.
            Calculates the transformation of a system by a given T matrix
            Requires:
                System object
                Numpy Array (T Matrix)
            Returns:
             Tuple with three values:
                0: Ah (Numpy array) -> value for the transformed A matrix
                1: Bh (Numpy array) -> value for the transformed B matrix
                2: Ch (Numpy array) -> value for the transformed C matrix
        """
        T_inv = la.inv(T)
        Ah = np.matmul(np.matmul(T_inv, system.A),T)
        Bh = np.matmul(T_inv, system.B)
        Ch = np.matmul(system.C, T)

        return Ah, Bh, Ch

    def get_gram_T_Matrix(gram):
        """
            Function used to obtain the T matrix of a gramian matrix.
            Calculates the range space and null space of the given matrix and
                organizes it in [orth , null]
            Requires:
                Numpy Array (Gramian matrix)
            Returns:
                Numpy Array: T Matrix
        """
        orth_ = spla.orth(gram)
        null_ = spla.null_space(gram.transpose())
        T = np.hstack([orth_,null_])
        return T

    def apply_decomposition(system,gram):
        """
            Function used for the decomposition of a Matrix.
            Calculates the T matrix of the given gramian and applies the sistem
                transformation using similarity tranformation
            Requires:
                System object
                Numpy Array (Gramian matrix)
            Returns:
                value of sim_transf (see funtion)
        """
        T = controls.get_gram_T_Matrix(gram)
        return controls.sim_transf(system,T)

    def minimalize(system):
        """
            Minimilize the system using kalman decomposition.
            Calculates the kalman decomposition and applies the similarity
                transformation.
            Requires:
                System object
                Numpy Array (Gramian matrix)
            Returns:
                value of sim_transf (see funtion)
        """
        if system.ctrb.type != 'CC' or system.obsv.type != 'CO':
            ctrb_orth = spla.orth(system.ctrb.gram)
            ctrb_null = spla.null_space(system.ctrb.gram.transpose())
            obsv_orth = spla.orth(system.obsv.gram)
            obsv_null = spla.null_space(system.obsv.gram.transpose())

            n = system.n
            T = np.zeros([n,n])

            range_n = range(n)
            for i in range_n:
                if np.abs(ctrb_orth[:,i] - obsv_orth[:,i]) < 0.0001: #CO
                    for j in range_n:
                        T[j][i] = ctrb_orth[:,i]
                elif np.abs(ctrb_orth[:,i] - obsv_null[:,i]) < 0.0001: #CO_
                    for j in range_n:
                        T[j][i] = ctrb_orth[:,i]
                elif np.abs(ctrb_null[:,i] - obsv_orth[:,i]) < 0.0001: #C_O
                    for j in range_n:
                        T[j][i] = ctrb_orth[:,i]
                elif np.abs(ctrb_null[:,i] - obsv_null[:,i]) < 0.0001: #C_O_
                    for j in range_n:
                        T[j][i] = ctrb_orth[:,i]

            return sim_transf(system,T)

    #____________________________________________________________
    # State Feedback
    #____________________________________________________________

    def place_poles(system,poles):
        """
            Wrapper of the scipy.signal.place_poles
            Calculates the K matrix of a system with desired poles
            Requires:
                System object
                Array (desired poles)
            Returns:
                Numpy array: K matrix with the placed poles
        """
        K = spsg.place_poles(system.A,system.B,poles)
        return K.gain_matrix

    def sf_pole_placement(system,poles,subs=False):
        """
            State feedback with pole placement function.
            Calculates the K matrix and applys the state feedback with pole
                to the system and returns the new system.
            Requires:
                System object
                Array (desired poles)
            Returns:
                System object: new system with the applied poles.
        """
        K = controls.place_poles(system,poles)
        #make hurwitz
        Abar = np.subtract(system.A, np.matmul(system.B,K))

        if subs == False:
            new_sys = copy.deepcopy(system)
        else:
            new_sys = system

        new_sys.A = Abar
        new_sys.ctrs.applied.append('poles_placement')
        new_sys.ctrs.K1 = K
        new_sys.initialize()

        return new_sys

    def sf_optimal_LQR(system,Q,R,subs=False):

        P = spla.solve_continuous_are(system.A,system.B,Q,R)

        K = np.matmul(np.matmul(la.inv(R),system.B.transpose()),P)

        Abar = np.subtract(system.A, np.matmul(system.B,K))
        if subs == False:
            new_sys = copy.deepcopy(system)
        else:
            new_sys = system

        new_sys.A = Abar
        new_sys.ctrs.applied.append('LQR')
        new_sys.ctrs.K1 = K
        new_sys.initialize()
        return new_sys

    def fbfw(system,subs=False):

        if system.ctrb.type == 'CC':
            a = np.matmul(-system.C, la.inv(system.A))
            a = np.matmul(a, system.B)

            if len(a) == 1:
                K2 = 1/a[0]
            else:
                K2 = la.solve(a,np.identity(system.n))

            Bbar = system.B*K2

            if subs == False:
                new_sys = copy.deepcopy(system)
            else:
                new_sys = system

            new_sys.B = Bbar
            new_sys.ctrs.applied.append('fbfw')
            new_sys.ctrs.K2 = K2
            new_sys.initialize()
            return new_sys
        else:
            raise ValueError('System must be CC. Please, check input system.')




    #____________________________________________________________
    # State transition matrix
    #____________________________________________________________

    def apply_method3(Amatrix):
        """
            Function to calculate the method3.
            Calculates system state transition matrix using the 3 method
            Requires:
                Numpy array (A Matrix)
            Returns:
                Numpy array: Phi matrix
        """
        w, v = la.eig(Amatrix)

        # Calculating gamma
        # T is based on eigvalues
        print('eigvalues')
        print(w)

        n = len(w)
        gamma = np.zeros([n,n])

        n_cpx = []
        n_real = []
        r_w = range(n)
        for i in r_w:
            val = w[i]
            if np.iscomplex(val):
                n_cpx.append(i)
            elif np.isreal(val):
                n_real.append(i)

        for idx in n_cpx:
            gamma[idx][idx] = w[idx]

        for idx in n_real:
            gamma[idx][idx] = w[idx]

        print('gamma')
        print(gamma)

        # Calculating T
        # T is based on eigVectors

        print('v')
        print(v)

        n = len(v)
        T = np.zeros([n,n])

        n_cpx = []
        n_real = []
        r_v = range(n)
        for i in r_v:
            u = v[:,i]
            if np.iscomplex(u.all()):
                n_cpx.append(i)
            elif np.isreal(u.all()):
                n_real.append(i)

        for idx in n_cpx:
            u = v[:,idx]
            for j in r_v:
                T[j][idx] = u[j]

        for idx in n_real:
            u = v[:,idx]
            for j in r_v:
                T[j][idx] = u[j]

        print('T')
        print(T)

        e_gamma = np.exp(gamma)
        print('e_gamma')
        print(e_gamma)

        phi = np.matmul(T,e_gamma)
        phi = np.matmul(phi,la.inv(T))

        print('phi')
        print(phi)

        return phi
