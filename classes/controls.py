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
        # Get eigvalues of matrix mat
        eigVals = la.eig(mat)[0]
        # Initiate counters and Booleans
        zeros = 0
        positives = 0
        repeated = False
        # Max error allowed for the given task
        eAllowed = 0.001
        # Iterate through the eigvalues
        idx = 0
        n = len(eigVals)
        for eig in eigVals:
            idx += 1
            # "Recursive" approach that checks other values in array with respect to the one being analyzed
            rest = abs(eigVals[idx:])
            # Applied logic to increase counters
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
        # Analysing counters and determining system stability
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

        # Creating object to hold stability nformation
        class stbly_obj:
            def __init__(self,type,reason):
                self.type = type
                self.reason = reason

        # Returning stability object
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
                    gram[i][j] = val[0][j] # KEEEP A LOOK HERE, MIGHT AFFECT STUFF LATER

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
        # Create a class to store ctrb information
        class ctrb_obj:
            def __init__(self,type,rank,val,gram):
                self.type = type
                self.rank = rank
                self.val = val
                self.gram = gram
        # Calculate the controlability gramian
        gram = controls.gramian(system,'ctrb')
        # Get the rank of the controlability gramian
        rank = la.matrix_rank(gram)
        # Determine system controlability
        if rank == system.n:
            return ctrb_obj('CC',rank,gram,gram)
        # If rank != system.n, apply decomposition
        else:
            val = controls.apply_decomposition(system,gram,'ctrb')
            # Check the value of the lowest right element of the upper triangle matrix
            # of the decomposed A matrix

            A22 = val[0][system.n-1][system.n-1]
            type = 'NCC. But, can make closed loop system AS.'
            # If the value or its eigvalues are larger than 0, system is not controllable
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

        # Create a class to store obsv information
        class obsv_obj:
            def __init__(self,type,rank,val,gram):
                self.type = type
                self.rank = rank
                self.val = val
                self.gram = gram

        # Calculate the observability gramian
        gram = controls.gramian(system,'obsv')
        # Check the rank of the obsv gramian
        rank = la.matrix_rank(gram)
        # Determine the observability of the system
        if rank == system.n:
            return obsv_obj('CO',rank,gram,gram)
        # If the system is not CO, apply decomposition and return its value
        else:
            val = controls.apply_decomposition(system,gram,'obsv')
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

    def get_gram_T_Matrix(gram,type):
        """
            Function used to obtain the T matrix of a gramian matrix.
            Calculates the range space and null space of the given matrix and
                organizes it in [orth , null]
            Requires:
                Numpy Array (Gramian matrix)
            Returns:
                Numpy Array: T Matrix
        """
        if type == 'ctrb':
            orth_ = spla.orth(gram)
            null_ = spla.null_space(gram.transpose())
        elif type == 'obsv':
            orth_ = spla.orth(gram.transpose())
            null_ = spla.null_space(gram)

        T = np.hstack([orth_,null_])
        return T

    def apply_decomposition(system,gram,type):
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
        T = controls.get_gram_T_Matrix(gram,type)
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

    def make_hurwitz(A,B,K):
        return np.subtract(A, np.matmul(B,K))

    #____________________________________________________________
    # Pole Placement
    #____________________________________________________________

    def place_poles(A,B,poles):
        """
            Wrapper of the scipy.signal.place_poles
            Calculates the K matrix of a system with desired poles
            Requires:
                System object
                Array (desired poles)
            Returns:
                Numpy array: K matrix with the placed poles
        """
        K = spsg.place_poles(A,B,poles)
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
        # Place the poles on the system and retrieve the K matrix
        K = controls.place_poles(system.A,system.B,poles)

        # Define the behavior of the controller u
        def u(self,x,c_point=0,dt=0.001):
            u = - np.matmul(self.K1,x)
            return u

        # Assign the values for the system.
        new_sys = system if subs else copy.deepcopy(system)
        new_sys.Abar = controls.make_hurwitz(system.A,system.B,K)

        new_sys.K1 = K
        new_sys.u = u

        new_sys.ctrs.applied.append('poles_placement')
        new_sys.ctrs.last = 'poles_placement'

        return new_sys

    #____________________________________________________________
    # LQR
    #____________________________________________________________

    def LQR(A,B,Q,R):
        # Solve the continuous-time algebraic Riccati equation
        P = spla.solve_continuous_are(A,B,Q,R)
        # Obtain K matrix
        K = np.matmul(np.matmul(la.inv(R),B.transpose()),P)
        return K

    def sf_optimal_LQR(system,Q,R,subs=False):

        # Compute K matrix given A, B, Q and R
        K = controls.LQR(system.A,system.B,Q,R)

        # Define the behavior of the controller u
        def u(self,x,c_point=0,dt=0.001):
            u = - np.matmul(self.K1,x)
            return u

        # Assign the values for the system.
        new_sys = system if subs else copy.deepcopy(system)
        new_sys.Abar = controls.make_hurwitz(system.A,system.B,K)

        new_sys.K1 = K
        new_sys.u = u

        new_sys.ctrs.applied.append('LQR')
        new_sys.ctrs.last = 'LQR'

        return new_sys

    #____________________________________________________________
    # Feedback / Feedforward SF and OB Control
    #____________________________________________________________

    def fbfw(system,values,type='SF',k_meth='LQR',ob_meth='LQR',subs=False):

        # Find K1 value based on chosen stabilization method
        if k_meth =='PP':
            K1 = controls.place_poles(system.A,system.B,values['PPK'])

        elif k_meth =='LQR':
            K1 = controls.LQR(system.A,system.B,values['QK'],values['RK'])
        else:
            raise ValueError('K Stabilization method not understood. Please verify.')

        # Make the system hurwitz
        Abar = controls.make_hurwitz(system.A,system.B,K1)
        a = np.matmul(np.matmul(-system.C, la.inv(Abar)),system.B)
        # Find K2 by solving the Hurwitz system and making it equal to I
        if len(a) == 1:
            K2 = 1/a[0]
        else:
            K2 = la.solve(a,np.identity(system.n))

        # finding L values (if it Observed Based Control)
        if type == 'OB':
            if ob_meth =='PP':
                L = controls.place_poles(system.A.transpose(),system.C.transpose(),values['PPL'])
            elif ob_meth =='LQR':
                L = controls.LQR(system.A.transpose(),system.C.transpose(),values['QL'],values['RL'])
            else:
                raise ValueError('L Stabilization method not understood. Please verify.')

            # Define the observed based control u signal behavior
            def u(self,x,c_point=0,dt=0.001):

                u = - np.matmul(self.K1,new_sys.ctrs.strg['x_ob']) + self.K2 * c_point
                new_sys.ctrs.strg['x_ob'] = new_sys.ctrs.strg['x_ob'] + dt*(np.matmul(self.A,new_sys.ctrs.strg['x_ob']) + np.matmul(self.B,u) + self.L * np.subtract(np.matmul(self.C,x),np.matmul(self.C,new_sys.ctrs.strg['x_ob'])))
                new_sys.ctrs.strg['e'] = x - new_sys.ctrs.strg['x_ob']

                # Saving response in case it is needed later or it needs to be plotted
                new_sys.ctrs.res['x_ob'].append(new_sys.ctrs.strg['x_ob'])
                new_sys.ctrs.res['e'].append(new_sys.ctrs.strg['e'])

                return u

        elif type == 'SF':
            # Define the state feedback control u signal behavior
            def u(self,x,c_point=0,dt=0.001):
                u = - np.matmul(self.K1,x) + self.K2 * c_point
                return u
        else:
            raise ValueError('FBFW type not understood. Please verify, options are: SF or OB')

        # Assign the values for the system.
        new_sys = system if subs else copy.deepcopy(system)
        new_sys.Abar = Abar
        new_sys.Bbar = system.B*K2

        new_sys.K = [K1, K2]
        new_sys.K1 = K1
        new_sys.K2 = K2
        new_sys.L = L.transpose() if type == 'OB' else None
        new_sys.u = u

        new_sys.ctrs.applied.append(['OB_FBBK'])
        new_sys.ctrs.strg = {'x_ob': 0, 'e': 0 }
        new_sys.ctrs.res = {'x_ob': [], 'e': []}
        new_sys.ctrs.last = 'OB_FBBK'


        return new_sys

    #____________________________________________________________
    # PI Control
    #____________________________________________________________

    def pi(system,values,type='SF',k_meth='LQR',ob_meth='LQR',subs=False):
        # def new system:
        F = np.c_[ np.vstack([system.A, system.C]), np.zeros(system.A.shape[0] + system.C.shape[0])]
        G = np.r_[ system.B, np.zeros([system.B.shape[1],system.B.shape[1]])]
        H = np.r_[ np.zeros([system.B.shape[0],system.B.shape[1]]), -np.identity(system.B.shape[1])]

        #check controlability
        c_matrix = np.vstack((np.hstack((system.A,system.B)), np.hstack((system.C, np.zeros([system.C.shape[0],system.B.shape[1]])))))

        if la.matrix_rank(c_matrix) == system.n + system.m:
            sys_ctrb = 'CC'
        else:
            sys_ctrb = 'NCC'

        # Find K values
        if k_meth =='PP':
            K = controls.place_poles(F,G,values['PPK'])
        elif k_meth =='LQR':
            K = controls.LQR(F,G,values['QK'],values['RK'])
        else:
            raise ValueError('Stabilization method not understood. Please verify.')

        K1 = np.array([K[0][:system.n]])
        K2 = np.array([K[0][system.n:]])

        # Calculate response based on required method type
        # If Obrserved based
        if type == 'OB':
            # Find L values
            if ob_meth =='PP':
                L = controls.place_poles(system.A.transpose(),system.C.transpose(),values['PPL'])
            elif ob_meth =='LQR':
                L = controls.LQR(system.A.transpose(),system.C.transpose(),values['QL'],values['RL'])
            else:
                raise ValueError('L Stabilization method not understood. Please verify.')

            # Define the observed based control u signal behavior
            def u(self,x,c_point=0,dt=0.001):

                u = - np.matmul(self.K1,new_sys.ctrs.strg['x_ob']) - self.K2 * new_sys.ctrs.strg['z']
                new_sys.ctrs.strg['z'] = new_sys.ctrs.strg['z'] + dt*(np.subtract(np.matmul(self.C,x), c_point)[0])
                new_sys.ctrs.strg['x_ob'] = new_sys.ctrs.strg['x_ob'] + dt*(np.matmul(self.A,new_sys.ctrs.strg['x_ob']) + np.matmul(self.B,u) + self.L * np.subtract(np.matmul(self.C,x),np.matmul(self.C,new_sys.ctrs.strg['x_ob'])))
                new_sys.ctrs.strg['e'] = x - new_sys.ctrs.strg['x_ob']

                # Saving response in case it is needed later or it needs to be plotted
                new_sys.ctrs.res['z'].append(new_sys.ctrs.strg['z'])
                new_sys.ctrs.res['x_ob'].append(new_sys.ctrs.strg['x_ob'])
                new_sys.ctrs.res['e'].append(new_sys.ctrs.strg['e'])

                return u

        # If State Feedback
        elif type == 'SF':
            # Define the state feedback control u signal behavior
            def u(self,x,c_point=0,dt=0.001):
                u = - np.matmul(self.K1,x) - self.K2 * new_sys.ctrs.strg['z']
                new_sys.ctrs.strg['z'] = new_sys.ctrs.strg['z'] + dt*(np.subtract(np.matmul(self.C,x), c_point)[0])
                # Saving response in case it is needed later or it needs to be plotted
                new_sys.ctrs.res['z'].append(new_sys.ctrs.strg['z'])
                return u

        else:
            raise ValueError('FBFW type not understood. Please verify, options are: SF or OB')

        new_sys = system if subs else copy.deepcopy(system)
        new_sys.Abar = controls.make_hurwitz(system.A,system.B,K1)
        new_sys.Bbar = system.B*K2

        new_sys.K = K
        new_sys.K1 = K1
        new_sys.K2 = K2
        new_sys.L = L.transpose() if type == 'OB' else None
        new_sys.u = u

        new_sys.ctrs.last = 'PI'
        new_sys.ctrs.applied.append(['poles_placement','PI'])
        new_sys.ctrs.strg = {'z': 0, 'x_ob': 0, 'e': 0 }
        new_sys.ctrs.res = {'z': [], 'x_ob': [], 'e': []}

        # new_sys.ctrb = sys_ctrb

        return new_sys

# OLD:
# def fbfw(system,subs=False):
#
#     if system.ctrb.type == 'CC':
#
#         Abar = controls.make_hurwitz(system.A,system.B,system.K1)
#         a = np.matmul(np.matmul(-system.C, la.inv(Abar)),system.B)
#
#         if len(a) == 1:
#             K2 = 1/a[0]
#         else:
#             K2 = la.solve(a,np.identity(system.n))
#
#         def u(self,x,c_point=0,dt=0.001):
#             u = -np.matmul(self.K1,x) + self.K2 * c_point
#             return u
#
#         new_sys = system if subs else copy.deepcopy(system)
#         new_sys.Bbar = system.B*K2
#         new_sys.K2 = K2
#         new_sys.u = u
#
#         new_sys.ctrs.applied.append('FBFW')
#         new_sys.ctrs.last = 'FBFW'
#         # new_sys.initialize()
#         return new_sys
#     else:
#         raise ValueError('System must be CC. Please, check input system.')
