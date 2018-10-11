#============================================================
# File solver.py
# 
# QMINOS as generic LP solver
# - including DQQ of Ma et al.
#
# Laurence Yang, SBRG, UCSD
#
# 27 Apr 2016:  first version
# 05 May 2016:  ported from polytope.py from cobrame
# 10 May 2016:  standalone version
# 11 Oct 2018:  ported to solveme package
#============================================================

import numpy as np
import warnings

from qminospy import qwarmLP
from qminospy import warmLP

import time

class QMINOS(object):
    """
    Generic LP solver.
    """
    def __init__(self):
        #----------------------------------------------------
        # Solver options
        self.opt_strwidth = {}
        self.opt_realwidth = {}
        self.opt_intwidth = {}
        self.opt_strlist = {}
        self.opt_intdict = {}
        self.opt_realdict = {}
        self.opt_stropts = {}
        self.opt_intopts = {}
        self.opt_intvals = {}
        self.opt_realopts = {}
        self.opt_realvals = {}

        #----------------------------------------------------
        # LP options
        #----------------------------------------------------
        # Width of characters allowed in each options for qMINOS
        self.opt_strwidth['lp'] = 72
        self.opt_realwidth['lp'] = 55
        self.opt_intwidth['lp'] = 55
        self.opt_strlist['lp'] = [
                'Maximize',     # Default obj sense is to maximize
                'Solution No'
                ]
        self.opt_intdict['lp'] = {
                'New basis file': 11,
                'Save frequency': 500000,
                'Print level': 0,
                'Print frequency': 100000,
                'Scale option': 2,
                'Iteration limit': 2000000,
                'Expand frequency': 100000
                }
        self.opt_realdict['lp'] = {
                'Penalty parameter':100.0,
                'LU factor tol': 10.0,
                'LU update tol': 10.0,
                'LU singularity tol': 1e-30,
                'Feasibility tol': 1e-20,
                'Optimality tol': 1e-20,
                'Unbounded step size': 1e+30
                }

        #----------------------------------------------------
        # LP options: quad-precision 2 for DQQ
        #----------------------------------------------------
        self.opt_strwidth['lp_q2'] = 72
        self.opt_realwidth['lp_q2'] = 55
        self.opt_intwidth['lp_q2'] = 55
        self.opt_strlist['lp_q2'] = [
                'Maximize',     # Default obj sense is to maximize
                'Solution No'
                ]
        self.opt_intdict['lp_q2'] = {
                'New basis file': 11,
                'Save frequency': 500000,
                'Print level': 0,
                'Print frequency': 100000,
                'Scale option': 0,
                'Iteration limit': 2000000,
                'Expand frequency': 100000
                }
        self.opt_realdict['lp_q2'] = {
                'Penalty parameter':100.0,
                'LU factor tol': 10.0,
                'LU update tol': 10.0,
                'LU singularity tol': 1e-30,
                'Feasibility tol': 1e-20,
                'Optimality tol': 1e-20,
                'Unbounded step size': 1e+30
                }

        #----------------------------------------------------
        # LP options: double-precision (can't set as strict tols)
        #----------------------------------------------------
        # Width of characters allowed in each options for qMINOS
        self.opt_strwidth['lp_d'] = 72
        self.opt_realwidth['lp_d'] = 55
        self.opt_intwidth['lp_d'] = 55
        self.opt_strlist['lp_d'] = [
                'Maximize',     # Default obj sense is to maximize
                'Solution No'
                ]
        self.opt_intdict['lp_d'] = {
                'New basis file': 11,
                'Save frequency': 500000,
                'Print level': 0,
                'Print frequency': 100000,
                'Scale option': 2,
                'Iteration limit': 2000000,
                'Expand frequency': 100000
                }
        self.opt_realdict['lp_d'] = {
                'Penalty parameter':100.0,
                'LU factor tol': 1.9,
                'LU update tol': 1.9,
                'LU singularity tol': 1e-12,
                'Feasibility tol': 1e-7,
                'Optimality tol': 1e-7,
                'Unbounded step size': 1e+18
                }


    def solvelp(self, A, b, c, xl, xu, csense, precision='quad',
            x0=None, warm=False, basis=None, verbosity=1):
        """
        x,stat,hs = solvelp(self, A, b, c, xl, xu, csense, precision='quad')

        precision: 'double', 'quad', 'dq', 'dqq' (double, quad, (quad))

        Returns:
        x: optimal solution
        stat: status
        hs: optimal basis

        stat:
        0     Optimal solution found.
        1     The problem is infeasible.
        2     The problem is unbounded (or badly scaled).
        3     Too many iterations.
        4     Apparent stall.  The solution has not changed
              for a large number of iterations (e.g. 1000).

        Generic LP interface (double or quad precision)
        """
        # Reshape xl & xu if necessary
        if len(xl.shape)==1:
            xl.shape = (xl.shape[0],1)
        if len(xu.shape)==1:
            xu.shape = (xu.shape[0],1)

        J, ne, P, I, V, bl, bu = makeME_LP(A,b,c,xl,xu,csense)  

        # Linear constraints and objective
        m,n = J.shape
        ha = I
        ka = P
        ad = V
        bld = [bi for bi in bl.flat]
        bud = [bi for bi in bu.flat]
        nb = m + n

        if basis is None:
            hs = np.zeros(nb, np.dtype('i4'))   # Basis
        else:
            if len(basis) == nb:
                hs = basis
                warm=True
            else:
                hs = np.zeros(nb, np.dtype('i4'))   # Basis
                warm=False
                warnings.warn('Basis length of %d is wrong. Should be %d'%(len(basis), nb))

        if x0 is None:
            x0 = np.zeros(nb)
            
        inform = np.array(0)
        probname = 'lp'

        if verbosity > 0:
            print('Getting MINOS parameters...')

        tic = time.time()
        if precision is 'quad':
            stropts,intopts,realopts,intvals,realvals,nStrOpts,nIntOpts,nRealOpts =\
                self.get_solver_opts('lp')
            x,pi,rc = qwarmLP.qwarmlp(inform, probname, m, ha, ka, ad, bld, bud, hs, warm,
                                      stropts, intopts, realopts, intvals, realvals,
                                      nstropts = nStrOpts,
                                      nintopts = nIntOpts,
                                      nrealopts = nRealOpts)
        elif precision is 'double':
            stropts,intopts,realopts,intvals,realvals,nStrOpts,nIntOpts,nRealOpts =\
                self.get_solver_opts('lp_d')
            x,pi,rc = warmLP.warmlp(inform, probname, m, ha, ka, ad, bld, bud, hs, warm,
                                      stropts, intopts, realopts, intvals, realvals,
                                      nstropts = nStrOpts,
                                      nintopts = nIntOpts,
                                      nrealopts = nRealOpts)
        elif 'dq' in precision:
            # D
            self.opt_intdict['lp_d']['Scale option'] = 2
            stropts,intopts,realopts,intvals,realvals,nStrOpts,nIntOpts,nRealOpts =\
                self.get_solver_opts('lp_d')
            x,pi,rc = warmLP.warmlp(inform, probname, m, ha, ka, ad, bld, bud, hs, warm,
                                      stropts, intopts, realopts, intvals, realvals,
                                      nstropts = nStrOpts,
                                      nintopts = nIntOpts,
                                      nrealopts = nRealOpts)
            # Q1: pass optimal basis hs and scale = 2
            warm = True
            self.opt_intdict['lp']['Scale option'] = 2
            stropts,intopts,realopts,intvals,realvals,nStrOpts,nIntOpts,nRealOpts =\
                self.get_solver_opts('lp')
            x,pi,rc = qwarmLP.qwarmlp(inform, probname, m, ha, ka, ad, bld, bud, hs, warm,
                                      stropts, intopts, realopts, intvals, realvals,
                                      nstropts = nStrOpts,
                                      nintopts = nIntOpts,
                                      nrealopts = nRealOpts)
            # Last Q2 if requested: pass optimal basis hs and scale = 0
            if precision is 'dqq':
                self.opt_intdict['lp']['Scale option'] = 0
                stropts,intopts,realopts,intvals,realvals,nStrOpts,nIntOpts,nRealOpts =\
                    self.get_solver_opts('lp')
                x,pi,rc = qwarmLP.qwarmlp(inform, probname, m, ha, ka, ad, bld, bud, hs, warm,
                                          stropts, intopts, realopts, intvals, realvals,
                                          nstropts = nStrOpts,
                                          nintopts = nIntOpts,
                                          nrealopts = nRealOpts)
                # Kindly reset scale option to default
                self.opt_intdict['lp']['Scale option'] = 2

        else:
            raise ValueError('precision must be quad, double, dq, dqq. Provided: %s',
                    str(precision))

        time_elapsed = time.time()-tic

        self.inform = inform
        self.hs = hs
        self.x = x
        # Save dual and reduced cost information
        self.pi = pi
        self.rc = rc

        if verbosity > 0:
            print('Done in %g seconds with status %s'%(time_elapsed, str(inform)))

        return x, inform, hs


    def get_solver_opts(self, prob='lp'):
        """
        Return options that will be passed as arguments to minoss
        """
        nStrOpts = len(self.opt_strlist[prob])
        nRealOpts = len(self.opt_realdict[prob].keys())
        nIntOpts = len(self.opt_intdict[prob].keys())

        stropts  = np.array(
                np.array([c for c in [s.ljust(self.opt_strwidth[prob]) for 
                                      s in self.opt_strlist[prob]]],
                    dtype='c').T)
        self.opt_stropts[prob] = stropts

        intkeys = self.opt_intdict[prob].keys()
        realkeys = self.opt_realdict[prob].keys()

        intopts  = np.array(
                np.array([c for c in [s.ljust(self.opt_intwidth[prob]) for s in intkeys]],
                    dtype='c').T)
        self.opt_intopts[prob] = intopts

        realopts = np.array(
                np.array([c for c in [s.ljust(self.opt_realwidth[prob]) for s in realkeys]],
                    dtype='c').T)
        self.opt_realopts[prob] = realopts

        intvals  = np.array([self.opt_intdict[prob][k] for k in intkeys], dtype='i4')
        realvals = np.array([self.opt_realdict[prob][k] for k in realkeys], dtype='d')

        self.opt_intvals[prob] = intvals
        self.opt_realvals[prob]= realvals

        return stropts, intopts, realopts, intvals, realvals, nStrOpts, nIntOpts, nRealOpts

    def set_realopts(self, model, param_dict):
        """
        Set real valued solver options
        model:  'lp', 'lp_d','lp_q2','qp', 'nlp'
        param_dict: {param: val}
                    param:  name of param value
                    val:    real-valued value
        """
        for param,val in param_dict.iteritems():
            if self.opt_realdict[model].has_key(param):
                if isinstance(val, float):
                    self.opt_realdict[model][param]=val
                else:
                    warnings.warn('val must be float')
            else:
                warnings.warn('Unrecognized solver option:'+param)

    def set_intopts(self, model, param_dict):
        """
        Set int-valued solver options
        model:  'lp', 'lp_d','lp_q2','qp', 'nlp'
        param_dict: {param: val}
                    param:  name of param value
                    val:    int-valued value
        """
        for param,val in param_dict.iteritems():
            if self.opt_intdict[model].has_key(param):
                if isinstance(val, int):
                    self.opt_intdict[model][param]=val
                else:
                    warnings.warn('val must be int')
            else:
                warnings.warn('Unrecognized solver option:'+param)


def makeME_LP(S,b,c,xl,xu,csense):
    """
    Create simple LP for qMINOS and MINOS

    12 Aug 2015: first version
    """
    import numpy as np
    import scipy as sp
    import scipy.sparse as sps
    import time

    # c is added as a free (unbounded slacks) row, 
    # so that MINOS treats problem as an LP - Ding Ma
    J = sps.vstack((
        S,
        c)
        ).tocsc()
    J.sort_indices()
    if hasattr(b,'tolist'):
        b = b.tolist()
    b2 = b + [0.0]
    m,n = J.shape
    ne = J.nnz
    # Finally, make the P, I, J, V, as well
    # Row indices: recall fortran is 1-based indexing
    I = [i+1 for i in J.indices]
    V = J.data
    # Pointers to start of each column
    # Just change to 1-based indexing for Fortran
    P = [pi+1 for pi in J.indptr]

    # Make primal and slack bounds
    bigbnd =   1e+40
    # For csense==E rows (equality)
    sl     =   np.matrix([bi for bi in b2]).transpose()
    su     =   np.matrix([bi for bi in b2]).transpose()
    for row,csen in enumerate(csense):
        if csen == 'L':
            sl[row] = -bigbnd
        elif csen == 'G':
            su[row] = bigbnd
        elif csen == 'E':
            pass
        else:
            raise ValueError('Unrecognized csense: %s'%csen)
    # Objective row has free bounds
    sl[m-1] = -bigbnd
    su[m-1] = bigbnd

    bl = sp.vstack([xl, sl])
    bu = sp.vstack([xu, su])

    return J, ne, P, I, V, bl, bu
