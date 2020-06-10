from __future__ import print_function, absolute_import
#============================================================
# File me2.py
#
#   class     ME_NLP
#
# Class and methods for creating LP/NLP problems from ME
# object, and solving using Quad MINOS
#
# Laurence Yang, SBRG, UCSD
#
# 20 Jul 2015:  first version. Python port of dumpME_NLP.m and
#               dumpMat.m
# 23 Oct 2015:  updated LP/FVA to use non-zero b (met._bound)
#               also updated LP/FVA to use met._constraint_sense
# 18 May 2016:  moved make_dilution_fluxes from ME_NLP1 to ME_NLP
#               (ME_NLP1 inherits)
# 05 Oct 2016:  renamed qnonlinme to me2. Keeping qnonlinme.py
#               for backwards compatibility, where we just import *
#               from me2.py
#============================================================

import numpy as np
import scipy.sparse as sps
from cobra.core.Solution import Solution
from cobra import DictList
from cobra import Reaction, Metabolite
from sympy import Basic
import time
import warnings
import cobrame
from qminospy import qwarmLP
from qminospy import warmLP
from cobrame import MEModel
import re
import six


class ME_NLP:
    """
    Contains the data matrices needed for solving ME as an NLP using qMINOS
    """
    def __init__(self, me, growth_key='mu', growth_rxn='biomass_dilution'):
        # The ME model object
        self.me = me
        self.growth_symbol = growth_key
        self.growth_key = growth_key
        self.growth_rxn = growth_rxn
        self.scaleUnits = False
        self.typeM = [cobrame.MetabolicReaction]
        self.typeE = [cobrame.TranscriptionReaction,
                      cobrame.TranslationReaction, 
                      cobrame.ComplexFormation]
        self.unitDict = {
                'e_mult': 1e-6,
                'typeE': self.typeE,
                'typeM': self.typeM,
                'rows_compl': get_rows_compl(me, self.typeM, self.typeE)
                } 
        # Reformulation of ME to NLP
        self.A = None
        self.B = None
        self.S = None
        self.b = None
        self.c = None
        self.xl= None
        self.xu= None
        # Inputs to qminos
        self.J = None
        self.nnCon = None
        self.nnJac = None
        self.neJac = None
        self.ne    = None
        self.ha    = None
        self.ka    = None
        self.ad    = None
        self.bld   = None
        self.bud   = None
        self.mu0   = None
        self.probname = "me_nlp"
        self.M        = None
        self.N        = None
        self.nb       = None
        # Solution and exit flag
        self.x     = None
        self.inform = np.array(0)
        # Hold LP results and options
        self.lp_inform  = None
        self.lp_hs      = None
        self.lp_x       = None
        # Initialize solver options
        self.init_solver_opts()

    def set_realopts(self, model, param_dict):
        """
        Set real valued solver options
        model:  'lp', 'nlp'
        param_dict: {param: val}
                    param:  name of param value
                    val:    real-valued value
        """
        for param, val in six.iteritems(param_dict):
            if param in self.opt_realdict[model]:
                if isinstance(val, float):
                    self.opt_realdict[model][param]=val
                else:
                    warnings.warn('val must be float')
            else:
                warnings.warn('Unrecognized solver option:'+param)

    def set_intopts(self, model, param_dict):
        """
        Set int-valued solver options
        model:  'lp', 'nlp'
        param_dict: {param: val}
                    param:  name of param value
                    val:    int-valued value
        """
        for param, val in six.iteritems(param_dict):
            if param in self.opt_intdict[model]:
                if isinstance(val, int):
                    self.opt_intdict[model][param] = val
                else:
                    warnings.warn('val must be int')
            else:
                warnings.warn('Unrecognized solver option:'+param)


    def init_solver_opts(self):
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
        # NLP solver options
        #----------------------------------------------------
        # Width of characters allowed in each options for qMINOS
        self.opt_strwidth['nlp'] = 72
        self.opt_realwidth['nlp'] = 55
        self.opt_intwidth['nlp'] = 55
        self.opt_strlist['nlp'] = [
                'Maximize',     # Default obj sense is to maximize
                'Completion full',
                'Print level (jflxb) 00001',
                'Solution No'
                ]
        self.opt_intdict['nlp'] = {
                'Major iterations': 1000,
                'Superbasics limit': 40,
                'Verify level': 0,
                'Scale option': 2,
                'Partial price': 1,
                'Iterations': 10000,
                'Print frequency': 100000,
                'Summary level': 0,
                'Summary frequency': 100,
                'Solution file': 9,
                'New basis file': 11,
                'Save frequency': 500000
                }
        self.opt_realdict['nlp'] = {
                'Penalty parameter':100.0,
                'LU factor tol': 1.1,
                'LU update tol': 1.1,
                'LU singularity tol': 1e-30,
                'Feasibility tol': 1e-15,
                'Optimality tol': 1e-15,
                'Unbounded step size': 1e+30
                }

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

    def make_matrices(self, max_mu = True):
        """Constructs NLP matrices from me model object"""
        A,B,d,S,b,c,xl,xu,csense,csense_nonlin = me2nlp(self.me,
                growth_symbol=self.growth_symbol,
                scaleUnits=self.scaleUnits,
                unitDict=self.unitDict,
                max_mu = max_mu, growth_rxn=self.growth_rxn)
        self.A = A
        self.B = B
        self.d = d
        self.S = S
        self.b = b
        self.c = c
        self.xl = xl
        self.xu = xu
        self.csense = csense
        self.csense_nonlin = csense_nonlin


    def make_lp_for_nlp(self, mu_fix, verbosity=0):
        """
        Construct LP whose basis is compatible with ME-NLP
        """
        from cobrame import mu

        me = self.me

        if self.A is None:
            self.make_matrices()
        # Nonlinear constraints
        # Substituting mu is as simple as mu*A
        A = self.A*mu_fix
        B = self.B
        # Linear constraints
        S = self.S
        b = self.b
        c = [r.objective_coefficient for r in me.reactions]
        # self.xl doesn't account for symbolic mu, so just build it anew here...
        xl = np.matrix([r.lower_bound for r in me.reactions]).transpose()
        xu = np.matrix([r.upper_bound for r in me.reactions]).transpose()
        # ... and substitute mu in bounds
        for j,rxn in enumerate(me.reactions):
            lb = rxn.lower_bound
            ub = rxn.upper_bound
            if hasattr(lb, 'subs'):
                xl[j] = float(lb.subs(mu, mu_fix))
            if hasattr(ub, 'subs'):
                xu[j] = float(ub.subs(mu, mu_fix))

        # This J has extra row added. Also, bl & bu have extra slack (unbounded) for
        # the "objective" row
        J, ne, P, I, V, bl, bu = makeME_LP_for_NLP(A,B,S,b,c,xl,xu) 

        # Solve a single LP
        m,n = J.shape
        ha = I
        ka = P
        ad = V
        bld = [bi for bi in bl.flat]
        bud = [bi for bi in bu.flat]
        nb = m + n
        hs = np.zeros(nb, np.dtype('i4'))
        return m, n, ha, ka, ad, bld, bud, hs


    def make_lp(self, mu_fix, verbosity=0):
        """
        Construct LP problem for qMINOS or MINOS. 
        """
        from cobrame import mu

        me = self.me
        S = me.construct_S(mu_fix).tocsc()
        xl = np.matrix([r.lower_bound for r in me.reactions]).transpose()
        xu = np.matrix([r.upper_bound for r in me.reactions]).transpose()
        # Also substitute mu in bounds
        for j,rxn in enumerate(me.reactions):
            lb = rxn.lower_bound
            ub = rxn.upper_bound
            if hasattr(lb, 'subs'):
                xl[j] = float(lb.subs(mu, mu_fix))
            if hasattr(ub, 'subs'):
                xu[j] = float(ub.subs(mu, mu_fix))

        #b = [0. for m in me.metabolites]
        b = [m._bound for m in me.metabolites]
        c = [r.objective_coefficient for r in me.reactions]
        # constraint sense eventually be in the metabolite...
        # csense = ['E' for m in me.metabolites]
        csense = [m._constraint_sense for m in me.metabolites]
        J, ne, P, I, V, bl, bu = makeME_LP(S,b,c,xl,xu,csense) 

        # Solve a single LP
        m,n = J.shape
        ha = I
        ka = P
        ad = V
        bld = [bi for bi in bl.flat]
        bud = [bi for bi in bu.flat]
        nb = m + n
        hs = np.zeros(nb, np.dtype('i4'))
        return m, n, ha, ka, ad, bld, bud, hs


    def solvelp(self, muf, quad=True, basis=None, nlp_compat=False, verbose=False,
            lpopt_file = 'fort.14', verbosity=0, precision='quad'):
        """
        x, status, hs = solvelp(self, muf, quad=True, basis=None, nlp_compat=False, verbose=False) 

        Solve LP at mu using qMINOS (quad=True) or MINOS (quad=False).
        Pass the basis (hs) back and forth with Fortran for warm-start.
        Inputs:
        muf: fixed growth rate
        basis: basis vector
        nlp_compat  If True, returns basis that is compatible with ME-NLP

        Outputs:
        x: primal solution
        status: solver status
        hs: basis

        me_nlp.pi (me.solution.y or me.solution.y_dict): dual values (shadow prices)
        me_nlp.rc: reduced costs
        """
        me = self.me

        hs = basis
        if nlp_compat:
            m,n,ha,ka,ad,bld,bud, hs0 = self.make_lp_for_nlp(muf, verbosity=verbosity)
        else:
            m,n,ha,ka,ad,bld,bud, hs0 = self.make_lp(muf, verbosity=verbosity)

        if hs is None:
            warm = False
            hs = hs0
        else:
            warm = True
            if verbose:
                print('Using provided basis of length %d' % len(hs))

        import os.path

        inform = np.array(0)
        probname = 'me_lp'

        if verbosity > 0:
            print('Getting MINOS parameters from ME_NLP...')

        precision = precision.lower()

        tic = time.time()
        if precision == 'quad':
            stropts,intopts,realopts,intvals,realvals,nStrOpts,nIntOpts,nRealOpts =\
                self.get_solver_opts('lp')
            x,pi,rc = qwarmLP.qwarmlp(inform, probname, m, ha, ka, ad, bld, bud, hs, warm,
                                      stropts, intopts, realopts, intvals, realvals,
                                      nstropts = nStrOpts,
                                      nintopts = nIntOpts,
                                      nrealopts = nRealOpts)
        elif precision == 'double':
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
            if precision == 'dqq':
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
        self.lp_hs = hs
        self.x = x
        # Save dual and reduced cost information
        self.pi = pi
        # Reduced cost: g - (A I)'*pi, where g is the gradient, l <= A*x <= u are constraints
        # including the objective function in the last row
        self.rc = rc

        # Write the solution to the ME model's solution for a consistent solve interface
        #f = x[0]
        # Aug 27, 2015: obj coeffs are not always mu (or x[0])
        f =sum([rxn.objective_coefficient * x[j] for j,rxn in enumerate(self.me.reactions)])
        x_primal = x[0:len(self.me.reactions)]   # The remainder are the slacks
        x_dict = {rxn.id: x[j] for j,rxn in enumerate(self.me.reactions)}
        y = pi
        # J = [S; c]
        y_dict = {met.id: y[i] for i,met in enumerate(me.metabolites)}
        y_dict['linear_objective'] = y[len(y)-1]
        status = self.inform
        if int(status) == 0:
            status = 'optimal'
        self.me.solution = Solution(f, x_primal, x_dict, y, y_dict, 'qminos', time_elapsed, status)

        return x, status, hs


    def bisectmu(self, precision=1e-3, mumin=0.0, mumax=2.0, maxIter=100, quad=True,
                 basis=None, nlp_compat=False, check_feas0=False, zero_mu=1e-3,
                 verbosity=2):
        """
        muopt, hs, xopt, cache = bisectmu(precision=1e-3, mumin=0.0, mumax=2.0)

        Bisection to maximize mu using qMINOS. Sequence of feasibility problems.
        More efficient than the more general binary search, which finds the extremum
        of a unimodal function.
        TODO: set quad=False if need fast basis for subsequent quadMINOS
        """
        import copy as cp

        if self.nb is None:
            self.make_nlp()

        me = self.me
        hs = basis

        # basis (hs) intent(inout). Will generate first basis from Cold-start if hs=None
        # Save solutions, recovered at golden ratios
        cache = {}
        # Store the final Solution object
        solution = None

        # Check feasibility at mu=zero?
        if check_feas0:
            x0, stat0, hs0 = self.solvelp(zero_mu, nlp_compat=nlp_compat, verbosity=verbosity,
                    basis=hs)
            if me.solution.status is not 'optimal':
                warnings.warn('Infeasible at mu=%g. Returning.'%zero_mu)
                return zero_mu, hs0, x0, cache 
            else:
                hs = hs0

        def checkmu(muf, hs):
            if muf not in cache:
                x_new, stat_new, hs_new = self.solvelp(
                    muf, basis=hs, nlp_compat=nlp_compat, verbosity=verbosity)
                if me.solution.status is 'optimal':
                    hs = hs_new
                stat = me.solution.status
                sol = cp.deepcopy(me.solution)
                cache[muf] = stat

            return cache[muf], hs, sol, x_new

        warm = False
        converged = False
        a = mumin
        b = mumax
        muopt = a
        xopt = None
        iter = 0

        if verbosity >= 2:
            print('iter\tmuopt    \ta     \tb     \tmu1       \tstat')
        while iter < maxIter and not converged:
            # Just a sequence of feasibility checks
            mu1 = (a+b)/2.
            # Retrieve evaluation from cache if it exists: golden section advantage
            stat1, hs, sol1, x1 = checkmu(mu1, hs)
            if stat1 is 'optimal':
                a = mu1
                muopt = mu1
                solution = sol1
                xopt = x1
            else:
                b = mu1

            converged = abs(b-a) <= precision
            warm = True
            iter = iter+1

            if verbosity >= 2:
                print(iter, muopt, a, b, mu1, stat1)

        # Save final solution
        me.solution = solution

        return muopt, hs, xopt, cache 


    def bisectme(self, precision=1e-3, mumin=0.0, mumax=2.0, maxIter=100, quad=True, golden=True, basis=None, nlp_compat=False, check_feas0=False, zero_mu=1e-3, verbosity=0):
        """
        [Deprecated] legacy bisection using qMINOS. 
        TODO: set quad=False if need fast basis for subsequent quadMINOS
        """
        import copy as cp

        if self.nb is None:
            self.make_nlp(verbosity=verbosity)

        me = self.me
        hs = basis

        # basis (hs) intent(inout). Will generate first basis from Cold-start if hs=None
        # Save solutions, recovered at golden ratios
        cache = {}
        # Store the final Solution object
        solution = None

        # Check feasibility at mu=zero?
        if check_feas0:
            x0, stat0, hs0 = self.solvelp(zero_mu, nlp_compat=nlp_compat)
            if me.solution.status is not 'optimal':
                warnings.warn('Infeasible at mu=%g. Returning.'%zero_mu)
                return zero_mu, hs0, x0, cache 
            else:
                hs = hs0

        def checkmu(muf, hs):
            if muf not in cache:
                x_new, stat_new, hs_new = self.solvelp(
                    muf, basis=hs, nlp_compat=nlp_compat, verbosity=verbosity)
                if me.solution.status is 'optimal':
                    hs = hs_new
                stat = me.solution.status
                sol = cp.deepcopy(me.solution)
                cache[muf] = stat

            return cache[muf], hs, sol, x_new


        phi = 0.75  # If using binary, instead of golden-section, search
        if golden:
            phi = (-1 + 5.**0.5)/2.
        warm = False
        converged = False
        a = mumin
        b = mumax
        muopt = a
        xopt = None
        iter = 0

        if verbosity > 1:
            print('iter\tmuopt    \ta     \tb     \tmu1       \tmu2       \tstat1  \tstat2')
        while iter < maxIter and not converged:
            mu1 = b - (b - a) * phi
            mu2 = a + (b - a) * phi
            # Retrieve evaluation from cache if it exists: golden section advantage
            stat1, hs, sol1, x1 = checkmu(mu1, hs)
            if stat1 is not 'optimal':
                # Implies stat2 is also infeasible
                #b = mu2
                b = mu1
                stat2 = 'infeasible'
            else:
                stat2, hs, sol2, x2 = checkmu(mu2, hs)
                if stat2 is not 'optimal':
                    a = mu1  # a is feasible
                    b = mu2
                    muopt = mu1
                    solution = sol1
                    xopt = x1
                else:        # both feasible
                    a = mu2 
                    muopt = mu2
                    solution = sol2
                    xopt = x2

            converged = abs(b-a) <= precision
            warm = True
            iter = iter+1

            if verbosity>1:
                print(iter, muopt, a, b, mu1, mu2, stat1, stat2)

        # Save final solution
        me.solution = solution

        return muopt, hs, xopt, cache 


    def make_nlp(self, verbosity=0):
        """
        Constructs NLP problem for qMINOS. Requires matrices from
        make_matrices. Thus, will run make_matrices if the
        required matrices do not exist yet.
        """
        if self.A is None:
            self.make_matrices()

        J,nnCon,nnJac,neJac,ne,P,I,V,bl,bu = makeME_NLP(self.A, self.B, 
                self.S, self.b, self.c, self.xl, self.xu)

        M,N = J.shape

        self.M = M
        self.N = N
        self.nnCon = nnCon
        self.nnJac = nnJac
        self.neJac = neJac
        self.nb = M+N
        self.ne = ne
        self.ha = I
        self.ka = [int(pi) for pi in P]
        self.ad = V
        self.bld = [bi for bi in bl.flat]
        self.bud = [bi for bi in bu.flat]

    def construct_S(self, growth_rate):
        """
        From cobrame--in case me does not have construct_S
        """
        me = self.me
        growth_key = self.growth_key

        # intialize to 0
        S = sps.dok_matrix((len(me.metabolites), len(me.reactions)))
        # populate with stoichiometry
        for i, r in enumerate(me.reactions):
            for met, value in six.iteritems(r._metabolites):
                met_index = me.metabolites.index(met)
                if hasattr(value, "subs"):
                    S[met_index, i] = float(value.subs(growth_key, growth_rate))
                else:
                    S[met_index, i] = float(value)
        return S

    def varyme(self, mu_fixed, rxns_fva0, basis=None, verbosity=0):
        """
        fva_result, fva_stats = varyme(self, mu_fixed, rxns_fva)

        rxns_fva:  list of reactions to be varied (Reaction objects or ID)

        High-level interface for qvaryME (quad-prec FVA)
        12 Aug 2015: first version. Must fix bugs.
        """
        from qminospy import qvaryME
        from cobrame import mu
        import time as time
        import six

        me = self.me
        hs = basis

        if isinstance(rxns_fva0[0], six.string_types):
            rxns_fva = [me.reactions.get_by_id(rid) for rid in rxns_fva0]
        else:
            rxns_fva = rxns_fva0

        if hasattr(me, 'construct_S'):
            S = me.construct_S(mu_fixed).tocsc()
        else:
            S = self.construct_S(mu_fixed).tocsc()

        xl = np.matrix([r.lower_bound for r in me.reactions]).transpose()
        xu = np.matrix([r.upper_bound for r in me.reactions]).transpose()
        # Also substitute mu in bounds
        for j,rxn in enumerate(me.reactions):
            lb = rxn.lower_bound
            ub = rxn.upper_bound
            if hasattr(lb, 'subs'):
                xl[j] = float(lb.subs(mu, mu_fixed))
            if hasattr(ub, 'subs'):
                xu[j] = float(ub.subs(mu, mu_fixed))

        b = [m._bound for m in me.metabolites]
        c = [r.objective_coefficient for r in me.reactions]

        obj_inds0 = [me.reactions.index(rxn) for rxn in rxns_fva for j in range(0, 2)]
        obj_coeffs = [ci for rxn in rxns_fva for ci in (1.0, -1.0)]
        # csense = ['E' for m in me.metabolites]
        csense = [m._constraint_sense for m in me.metabolites]

        J,ne,P,I,V,bl,bu, obj_inds = makeME_VA(S, b, xl, xu, csense, obj_inds0, obj_coeffs)

        m,n = J.shape
        ha = I
        ka = P
        ad = V
        bld = [bi for bi in bl.flat]
        bud = [bi for bi in bu.flat]
        nb = m + n
        if hs is None:
            warm = False
            hs = np.zeros(nb, np.dtype('i4'))
        else:
            warm = True
            if verbosity > 0:
                ('Warm-starting first run using basis of length %d' % len(hs))

        # Get MINOS options
        if verbosity > 0:
            print('Getting MINOS parameters for LP')
        stropts,intopts,realopts,intvals,realvals,nStrOpts,nIntOpts,nRealOpts =\
            self.get_solver_opts('lp')

        nVary = len(obj_inds)
        obj_vals = np.zeros(nVary)
        fva_stats = np.zeros(nVary, np.dtype('i4'))
        probname = 'varyme'

        tic = time.time()
#         qvaryME.qvaryme(fva_stats, probname, m, ha, ka, ad, bld, bud,
#                 obj_inds, obj_coeffs, obj_vals)
        qvaryME.qvaryme(fva_stats, probname, m, ha, ka, ad, bld, bud, hs, warm,
                obj_inds, obj_coeffs, obj_vals,
                stropts, intopts, realopts, intvals, realvals,
                nstropts = nStrOpts,
                nintopts = nIntOpts,
                nrealopts = nRealOpts)

        t_elapsed = time.time()-tic
        if verbosity>0:
            print('Finished varyME in %f seconds for %d rxns (%d quadLPs)' %
                  (t_elapsed, len(rxns_fva), len(obj_inds)))

        # Return result consistent with cobrame fva
        fva_result = {
            (self.me.reactions[obj_inds0[2*i]].id):{
                'maximum':obj_vals[2*i], 
                'minimum':obj_vals[2*i+1] } for i in range(0, nVary/2) }

        # Save updated basis
        self.hs = hs
        self.lp_hs = hs

        return fva_result, fva_stats

    def reorder_rxns(self, rxn_list, where='first'):
        """
        rxnids_old = reorder_rxns(self, rxns, where='first') 

        Inputs:
        rxns        [cobra.Reaction] or [str] 
        where       push given rxns 'first' or 'last'

        Values:
        rxns_old      Reactions in original order

        Reorder rxns by placing given rxns to 'where'
        E.g., if where=first, push given rxns to top of DictList
        """
        me = self.me

        if isinstance(rxn_list[0], six.string_types):
            rxns = [me.reactions.get_by_id(r) for r in rxn_list]
        else:
            rxns = rxn_list

        rxns0 = DictList([r for r in me.reactions])        # Original order
        #nRxn = len(me.reactions)
        #order_new = range(0,nRxn)
        #for i,rxn in enumerate(rxns):
        #    ind_old = me.reactions.index(rxn)
        #    if where is 'first':
        #        i_new = i
        #    elif where is 'last':
        #        i_new = nRxn-i-1
        #    else:
        #        raise ValueError("where must be 'first' or 'last'. Was given:"+str(where))
        #    order_new[ind_old] = i_new
        #    order_new[i_new] = ind_old
        #me.reactions = DictList( [ me.reactions[i] for i in order_new] )

        # Use built-in sorted method instead
        rxns_new_order = sorted(me.reactions, key=lambda r: 0 if r in rxns else 1)
        me.reactions = DictList(rxns_new_order)

        return rxns0

    def solve(self, x0=None, basis=None, param_file='fort.4'):
        """
        High-level interface for solving a ME model. Cleaner
        interface since the nonlinear program matrices do not
        need to be passed around.

        30 Jul 2015: first version
        14 Aug 2015: can warm-start by providing basis (hs)
        """
        from qminospy import qsolveME as qminos
        import time as clock
        import os.path

        x = self.x
        hs = basis

        time_elapsed = 0

        if self.mu0 is None and x0 is None:
            raise ValueError('Must set self.mu0 first!')
        else:
            if hs is None:
                warm = False
                nb = self.nb
                hs = np.zeros(nb, np.dtype('i4'))
            else:
                warm = True
                print('Using provided basis of length %d' % len(hs))
            if x0 is None:
                x0 = np.zeros(nb)
                x0[0] = self.mu0


            ### Set solver params
            import os.path
            # if param_file is not None and os.path.isfile(param_file):
            #     print('Getting MINOS parameters from specs file:', param_file)
            # else:

            print('Getting MINOS parameters from ME_NLP...')
            stropts,intopts,realopts,intvals,realvals,nStrOpts,nIntOpts,nRealOpts =\
                self.get_solver_opts('nlp')

            tic = clock.time()

            self.inform = np.array(0)
            qminos.qsolveme(x0, self.inform, self.mu0, self.probname, self.M,
                    self.nnCon, self.nnJac, self.neJac, self.ha, self.ka, self.ad,
                    self.bld, self.bud, hs, warm,
                    stropts, intopts, realopts, intvals, realvals,
                    nstropts = nStrOpts,
                    nintopts = nIntOpts,
                    nrealopts = nRealOpts)

            # x0 is intent(inout), so qsolveme has updated x0 with x
            x = x0

            time_elapsed = clock.time() - tic
            # Update solution in self
            self.x = x
            print("Time elapsed: ", time_elapsed, "seconds")
            print("Status: ", self.inform)

        # Write the solution to the ME model's solution for a consistent solve interface
        f = x[0]
        x_primal = x[0:len(self.me.reactions)]   # The remainder are the slacks
        x_dict = {rxn.id: x[j] for j,rxn in enumerate(self.me.reactions)}
        y = None
        y_dict = None
        status = self.inform
        if int(status) == 0:
            status = 'optimal'
        self.me.solution = Solution(f, x_primal, x_dict, y, y_dict, 'qminos', time_elapsed, status)

        # Return
        return x, self.inform, hs

    def update_obj(self):
        """
        Update matrices and NLP based on obj_coeffs in self.me
        from ME model
        """
        for j,rxn in enumerate(self.me.reactions):
            self.c[j] = rxn.objective_coefficient
        # Remaking nlp is fast now
        self.make_nlp()

    def update_bounds(self):
        """
        Update matrices and NLP based on bounds in self.me
        from ME model
        """
        # J,nnCon,nnJac,neJac,ne,P,I,V,bl,bu = makeME_NLP(self.A, self.B, 
        #         self.S, self.b, self.c, self.xl, self.xu)
        for j,rxn in enumerate(self.me.reactions):
            # If mu in bounds, warn and set to unbounded
            lb = rxn.lower_bound 
            ub = rxn.upper_bound 
            if hasattr(lb, 'subs'):
                warnings.warn('lb for %s is mu-dependent. Setting to 0.0'%(rxn.id))
                lb = 0.0
            if hasattr(ub, 'subs'):
                warnings.warn('ub for %s is mu-dependent. Setting to 1000.0'%(rxn.id))
                ub = 1000.0

            self.xl[j] = lb
            self.xu[j] = ub

        # Remaking nlp is fast now
        self.make_nlp()

    def solvenlp(self, precision=0.01, max_iter=20, check_feas0=False, zero_mu=1e-3, basis=None,
                 auto_update_bounds=True, auto_update_obj=True, verbosity=0):
        """
        x, stat, hs = solvenlp(self, precision=0.01, max_iter=20, check_feas0=False,
                               zero_mu=1e-3, basis=None) 

        The combined solution procedure: bisection (golden section) to
        get initial mu0, followed by qsolveME to find mu*
        15 Aug 2015: first version
        """
        if self.nb is None:
            self.make_nlp()

        hs = basis
        # Check feasibility at mu0 = zero_mu?
        if check_feas0:
            x0, stat0, hs0 = self.solvelp(zero_mu, nlp_compat=True, basis=None)
            if stat0 is not 'optimal':
                #raise ValueError('Infeasible at mu=0.0. Stopping.')
                warnings.warn('Infeasible at mu=%g. Returning.'%zero_mu)
                return x0, stat0, hs0
            else:
                hs = hs0

        # Bisection (golden section)
        tic1 = time.time()
        mu_bs, hs_bs, x_bs, cache_bs = self.bisectmu(precision=precision,
                maxIter=max_iter, nlp_compat=True, basis=hs)
        time_bs = time.time()-tic1

        # NLP
        if hs_bs is None or x_bs is None:
            warnings.warn('Feasible mu0 not found with bisectME. Returning.')
            return x_bs, 'infeasible', hs_bs
        else:
            if auto_update_bounds:
                if verbosity>1:
                    print('Updating bounds to me')
                self.update_bounds()
            if auto_update_obj:
                if verbosity>1:
                    print('Updating objective to me')
                self.update_obj()

            tic2 = time.time()
            self.mu0 = mu_bs
            x, stat, hs = self.solve(x0=x_bs[0:self.nb], basis=hs_bs[0:self.nb])
            time_nlp = time.time()-tic2

            t_elapsed = time.time()-tic1

            if verbosity>0:
                print('Finished in %f seconds (%f bisectME, %f ME-NLP)' %
                      (t_elapsed, time_bs, time_nlp))
            # Return the basis from the LP, since that is what will be used to
            # warm-start solvenlp. We could return the NLP basis, too.

            return x, stat, hs_bs

    def solve_qclp(self):
        """
        Solve quadratically constrained LP (linear objective)
        01 Dec 2015
        """
        pass


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

    def update_dilution_fluxes(self, csense='G'):
        """
        Update dilution fluxes, especially keffs that were updated by
        setting a reaction.keff to a new value.
        """
        return self.make_dilution_fluxes(csense)


    def make_dilution_fluxes(self, csense='G', macromolecules=['complex'], verbosity=1,
            LB_DIL=0.0, UB_DIL=1e6, mols_exclude=[]):
        """
        rxns_dil, constraints_dil, macromol_error = make_dilution_fluxes(self, csense='G', macromolecules=['complex'])

        Allow dilution to be >= sum_j mu/keffj v_usagej (i.e., that accounted for by
        sum of complex usage fluxes). Useful when forcing macromolecule synthesis
        beyond usable capacity.
        
        Simply add EXTRA_DILUTION flux for every diluted macromolecule. 
        Works because we would only allow >= constraints on dilution flux anyway.

        Inputs:
        csense          constraint sense on dilution fluxes (if G: v_dili >= mu*sum_j/keffij*v_usej
        macromolecules  list of macromolecule classes that will get dilution fluxes.
                        One or more of 'complex','protein','rna','all' (default: ['complex'])
        mols_exclude    macromolecules to exclude from adding extra dilution fluxes

        Make separate dilution fluxes for each complex, with specified
        constraint-sense.
        Especially useful for cobrame (ME 2.0) models, where dilution fluxes
        are often constrained by equalities, and we sometimes want to work
        with inequalities (e.g., sampling keffs over a wider range without
        making model infeasible).
        """
        me = self.me
        rxns_dil = []           # return dilution rxns
        constraints_dil = []    # return dilution coupling constraints
        macromol_error = []     # report macromolecules where failed to add dilution

        if isinstance(me, MEmodel):
            """
            For now, should only support this function for cobrame models
            since ME 1.0 models already have separate dilution fluxes
            """
            #================================================
            # 1. Get all the diluted macromolecules and add dilution fluxes
            #================================================
            macromols = []
            if 'complex' in macromolecules:
                macromols = macromols + [d.complex if isinstance(d.complex, six.string_types)
                        else d.complex.id for d in me.complex_data]
            if 'protein' in macromolecules:
                macromols = macromols + [d.protein if isinstance(d.protein, six.string_types)
                        else d.protein.id for d in me.translation_data]
            if 'rna' in macromolecules:
                macromols = macromols + list(set(
                    [r for d in me.transcription_data for r in d.RNA_products]))

            # Exclude specified macromolecules
            if len(mols_exclude)>0:
                macromols = [m for m in macromols if m not in mols_exclude]

            for mol_id in macromols:
                try:
                    mol = me.metabolites.get_by_id(mol_id)
                except KeyError:
                    macromol_error.append(mol_id)
                    if verbosity > 1:
                        warnings.warn('Unable to retrieve macromolecule: %s'%mol_id)

                rxn_dil_id = 'extra_dilution_'+mol.id
                # Usually, dilution rxn will not be in model, so try adding first
                try:
                    rxn_dil = Reaction(rxn_dil_id)
                    me.add_reaction(rxn_dil)
                except:
                    rxn_dil = me.reactions.get_by_id(rxn_dil_id)

                #============================================
                # 2. Add macromolecule sink for extra dilution
                #    beyond concentration accounted for by complex usage flux
                #
                #    macromolecule: v_synth - v_dil - v_dil_extra = 0
                #============================================
                rxn_dil.add_metabolites({mol: -1.}, combine=False)
                rxn_dil.lower_bound = LB_DIL
                rxn_dil.upper_bound = UB_DIL
                rxns_dil.append(rxn_dil)

        else:
            warnings.warn('Dilution fluxes should already be inequalities for ME 1.0 model. Quitting.')

        if verbosity > 0:
            print('Done adding dilution fluxes.')
            if len(macromol_error)>0:
                print('Failed to retrieve %d macromolecules from model' %
                      len(macromol_error))

        return rxns_dil, constraints_dil, macromol_error

## End of ME_NLP
##//============================================================


def writeME_NLP(me, outname=None):
    """
    Builds ME NLP data matrix from a ME 2.0 model, to be passed to
    qsolveME
    [LY] 21 Jul 2015: first version
    """
    A,B,d,S,b,c,xl,xu,csense,csense_nonlin = me2nlp(me)
    J, nnCon, nnJac, neJac, nzS, P, I, V, bl, bu  = makeME_NLP(A,B,S,b,c,xl,xu)
    if outname is not None:
        print('Writing to file: ', outname)
        dumpMat(J, outname, nnCon, nnJac, neJac, nzS, P, I, V, bl, bu)
        print('Done!')

    return J, nnCon, nnJac, neJac, bl, bu


def me2nlp(me, growth_symbol='mu', scaleUnits=False, LB=0.0, UB=1000.0,
           growth_rxn='biomass_dilution', unitDict={'e_mult': 1e-6,
              'typeE': [cobrame.TranscriptionReaction,
                        cobrame.TranslationReaction]},
              max_mu=True):
    """
    From ME model object, create NLP data matrices of the form:
    max mu = c'*x
    mu,x
    s.t.  mu*A*x + B*x = 0
          S*x          = b
          xl <= x <= xu

    21 Jul 2015: [LY] first version. Must reorder rxns such that mu (biomass dilution) is x[0]
    """
    import numpy as np
    from cobra.core.DictList import DictList

    # Reorder rxns such that growth_rxn is the first rxn
    nRxn = len(me.reactions)
    rxn_mu = me.reactions.get_by_id(growth_rxn)
    ind_mu = me.reactions.index(rxn_mu)
    order_new = list(range(0, nRxn))
    order_new[ind_mu] = 0
    order_new[0] = ind_mu
    me.reactions = DictList([me.reactions[i] for i in order_new])

    # Identify mets in nonlinear vs. linear constraints
#     mets_nonlin = [met for met in me.metabolites if sum(
#         [growth_symbol in str(rxn.metabolites[met]) for rxn in met.reactions ])
#     ]
#     mets_linear = list(set(me.metabolites).difference(set(mets_nonlin)))
#     mets_nonlin = DictList(mets_nonlin)
#     mets_linear = DictList(mets_linear)

    mets_linear, mets_nonlin = get_lin_nonlin_mets(me, growth_symbol)

    # Create A, B: nonlinear constraints
    A,B,d,csense_nonlin = make_nonlin_constraints(me, mets_nonlin, growth_symbol, scaleUnits, unitDict)
    # Create S & b: linear constraints
    S,b,csense = make_linear_constraints(me, mets_linear, scaleUnits, unitDict)
    # Create c: objective coefficient vector
    # c = [r.objective_coefficient for r in me.reactions]
    # Linear objective
    if max_mu:
        c = [0.0 for r in me.reactions]
        c[me.reactions.index(rxn_mu)] = 1.0
    else:
        c = [r.objective_coefficient for r in me.reactions]

    # Create xl and xu: primal bounds
    # bounds might have mu (growth_symbol) in it, which should be corrected
    xl = np.matrix([r.lower_bound for r in me.reactions]).transpose()
    xu = np.matrix([r.upper_bound for r in me.reactions]).transpose()

#   Sept 04, 2015: skip because comparing mu in numpy arrays seems to 
#   cause "cannot compare the truth value of ..." error in some installations
#     print('xl: ', xl.min(), xl.max())
#     print('xu: ', xu.min(), xu.max())

    # Replace mu with default bounds
    from sympy.core.symbol import Symbol

    for i in range(0, len(xl)):
        if isinstance(xl.item(i), Symbol):
            xl[i] = LB
        if isinstance(xu.item(i), Symbol):
            xu[i] = UB

    return A, B, d, S, b, c, xl, xu, csense, csense_nonlin


def me2nlp_general(me, growth_symbol='mu', scaleUnits=False, LB=0.0, UB=1000.0,
                   growth_rxn='biomass_dilution',
                   unitDict={'e_mult': 1e-6,
                             'typeE': [
                                 cobrame.TranscriptionReaction,
                                 cobrame.TranslationReaction]},
                   max_mu=True):
    """
    From ME model object, create NLP data matrices of the form:
    max mu = c'*x
    mu,x
    s.t.  gk(mu)*A*x + B*x = 0
          S*x          = b
          xl <= x <= xu

    where gk(mu) are a family of functions
    """
    import numpy as np
    from cobra.core.DictList import DictList

    # Reorder rxns such that growth_rxn is the first rxn
    nRxn = len(me.reactions)
    rxn_mu = me.reactions.get_by_id(growth_rxn)
    ind_mu = me.reactions.index(rxn_mu)
    order_new = list(range(0, nRxn))
    order_new[ind_mu] = 0
    order_new[0] = ind_mu
    me.reactions = DictList([me.reactions[i] for i in order_new])

    mets_linear, mets_nonlin = get_lin_nonlin_mets(me, growth_symbol)

    # Create A, B: nonlinear constraints
    A,B,d,csense_nonlin = make_nonlin_constraints(me, mets_nonlin, growth_symbol, scaleUnits, unitDict)
    # Create S & b: linear constraints
    S,b,csense = make_linear_constraints(me, mets_linear, scaleUnits, unitDict)
    # Create c: objective coefficient vector
    # Linear objective
    if max_mu:
        c = [0.0 for r in me.reactions]
        c[me.reactions.index(rxn_mu)] = 1.0
    else:
        c = [r.objective_coefficient for r in me.reactions]

    # Create xl and xu: primal bounds
    # bounds might have mu (growth_symbol) in it, which should be corrected
    xl = np.matrix([r.lower_bound for r in me.reactions]).transpose()
    xu = np.matrix([r.upper_bound for r in me.reactions]).transpose()

    # Replace mu with default bounds
    from sympy.core.symbol import Symbol

    for i in range(0,len(xl)):
        if isinstance(xl.item(i), Symbol):
            xl[i] = LB
        if isinstance(xu.item(i), Symbol):
            xu[i] = UB

    return A,B,d,S,b,c,xl,xu,csense,csense_nonlin


def make_nonlin_constraints(me, mets_nonlin, growth_symbol='mu', scaleUnits=False, unitDict=None,
                            verbose=False):
    """
    A,B,d,csense_nonlin = make_nonlin_constraints(me, mets_nonlin, growth_symbol, scaleUnits, unitDict)

    17 Nov 2015: [LY] reordered rxns so that nonlin variables also appear first
    08 Mar 2016: accept polynomials of mu of any order
    """
    # Return sparse matrices, A and B, where mu*A*x + B*x = 0
    import scipy.sparse as sps
    Am = len(mets_nonlin)
    Bm = Am
    N = len(me.reactions)
    # Find columns involved in nonlinear terms
#     rxns_nonlin = set (
#             [rxn for sublist in ([met.reactions for met in mets_nonlin ]) for rxn in sublist ]
#             ) 
    # Make constraints
    cons_nl = get_nonlin_triplets(me, mets_nonlin, growth_symbol, scaleUnits, unitDict)
    dataA = [a for m,r,a,b in cons_nl if a != 0]
    dataB = [b for m,r,a,b in cons_nl if b != 0]
    # Re-order rxns so that nonlin variables appear first (vars in A)
    rxnsnl = list(set([r for m,r,a,b in cons_nl if a != 0]))
    if verbose:
        print(len(rxnsnl), 'variables in nonlinear constraints')
    #print('Reordering rxns so nonlinear variables come first')
    #me.reactions = DictList( rxnsnl + [r for r in me.reactions if r not in rxnsnl] )

    row_indA = [mets_nonlin.index(m) for m,r,a,b in cons_nl if a != 0]
    col_indA = [me.reactions.index(r) for m,r,a,b in cons_nl if a != 0]
    row_indB = [mets_nonlin.index(m) for m,r,a,b in cons_nl if b != 0]
    col_indB = [me.reactions.index(r) for m,r,a,b in cons_nl if b != 0]

    A = sps.csc_matrix( (dataA, (row_indA, col_indA)), shape=(Am,N) )
    B = sps.csc_matrix( (dataB, (row_indB, col_indB)), shape=(Bm,N) )
    d = [met._bound for met in mets_nonlin]
    csense_nonlin = [met._constraint_sense for met in mets_nonlin]

    return A,B,d,csense_nonlin


def make_nonlin_constraints_general(me, mets_nonlin, growth_symbol='mu', scaleUnits=False, unitDict=None):
    """
    A,B,d,csense_nonlin,g = make_nonlin_constraints(me, mets_nonlin, growth_symbol, scaleUnits, unitDict)

    17 Nov 2015: [LY] reordered rxns so that nonlin variables also appear first
    08 Mar 2016: accept polynomials of mu of any order
    """
    # Return sparse matrices, A and B, where mu*A*x + B*x = 0
    import scipy.sparse as sps
    Am = len(mets_nonlin)
    Bm = Am
    N = len(me.reactions)
    # Make constraints
    cons_nl = get_nonlin_triplets_general(me, mets_nonlin, growth_symbol, scaleUnits, unitDict)
    dataA = [a for m,r,a,b in cons_nl if a != 0]
    dataB = [b for m,r,a,b in cons_nl if b != 0]
    # Re-order rxns so that nonlin variables appear first (vars in A)
    rxnsnl = list(set([r for m,r,a,b in cons_nl if a != 0]))
    print(len(rxnsnl), 'variables in nonlinear constraints')

    row_indA = [mets_nonlin.index(m) for m,r,a,b in cons_nl if a != 0]
    col_indA = [me.reactions.index(r) for m,r,a,b in cons_nl if a != 0]
    row_indB = [mets_nonlin.index(m) for m,r,a,b in cons_nl if b != 0]
    col_indB = [me.reactions.index(r) for m,r,a,b in cons_nl if b != 0]

    A = sps.csc_matrix( (dataA, (row_indA, col_indA)), shape=(Am,N) )
    B = sps.csc_matrix( (dataB, (row_indB, col_indB)), shape=(Bm,N) )
    d = [met._bound for met in mets_nonlin]
    csense_nonlin = [met._constraint_sense for met in mets_nonlin]

    return A,B,d,csense_nonlin,g


def get_nonlin_triplets(model, mets_nonlin, growth_symbol='mu',
        scaleUnits=False, unitDict=None,
        param_values=None):
    import numpy as np
    import sympy as sp
    from sympy import symbols, Poly

    mu = symbols(growth_symbol)
    if param_values is None:
        param_values = {}
        param_values[growth_symbol] = mu

    cons_nl = []

    # Constraint are in the form, sum_j(mu*A[i,j]*x[j] + B[i,j]*x[j]) = 0   forall i in NLrows
    for met in mets_nonlin:
        for rxn in met.reactions:
            s = str(rxn.metabolites[met])
            if scaleUnits:
                #if met in unitDict['rowsB'] and rxn in unitDict['colsE']:
                #if met in unitDict['rowsB'] and type(rxn) in unitDict['typeE']:
                if met in unitDict['rows_compl'] and type(rxn) in unitDict['typeE']:
                    s = "(" + s + ")*%g" % unitDict['e_mult']

            expr = eval(s,{},param_values)
            if np.isscalar(expr) or expr.is_Number:
                if not expr==0:
                    coeffs = [0.0, expr]
                    cons_nl.append([met, rxn] + [float(c) for c in coeffs])
            else:
                p = Poly(expr)
                coeffs = p.all_coeffs()
                if len(coeffs)==1:
                    coeffs = [0.0] + coeffs
                cons_nl.append([met, rxn] + [float(c) for c in coeffs])

    return cons_nl


def get_nonlin_triplets_general(model, mets_nonlin, growth_symbol='mu',
        scaleUnits=False, unitDict=None,
        param_values=None):
    """
    cons_nl = get_nonlin_triplets_general(model, mets_nonlin, growth_symbol='mu',
        scaleUnits=False, unitDict=None, param_values=None) 

    Generalized functions for nonlinear constraints.
    """
    import numpy as np
    import sympy as sp
    from sympy import symbols, Poly

    mu = symbols(growth_symbol)
    if param_values is None:
        param_values = {}
        param_values[growth_symbol] = mu

    #cons_nl = []
    cons_nl = {}

    # Constraint are in the form, sum_j(gk(mu)*A[i,j]*x[j] + B[i,j]*x[j]) = 0   forall i in NLrows
    # where gk(mu) is a functional family for mu
    for met in mets_nonlin:
        for rxn in met.reactions:
            s = str(rxn.metabolites[met])
            if scaleUnits:
                if met in unitDict['rows_compl'] and type(rxn) in unitDict['typeE']:
                    s = "(" + s + ")*%g" % unitDict['e_mult']

            expr = eval(s,{},param_values)
            if np.isscalar(expr) or expr.is_Number:
                if not expr==0:
                    coeffs = [0.0, expr]
                    cons_nl.append([met, rxn] + [float(c) for c in coeffs])
            else:
                # Not always Poly
                try:
                    p = Poly(expr)
                    coeffs = p.all_coeffs()
                    if len(coeffs)==1:
                        coeffs = [0.0] + coeffs
                    #cons_nl.append([met, rxn] + [float(c) for c in coeffs])
                    cons_nl[(met,rxn)] = [float(c) for c in coeffs]
                except:
                    print('Not expressable as Poly:')
                    print(expr)

    return cons_nl


def make_linear_constraints(me, mets_linear, scaleUnits=False, unitDict=None):
    # Return sparse stoich matrix
    import scipy.sparse as sps
    #from cobra import DictList
    Sm = len(mets_linear)
    N = len(me.reactions)
    # Returns met, rxn, stoich
    mrs = get_stoich_triplets(me, mets_linear)
    if scaleUnits:
        for ind, mrsi in enumerate(mrs):
            if mrsi[2] in unitDict['rows_compl'] and type(r) in unitDict['typeE']:
                mrs[ind] = s*unitDict['e_mult'] 

    data = [s for m,r,s in mrs]
    # Major bottleneck: list.index(m)
    # Overcome by using cobra.DictList !
    # mets_linear_dl = DictList(mets_linear)
    # row_ind = [mets_linear_dl.index(m) for m,r,s in mrs]
    row_ind = [mets_linear.index(m) for m,r,s in mrs]
    col_ind = [me.reactions.index(r) for m,r,s in mrs]
    S = sps.csc_matrix( (data, (row_ind, col_ind) ), shape=(Sm,N))
    b = [met._bound for met in mets_linear]
    csense = [met._constraint_sense for met in mets_linear]

    return S,b,csense



def get_stoich_triplets(model,metabolites):
    return [ (met, rxn, rxn.metabolites[met]) for met in metabolites for rxn in met.reactions ]



def writeME_NLP_from_file(infile = 'scaledME_NLP.mat', outname=None):
    # Load ME model from matfile
    import numpy as np
    import scipy as sp
    import scipy.sparse as sps
    import scipy.io as sio

    mat = sio.loadmat(infile)
    model = mat['model']
    A  = model['A'][0][0]    # [0 Abar]    1073 x 1755
    B  = model['B'][0][0]    # [0 Bbar]    1073 x 1755
    S  = model['S'][0][0]    # [s Sbar],    366 x 1755     s = -e52
    b  = model['b'][0][0]    # 1439 vector  0
    c  = model['c'][0][0]    # 1755 vector  e1
    xl = model['lb'][0][0]   # 1755 vector
    xu = model['ub'][0][0]   # 1755 vector
    J, nnCon, nnJac, neJac, nzS, P, I, V, bl, bu  = makeME_NLP(A,B,S,b,c,xl,xu)
    if outname is not None:
        print('Writing to file: ', outname)
        dumpMat(J, outname, nnCon, nnJac, neJac, nzS, P, I, V, bl, bu)
        print('Done!')

    return J, nnCon, nnJac, neJac, bl, bu



def makeME_LP(S, b, c, xl, xu, csense):
    """
    Create simple LP for qMINOS and MINOS
    Inputs:
    nlp_compat  Make matrices compatible with NLP so that basis can
                be used to warm start NLP by setting 
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
    # Objective row has free bounds
    sl[m-1] = -bigbnd
    su[m-1] = bigbnd

    bl = sp.vstack([xl, sl])
    bu = sp.vstack([xu, su])

    return J, ne, P, I, V, bl, bu



def makeME_LP_for_NLP(A,B,S,b,c,xl,xu):
    """
    Create LP whose basis is compatible with ME-NLP
    14 Aug 2015: first version
    """
    import numpy as np
    import scipy as sp
    import scipy.sparse as sps
    import time

    # Unlike the LP, NLP-compatible version includes slack variables
    # for linear and nonlinear constraints
    # Also, the first column is mu (x[0] = mu)
    #
    #          mu*A*x + w = 0
    #             B*x - w = 0
    #             S*x     = b
    #  -inf <=    c*x    <= inf  (last row so MINOS treats problem as LP)

    # Seems to be infeasible... thus, create from scratch
    #J,nnCon,nnJac,neJac,ne,P,I,V,bl,bu = makeME_NLP(A,B,S,b,c,xl,xu)
    #//--------------------------------------------------------
    mA,nA = A.shape
    mB,nB = B.shape
    mS,nS = S.shape
    nnCon = mA
    nlCon = mB + mS
    mCon  = nnCon + nlCon
    # These values are for NLP --------------------------------------------------------
    nnObj = 1
    nnJac = nA
    neJac = nnCon + A.nnz
    #//--------------------------------------------------------
    n = nA + mA
    e = sps.csc_matrix(np.ones((mA,1)) )
    z = sps.csc_matrix((mB,1))
    s = S[:,0]
    Z = sps.csc_matrix((mS,mA))
    Iw = sps.eye(nnCon).tocsc()
    # What was the Jacobian matrix for NLP must be constraint matrix for LP
    J = sps.vstack((
        sps.hstack((A, Iw)),
        sps.hstack((B,-Iw)),
        sps.hstack((S, Z ))
        )).tocsc()
    J.sort_indices()

    bigbnd  = 1e+40
    wl      = -bigbnd*np.ones((mA,1))
    wu      =  bigbnd*np.ones((mA,1))
    sl      =  np.zeros((mCon,1))
    su      =  np.zeros((mCon,1))
    bl      =  sp.vstack([xl, wl, sl])
    bu      =  sp.vstack([xu, wu, su])

    m,n = J.shape
    ne = J.nnz
    # 1-based indexing for Fortran
    I = [i+1 for i in J.indices]
    V = J.data
    P = [pi+1 for pi in J.indptr]

    #//--------------------------------------------------------
    # Need to add one last free row (slacks unbounded) so that
    # MINOS treats problem as LP
    rowc = sps.hstack( (c, sps.csc_matrix( (1,nnCon) )) )
    J = sps.vstack((J, rowc)).tocsc()
    bigbnd =   1e+40
    bl = np.vstack( (bl, -bigbnd) )
    bu = np.vstack( (bu,  bigbnd) )

    m,n = J.shape
    ne = J.nnz
    I = [i+1 for i in J.indices]
    V = J.data
    P = [pi+1 for pi in J.indptr]

    return J, ne, P, I, V, bl, bu



def makeME_VA(S,b,xl,xu,csense,obj_inds,obj_coeffs):
    """
    Creates ME_LP data for qvaryME, solved using qMINOS with warm-start
    obj_inds: obj column indices
    obj_coeffs: explicitly state -1 or 1 for min or max
    Thus, obj_inds & obj_coeffs have 2*n elements
    [LY] 11 Aug 2015: first version
    """
    import numpy as np
    import scipy as sp
    import scipy.sparse as sps
    import time

    # qMINOS requires an extra row holding the objective
    # Put ANY non-zero for all columns that will be min/maxed 
    Sm,Sn = S.shape
    c = [0. for j in range(0,Sn)]
    for j,v in zip(obj_inds, obj_coeffs):
        c[j] = 1.0

    tic = time.time()
    J = sps.vstack((
        S,
        c)
        ).tocsc()
    toc = time.time() - tic
    print('Stacking J took %f seconds' % toc)

    # Sort indices
    J.sort_indices()

    b2 = b + [0.0]
    m,n = J.shape
    ne = J.nnz
    # Finally, make the P, I, J, V, as well
    # Row indices: recall fortran is 1-based indexing
    tic = time.time()
    I = [i+1 for i in J.indices]
    V = J.data
    toc = time.time()-tic
    print('Making I & V took %f seconds' % toc)

    # Pointers to start of each column
    tic = time.time()
    # Just change to 1-based indexing for Fortran
    P = [pi+1 for pi in J.indptr]
    toc = time.time() - tic
    print('Making P took %f seconds' % toc)

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
    # Objective row has free bounds
    sl[m-1] = -bigbnd
    su[m-1] = bigbnd

    tic = time.time()
    bl = sp.vstack([xl, sl])
    bu = sp.vstack([xu, su])
    toc = time.time()-tic
    print('Stacking bl & bu took %f seconds' % toc)

    obj_indsf = [i+1 for i in obj_inds]

    return J, ne, P, I, V, bl, bu, obj_indsf


def makeME_NLP(A, B, S, b, c, xl, xu):
    """
    Creates ME NLP data matrix to be passed to qsolveME
    Basically the same matrix as dumpMat, but returns object 
    in RAM instead of writing to a file
    [LY] 21 Jul 2015: first version
    """
    import numpy as np
    import scipy as sp
    import scipy.sparse as sps

    mA, nA = A.shape
    mB, nB = B.shape            # 1073 1755
    mS, nS = S.shape            #  366 1755

    nnCon   = mA;                 # 1073
    nlCon   = mB + mS;            # 1439
    m       = nnCon + nlCon;      # 1073 + 1439 = 2512
    nnObj   = 1;                  #    1   mu = x(1)
    nnJac   = nA;                 # 1755    x = [mu; xbar]
    neJac   = nnCon + A.nnz;     # dense col 1 + nnz(Abar)
    n       = nA + mA;            # 2828    [x; w]
    nb      = n + m;              # 3901
    ne      = neJac + B.nnz + S.nnz + mA + mB;

    # Let x = [mu; xbar].
    # The nonlinear constraint function f(mu,xbar) = mu*A*x = mu*Abar*xbar.
    # The nonlinear Jacobian is [Abar*xbar mu*Abar], where we will treat
    # the first column Abar*xbar as dense.
    # Define the full Jacobian structure for the contraints
    #
    #    mu*A*x + w = 0
    #       B*x - w = 0
    #       S*x     = b,
    #
    # where mu is the first variable and A, B both have empty first column.

    e = sps.csc_matrix(np.ones((mA,1)))
    z = sps.csc_matrix((mB,1))
    s = S[:,0]
    Z = sps.csc_matrix((mS,mA))
    Iw = sps.eye(nnCon).tocsc()

    # Work with the J matrix (S from matlab version is J)
    J = sps.vstack( (
            sps.hstack([e, A[:,1:], Iw]),
            sps.hstack([z, B[:,1:], -Iw]),
            sps.hstack([s, S[:,1:], Z]))
        ).tocsc()

    J.sort_indices()
    # Construct bounds for expanded problem
    # xl and xu are already correct for x.
    # We need bounds for w and the slacks s.
    # For now, treat all constraints as equalities (E)
    # and assume b=0.  This makes it easy.

    bigbnd =   1e+40
    wl     = - bigbnd*np.ones((mA,1))
    wu     =   bigbnd*np.ones((mA,1))
    sl     =   np.zeros((m,1))
    su     =   np.zeros((m,1))

    bl     = sp.vstack([xl, wl, sl])
    bu     = sp.vstack([xu, wu, su])

    # Finally, make the P, I, J, V, as well
    m,n = J.shape
    ne = J.nnz
    # Need column-wise elements but nonzero returns row-wise
    # Row indices: recall fortran is 1-based indexing
    # 1-based for Fortran
    I = [i+1 for i in J.indices]
    V = J.data
    # Pointers to start of each column
    P = [pi+1 for pi in J.indptr]

    return J, nnCon, nnJac, neJac, ne, P, I, V, bl, bu



def dumpMat(S, name, nnCon, nnJac, neJac, nzS, P, I, V, bl=None, bu=None):
    """
    Outputs sparse matrix S column-wise
            followed by bl
            followed by bu
    as a flat file 'tinyME.txt' (or outname)
    for input to a Fortran program (say).
    First x is assumed to be mu
    S,bl,bu are assumed to be for constraints of the form
            Sx - s = 0,   bl <= [x;s] <= bu.
    nnCon and nnJac are the dimensions and nnz of a Jacobian
    that is assumed to be in the top-left corner of "S".

    06 Nov 2014: First version derived from part of dumpLP.m.
    30 Dec 2014: Also output nnCon, nnJac (first) and bounds, bl, bu (last).
    20 Jul 2015: First Python port

    name must not have trailing blanks.
    Grab at most 8 characters.
    """
    import numpy as np
    import scipy as sp
    import scipy.sparse as sps

    # print('Testing: you made it to dumpMat call')
    fname = name + '.txt'

    m, n = S.shape
#     nzS = S.nnz
#     #I,J = S.nonzero()
#     # Need column-wise elements but nonzero returns row-wise
#     IJ = [(row,col) for col in xrange(0,n) for row in S[:,col].nonzero()[0]]
#     I = [i for i,j in IJ]
#     J = [j for i,j in IJ]
#     V = [S[i,j] for i,j in IJ]
#     # Pointers to start of each column
#     #************************************************************
#     # LY: a direct port from matlab. May revise into more efficient way later.
#     # NOTE: qsolveME accepts 1-based pointers
#     P = np.zeros((n+1,1))
#     p = 1
#     for j in xrange(0,n):
#         P[j] = p
#         p = p + S[:,j].nnz
# 
#     P[n] = p
    #************************************************************
    # Write to file
    with open(fname, 'w') as f:
        f.write('%.8s\n' % name)      # First line, up to 8 chars
        f.write('%8i\n' % nnCon)      # Num of rows in Jacobian
        f.write('%8i\n' % nnJac)      # Num of cols in Jacobian
        f.write('%8i\n' % neJac)      # Num of nonzeros in Jacobian
        f.write('%8i\n' % m)          # Num of rows in S (including obj row)
        f.write('%8i\n' % n)          # Num of cols in S
        f.write('%8i\n' % nzS)        # Num of nonzeros in S
        # Write as 1-based indexing
        for Pi in P:
            f.write('%8i\n'%Pi)     # Pointers      n+1 
        for Ii in I:
            f.write('%8i\n'%(Ii))      # Row indices   nzS 
        for Vi in V:
            f.write('%19.12e\n'%Vi)  # Values nzS; roughly double precision
        if bl is not None:
            for bli in bl:
                f.write('%19.12e\n'%bli)
        if bu is not None:
            for bui in bu:
                f.write('%19.12e\n'%bui)

    print('Created file: %s'%fname)
    print('File successfully closed? ', f.closed)



def get_lin_nonlin_mets(me, growth_symbol='mu'):
    """
    mets_linear, mets_nonlin = get_lin_nonlin_mets(me, growth_symbol='mu')

    Return linear and nonlinear mets (stoichiometry is function of mu).
    Growth symbol is typicall 'growth_rate_in_per_hour' for ME 1.0
    and 'mu' for ME 2.0 models.
    """
    from sympy import Basic

    mets_nonlin = [met for met in me.metabolites if \
            any([isinstance(rxn.metabolites[met],Basic) for rxn in met.reactions])]
    mets_linear = list(set(me.metabolites).difference(set(mets_nonlin)))
    mets_nonlin = DictList(mets_nonlin)
    mets_linear = DictList(mets_linear)

    return (mets_linear, mets_nonlin)


def split_lin_nonlin(me, growth_symbol='mu'):
    """
    mets_lin, mets_nonlin, rxns_lin, rxns_nonlin = split_lin_nonlin(me, growth_symbol='mu') 
    Split constraints and variables into linear and nonlinear sets
    """
    from sympy import Basic

    mets_lin, mets_nonlin = get_lin_nonlin_mets(me, growth_symbol)
    rxns_nonlin = [rxn for rxn in me.reactions if \
            any([isinstance(s, Basic) for s in rxn.metabolites.values()])]
    rxns_lin = list(set(me.reactions).difference(set(rxns_nonlin)))

    return (mets_lin, mets_nonlin, rxns_lin, rxns_nonlin)


def decomposeME(me):
    """
    Split M and E rows and columns. Can then do variable scaling.
    Aug 07 2015
    """
    pass


def get_rows_compl(me, typeM, typeE):
    """Get complicating constraints: i.e., mixed M and E fluxes"""

    def is_complicating(met):
        rtypes = [type(r) for r in met.reactions]
        hasM = sum( [(ri in typeM) for ri in rtypes ] ) > 0
        hasE = sum( [(ri in typeE) for ri in rtypes] ) > 0
        hasME = hasM and hasE
        return  hasME

    mets_compl = [met for met in me.metabolites if is_complicating(met)]

    return mets_compl


def get_gene(rxn):
    res = re.findall('[a-z]\d+',rxn)
    gene = ''
    if len(res)>0:
        gene = res[0]
    return gene


def swap_rows(mat, a, b):
    mat_csc = sps.csc_matrix(mat)
    a_idx = np.where(mat_csc.indices == a)
    b_idx = np.where(mat_csc.indices == b)
    mat_csc.indices[a_idx] = b
    mat_csc.indices[b_idx] = a
    return mat_csc.asformat(mat.format)

def swap_cols(mat, a, b):
    mat_csr = sps.csr_matrix(mat)
    a_idx = np.where(mat_csr.indices == a)
    b_idx = np.where(mat_csr.indices == b)
    mat_csr.indices[a_idx] = b
    mat_csr.indices[b_idx] = a
    return mat_csr.asformat(mat.format)
