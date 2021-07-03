from __future__ import print_function, absolute_import, division
#============================================================
# File me1.py
#
#   class   ME_NLP1
#   def     me1nlp
#
# Methods for handling ME 1.0 model files and make them
# compatible with solveME
#
# Laurence Yang, SBRG, UCSD
#
# 28 Sep 2015:  first version
# 29 Sep 2015:  cdecimal for sig_digs
# 04 Dec 2015:  corrected to initialize values in init
# 17 May 2016:  try not using cDecimal again. Substituting is
#               2x faster but cannot specifiy significant digits
#               in stoich coefficients.
# 18 May 2016:  moved make_dilution_fluxes to ME_NLP (ME_NLP1 inherits)
#============================================================

from qminospy.me2 import ME_NLP, makeME_LP, makeME_LP_for_NLP, makeME_NLP, makeME_VA
from qminospy.me2 import make_nonlin_constraints, make_linear_constraints
from qminospy import qwarmLP
import cobra
if cobra.__version__ < '0.6.0':
    from cobra.core.Solution import Solution
else:
    from cobra.core.solution import LegacySolution as Solution
from cobra import DictList
from cobra import Reaction, Metabolite
from sympy import lambdify, Basic, Symbol

# need this again for me1 backwards compatibility??
#from cdecimal import Decimal, getcontext
import qminospy.me2 as qme
import scipy.sparse as sps
import numpy as np
import time
import warnings
import six
import cobrame
import re

class ME_NLP1(ME_NLP):
    """
    Subclass of ME_NLP for handling ME 1.0 models.
    """
    def __init__(self, me, growth_key='mu'):
        # growth_key: 'mu' for ME 2.0; 'growth_rate_in_per_hour' for ME 1.0
        self.me = me
        # 04 Dec 2015:  corrected to initialize values in init
        self.growth_key = growth_key
        self.substitution_dict = {
#                            "growth_rate_in_per_hour": 1.0,
                            "mRNA_mean_lifetime_minutes":5.0,
                            "lol_efficiency_in_per_second":0.9,
                            "bam_efficiency_in_per_second":0.0267,
                            "tat_translocation_efficiency_per_second":0.0125,
                            "secA_translocation_efficiency_aa_per_second":4.0,
                            "proportion_of_rna_that_is_mrna":0.02,
                            "proportion_of_rna_that_is_trna":0.12,
                            "proportion_of_rna_that_is_rrna":0.86,
                            "average_rna_nt_mw":324.,
                            "average_aa_mw":109.,
                            "mass_rrna_per_ribosome":1700.*1000,
                            "average_trna_mw":25000.,
                            "enzyme_efficiency_scaling": 1.0,
                            "percent_nadh_dehydrogenase_that_ndh1": 0.5,
                            "unmodeled_protein_proportion_of_proteome": 0.36,
                            "gam_value" : 34.98,
                            "proportion_of_inner_membrane_sa_that_is_protein":0.75,
                            "proportion_of_outer_membrane_sa_that_is_protein":0.75,
                            "ngam_value":1.}
        self.substitution_dict[self.growth_key] = 1.0
        # Hold the keys in one order for lambdify
        self.subs_keys_ordered = self.substitution_dict.keys()
        # Initially have None for compiled_expressions
        self.compiled_expressions = None
        # Set sigdigs for constraint, bound values
        self.sig_digs = 4
        #----------------------------------------------------
        # Hold LP results and options
        self.lp_inform  = None
        self.lp_hs      = None
        self.hs         = None
        self.feas_basis = None
        self.lp_x       = None
        #----------------------------------------------------
        # Solver options
        self.init_solver_opts()

    def compile_expr(self, expr):
        """Compile expressions with all parameter keys in ME 1.0"""
        # Note that ufuncify too slow to compile
        #f = lambdify(self.subs_keys_ordered, expr) if isinstance(expr, Basic) else expr
        # 19 Jan 2017:  already checked isinstance(expr, Basic) before calling this method
        f = lambdify(self.subs_keys_ordered, expr)
        return f

    def compile_expressions(self, verbosity=0):
        """
        Compile expressions for ME 1.0.
        Use format consistent with cobrame:
        (met_index, rxn_index): stoichiometry, (None, rxn_index): (lower_bound, upper_bound)
        (met_index, None): (met_bound, met_constraint_sense)
        """
        tic = time.time()
        expressions = {}
        me = self.me
        for i, r in enumerate(me.reactions):
            for met, stoic in six.iteritems(r._metabolites):
                if isinstance(stoic, Basic):
                    expressions[(me.metabolites.index(met), i)] = self.compile_expr(stoic)
            # If lower or upper bound symbolic:
            if isinstance(r.lower_bound, Basic) or isinstance(r.upper_bound, Basic):
                expressions[(None,i)] = (self.compile_expr(r.lower_bound),
                                         self.compile_expr(r.upper_bound))
        # Metabolite bound
        for i, metabolite in enumerate(me.metabolites):
            if isinstance(metabolite._bound, Basic):
                expressions[(i, None)] = (self.compile_expr(metabolite._bound),
                        metabolite._constraint_sense)

        toc = time.time() - tic
        if verbosity > 0:
            print('Finished compiling expressions in %f seconds' % toc)

        return expressions

    def varyme(self, mu_fixed, rxns_fva0, basis=None, verbosity=0):
        """
        fva_result, fva_stats = varyme(self, mu_fixed, rxns_fva)

        rxns_fva:  list of reactions to be varied (Reaction objects or ID)

        High-level interface for qvaryME (quad-prec FVA)
        12 Aug 2015: first version. Must fix bugs.
        04 Aug 2016: ME 1.0 version
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

        S, b, xl, xu, csense, c = self.substitute_mu(mu_fixed)

        obj_inds0 = [ me.reactions.index(rxn) for rxn in rxns_fva for j in range(0,2)]
        obj_coeffs = [ci for rxn in rxns_fva for ci in (1.0, -1.0)]
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
                print('Warm-starting first run using basis of length %d'%len(hs))

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
            print('Finished varyME in %f seconds for %d rxns (%d quadLPs)'%(t_elapsed,
                    len(rxns_fva), len(obj_inds)))

        # Return result consistent with cobrame fva
        fva_result = {
            (self.me.reactions[obj_inds0[2*i]].id):{
                'maximum':obj_vals[2*i],
                'minimum':obj_vals[2*i+1] } for i in range(0, nVary//2) }

        # Save updated basis
        self.hs = hs
        self.lp_hs = hs


        return fva_result, fva_stats

    def substitute_mu(self, mu_fix, verbosity=0):
        """
        Use when me is missing the construct_S method
        """
        me = self.me
        if self.compiled_expressions is None:
            # 16 Sep 2016: [LY] should update ordered keys, too, in case the reason
            # compiled_expressions is None is because new symbols were added
            self.subs_keys_ordered = self.substitution_dict.keys()
            self.compiled_expressions = self.compile_expressions(verbosity=verbosity)
        compiled_expressions = self.compiled_expressions

        self.substitution_dict[self.growth_key] = mu_fix
        # Substitution keys need to be the same as when the lambdas were originally
        # compiled
        #
        # Get the subtitution values in the right order for lambdfiy
        sub_vals = [self.substitution_dict[k] for k in self.subs_keys_ordered]

        # Initialize S, lb, ub, b
        S = sps.dok_matrix((len(me.metabolites), len(me.reactions)))
        xl = np.matrix([r.lower_bound for r in me.reactions]).transpose()
        xu = np.matrix([r.upper_bound for r in me.reactions]).transpose()
        b = [0. for m in me.metabolites]    # Gets filled with actual rhs below
        # Fill in all matrix & constraint rhs entries (incl. not mu-dependent)
        tic = time.time()
        #getcontext().prec = self.sig_digs      # from cdecimal
        for mind,met in enumerate(me.metabolites):
            # Fill in constraint bounds: MOSTLY just float
            if hasattr(met._bound, 'subs'):
                expr = self.compiled_expressions[(mind,None)]
                b[mind] = float(expr[0](*sub_vals))
                # Decimal() for backwards compat with ME 1.0
                #b[mind] = float(Decimal(expr[0](*sub_vals)).normalize())   # decimal reallly?
                #b[mind] = float(Decimal(expr[0](*sub_vals)))  # slower than normalize?
            else:
                b[mind] = met._bound
            # Fill in stoichiometries: MOSTLY symbolic, or float? Hard to say.
            for rxn in met.reactions:
                rind = me.reactions.index(rxn)
                if (mind, rind) in self.compiled_expressions:
                    expr = self.compiled_expressions[(mind,rind)]
                    s = float(expr(*sub_vals))
                    # Decimal() for backwards compat with ME 1.0
                    #s = float(Decimal(expr(*sub_vals)).normalize())
                    #s = float(Decimal(expr(*sub_vals)))
                    if not np.isinf(s):
                        S[mind, rind] = s
                else:
                    S[mind,rind] = rxn.metabolites[met]
        # Fill in var bounds: MOSTLY just float
        for rind,rxn in enumerate(me.reactions):
            if hasattr(rxn.lower_bound, 'subs'):
                # Then, there must be a compiled expression
                expr = self.compiled_expressions[(None, rind)]
                xl[rind] = float(expr[0](*sub_vals))
                # Decimal() for backwards compat with ME 1.0
                #xl[rind] = float(Decimal(expr[0](*sub_vals)).normalize())
                #xl[rind] = float(Decimal(expr[0](*sub_vals)))
            else:
                xl[rind] = rxn.lower_bound
            if hasattr(rxn.upper_bound, 'subs'):
                expr = self.compiled_expressions[(None, rind)]
                xu[rind] = float(expr[1](*sub_vals))
                # Using Decimal for sig_digs, even if slower
                # Decimal() for backwards compat with ME 1.0
                #xu[rind] = float(Decimal(expr[1](*sub_vals)).normalize())
                #xu[rind] = float(Decimal(expr[1](*sub_vals)))
            else:
                xu[rind] = rxn.upper_bound
        toc = time.time() - tic

        if verbosity > 0:
            print('Finished substituting S,lb,ub in %f seconds' % toc)

        c = [r.objective_coefficient for r in me.reactions]
        csense = [m._constraint_sense for m in me.metabolites]

        return S,b,xl,xu,csense,c


    def make_lp(self, mu_fix, compiled_expressions=None, verbosity=0):
        """
        Construct LP problem for qMINOS or MINOS. For ME 1.0.
        """
        me = self.me
        if self.compiled_expressions is None:
            # 16 Sep 2016: [LY] should update ordered keys, too, in case the reason
            # compiled_expressions is None is because new symbols were added
            self.subs_keys_ordered = self.substitution_dict.keys()
            self.compiled_expressions = self.compile_expressions(verbosity=verbosity)
        compiled_expressions = self.compiled_expressions

        self.substitution_dict[self.growth_key] = mu_fix
        # Substitution keys need to be the same as when the lambdas were originally
        # compiled
        #
        # Get the subtitution values in the right order for lambdfiy
        sub_vals = [self.substitution_dict[k] for k in self.subs_keys_ordered]

        # Initialize S, lb, ub, b
        S = sps.dok_matrix((len(me.metabolites), len(me.reactions)))
        xl = np.matrix([r.lower_bound for r in me.reactions]).transpose()
        xu = np.matrix([r.upper_bound for r in me.reactions]).transpose()
        b = [0. for m in me.metabolites]
        # Fill in all matrix & constraint rhs entries (incl. not mu-dependent)
        tic = time.time()
        #getcontext().prec = self.sig_digs      # from cdecimal
        for mind,met in enumerate(me.metabolites):
            # Fill in constraint bounds: MOSTLY just float
            if hasattr(met._bound, 'subs'):
                expr = self.compiled_expressions[(mind,None)]
                b[mind] = float(expr[0](*sub_vals))
                # Decimal() for backwards compat with ME 1.0
                #b[mind] = float(Decimal(expr[0](*sub_vals)).normalize())   # decimal reallly?
                #b[mind] = float(Decimal(expr[0](*sub_vals)))  # slower than normalize?
            else:
                b[mind] = met._bound
            # Fill in stoichiometries: MOSTLY symbolic, or float? Hard to say.
            for rxn in met.reactions:
                rind = me.reactions.index(rxn)
                if (mind, rind) in self.compiled_expressions:
                    expr = self.compiled_expressions[(mind,rind)]
                    #****************************************
                    # DEBUG
                    try:
                        s = float(expr(*sub_vals))
                    except TypeError as e:
                        # Just indicate which rxn,met,stoich had issues
                        print(repr(e))
                        print('rxn=%s \t met=%s' % (rxn.id, met.id))
                        print('stoich=', rxn.metabolites[met])
                        raise Exception('Failed to convert symbolic stoichiometry to float')

                        #print('Trying float(Decimal())')
                        #s = float(Decimal(expr(*sub_vals)).normalize())

                    #****************************************
                    # Decimal() for backwards compat with ME 1.0
                    #s = float(Decimal(expr(*sub_vals)).normalize())
                    #s = float(Decimal(expr(*sub_vals)))
                    if not np.isinf(s):
                        S[mind, rind] = s
                else:
                    S[mind,rind] = rxn.metabolites[met]
        # Fill in var bounds: MOSTLY just float
        for rind,rxn in enumerate(me.reactions):
            if hasattr(rxn.lower_bound, 'subs'):
                # Then, there must be a compiled expression
                expr = self.compiled_expressions[(None,rind)]
                xl[rind] = float(expr[0](*sub_vals))
                # Decimal() for backwards compat with ME 1.0
                #xl[rind] = float(Decimal(expr[0](*sub_vals)).normalize())
                #xl[rind] = float(Decimal(expr[0](*sub_vals)))
            else:
                xl[rind] = rxn.lower_bound
            if hasattr(rxn.upper_bound, 'subs'):
                expr = self.compiled_expressions[(None,rind)]
                xu[rind] = float(expr[1](*sub_vals))
                # Using Decimal for sig_digs, even if slower
                # Decimal() for backwards compat with ME 1.0
                #xu[rind] = float(Decimal(expr[1](*sub_vals)).normalize())
                #xu[rind] = float(Decimal(expr[1](*sub_vals)))
            else:
                xu[rind] = rxn.upper_bound
        toc = time.time() - tic

        if verbosity > 0:
            print('Finished substituting S,lb,ub in %f seconds' % toc)

        c = [r.objective_coefficient for r in me.reactions]
        csense = [m._constraint_sense for m in me.metabolites]
        tic = time.time()
        J, ne, P, I, V, bl, bu = makeME_LP(S,b,c,xl,xu,csense)
        toc = time.time()-tic
        if verbosity > 0:
            print('Finished makeME_LP in %f seconds' % toc)

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

    def bisectmu(self, precision=1e-3, mumin=0.0, mumax=2.0, maxIter=100, quad=True,
                 basis=None, nlp_compat=False, check_feas0=False, zero_mu=1e-3,
                 verbosity=2, solver_precision='quad'):
        """
        muopt, hs, xopt, cache = bisectmu(self, precision=1e-3, mumin=0.0, mumax=2.0,
                maxIter=100, quad=True, basis=None, nlp_compat=False, check_feas0=False,
                zero_mu=1e-3, verbosity=2)

        Bisection to maximize mu using qMINOS.
        Sequence of feasibility problems.
        TODO: set quad=False if need fast basis for subsequent quadMINOS
        """
        import copy as cp

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
                    precision=solver_precision, basis=hs)
            if me.solution.status is not 'optimal':
                warnings.warn('Infeasible at mu=%g. Returning.'%zero_mu)
                return zero_mu, hs0, x0, cache
            else:
                hs = hs0


        def checkmu(muf, hs):
            if muf not in cache:
                x_new, stat_new, hs_new = self.solvelp(
                    muf, basis=hs, nlp_compat=nlp_compat, verbosity=verbosity,
                    precision=solver_precision)
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

        tic = time.time()

        if verbosity >= 2:
            print('iter\tmuopt    \ta     \tb     \tmu1       \tstat1')
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

        toc = time.time()-tic
        if verbosity >= 2:
            print('Bisection done in %g seconds' % toc)

        # Save final solution
        me.solution = solution
        # Save feasible basis
        self.feas_basis = hs

        return muopt, hs, xopt, cache


    def bisectme(self, precision=1e-3, mumin=0.0, mumax=2.0, maxIter=100, quad=True, golden=True, basis=None, nlp_compat=False, check_feas0=False, zero_mu=1e-3, verbosity=0,
            solver_precision='quad'):
        """
        Bisection using qMINOS.
        TODO: set quad=False if need fast basis for subsequent quadMINOS
        """
        import copy as cp

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
                    precision=solver_precision)
            if me.solution.status is not 'optimal':
                warnings.warn('Infeasible at mu=%g. Returning.'%zero_mu)
                return zero_mu, hs0, x0, cache
            else:
                hs = hs0

        def checkmu(muf, hs):
            if muf not in cache:
                x_new, stat_new, hs_new = self.solvelp(
                    muf, basis=hs, nlp_compat=nlp_compat, verbosity=verbosity,
                    precision=solver_precision)
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

        tic = time.time()

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

            if verbosity > 1:
                print(iter, muopt, a, b, mu1, mu2, stat1, stat2)

        toc = time.time() - tic
        if verbosity > 0:
            print('Bisection done in %g seconds' % toc)

        # Save final solution
        me.solution = solution

        return muopt, hs, xopt, cache

    def make_nlp(self, growth_rxn='biomass_dilution', growth_symbol='mu'):
        """
        make_nlp()

        Make NLP problem for qMINOS. Uses generalized formulation that is
        compatible with ME 1.0 models:

        max  mu
        s.t. gk(mu)*A*v + B*v = d
             Sv = b
             l <= v <= u

        where gk(mu) are known families of functions of mu.
        """
        import qminospy.me2 as qme

        if self.A is None:
            self.make_matrices()

        #----------------------------------------------------
        # Supported gk(mu):
        # - Poly of any order
        # - ?
        #----------------------------------------------------
        rxn_mu = me.reactions.get_by_id(growth_rxn)
        ind_mu = me.reactions.index(rxn_mu)
        mets_linear, mets_nonlin = qme.get_lin_nonlin_mets(me, growth_symbol)



        #----------------------------------------------------
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


    def solvenlp(self, precision=0.1, max_iter=20, check_feas0=False, zero_mu=1e-3, verbosity=0):
        """
        The combined solution procedure: bisection (golden section) to
        get initial mu0, followed by qsolveME to find mu*.
        Call:
        x, stat, hs = solvenlp(precision, max_iter, check_feas0, zero_mu)
        Parameters:
        precision: precision of the bisectme phase (i.e., epsilon for |b-a|<=epsilon to converge)
        max_iter
        check_feas0: check feasibility at mu = 0?
        zero_mu: terminate if mu is infeasible for mu <= zero_mu,
        """
        if self.nb is None:
            # Must allow formulation more general than mu*A*v + B*v = 0
            # Allow:
            # f(mu)*A*v + B*v = 0
            self.make_nlp()

        hs = None
        # Check feasibility at 0.0?
        if check_feas0:
            x0, stat0, hs0 = self.solvelp(zero_mu, nlp_compat=True, verbosity=verbosity)
            if stat0 is not 'optimal':
                #raise ValueError('Infeasible at mu=0.0. Stopping.')
                warnings.warn('Infeasible at mu=%g. Returning.'%zero_mu)
                return x0, stat0, hs0
            else:
                hs = hs0

        # Bisection (golden section)
        tic1 = time.time()
        mu_bs, hs_bs, x_bs, cache_bs = self.bisectme(precision=precision,
                maxIter=max_iter, nlp_compat=True, basis=hs, verbosity=verbosity)
        time_bs = time.time()-tic1

        # NLP
        if hs_bs is None:
            warnings.warn('Feasible mu0 not found with bisectME. Returning.')
            return x_bs, 'infeasible', hs_bs
        else:
            tic2 = time.time()
            self.mu0 = mu_bs
            x, stat, hs = self.solve(x0=x_bs[0:self.nb], basis=hs_bs[0:self.nb])
            time_nlp = time.time()-tic2

            t_elapsed = time.time()-tic1

            if verbosity > 0:
                print('Finished in %f seconds (%f bisectME, %f ME-NLP)' %
                      (t_elapsed, time_bs, time_nlp))

            return x, stat, hs

    def make_matrices(self):
        """Constructs NLP matrices from me model object"""
        A,B,S,b,c,xl,xu = me1nlp(self.me,
                growth_symbol=self.growth_symbol,
                scaleUnits=self.scaleUnits,
                unitDict=self.unitDict)
        self.A = A
        self.B = B
        self.S = S
        self.b = b
        self.c = c
        self.xl = xl
        self.xu = xu


##//============================================================
## End of ME_NLP1
##//============================================================



def me1nlp(me, growth_symbol='growth_rate_in_per_hour', scaleUnits=False, LB=0.0, UB=1000.0, growth_rxn='growth_rate_in_per_hour', unitDict=None):
    """
    From ME 1.0 model object, create NLP data matrices of the form:
    max mu = c'*x
    mu,x
    s.t.  g(mu)*A*x + B*x = d
          S*x          = b
          xl <= x <= xu
    Should support more general NLP constraints, too.
    """
    import numpy as np
    from cobra.core.DictList import DictList

#    # Reorder rxns such that growth_rxn is the first rxn
#     nRxn = len(me.reactions)
#     rxn_mu = me.reactions.get_by_id(growth_rxn)
#     ind_mu = me.reactions.index(rxn_mu)
#     order_new = range(0,nRxn)
#     order_new[ind_mu] = 0
#     order_new[0] = ind_mu
#     me.reactions = DictList( [ me.reactions[i] for i in order_new] )

    # Identify mets in nonlinear vs. linear constraints
#     mets_nonlin = [met for met in me.metabolites if sum(
#         [growth_symbol in str(rxn.metabolites[met]) for rxn in met.reactions ])
#     ]
    # Oct 1, 2015: faster? a lot faster
#     mets_nonlin = [met for met in me.metabolites if \
#             any([isinstance(rxn.metabolites[met],Basic) for rxn in met.reactions])]
#     mets_linear = list(set(me.metabolites).difference(set(mets_nonlin)))
#     mets_nonlin = DictList(mets_nonlin)
#     mets_linear = DictList(mets_linear)
    mets_linear, mets_nonlin = qme.get_lin_nonlin_mets(me, growth_symbol)

    # Create A, B: nonlinear constraints
    A,B = make_nonlin_constraints(me, mets_nonlin, growth_symbol, scaleUnits, unitDict)
    # Create S & b: linear constraints
    S,b = make_linear_constraints(me, mets_linear, scaleUnits, unitDict)
    # Create c: objective coefficient vector
    # c = [r.objective_coefficient for r in me.reactions]
    c = [0.0 for r in me.reactions]
    c[me.reactions.index(rxn_mu)] = 1.0
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

    return A,B,S,b,c,xl,xu
