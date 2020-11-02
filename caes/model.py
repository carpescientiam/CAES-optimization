"""Model for diabatic compressed air energy storage (CAES)."""

import numpy as np
import pandas as pd
import pyomo.environ as po
from pyomo.opt import SolverFactory
import rules as ru
import scipy.linalg
import CoolProp.CoolProp as CP


class CAES:

    def __init__(self, **kwargs):
        """Contruct class instance."""
        self.T0 = kwargs.get('T0', 288.15)
        self.p0 = kwargs.get('p0', 1.013e5)
        self.R0 = 287.05
        self.h0 = kwargs.get('h0', CP.PropsSI('HMASS', 'T', self.T0, 'P',
                                              self.p0, 'air'))
        self.s0 = kwargs.get('s0', CP.PropsSI('SMASS', 'T', self.T0, 'P',
                                              self.p0, 'air'))
        self.z = kwargs.get('z', 1045)
        self.g = kwargs.get('g', 9.81)
        self.eta_cmp = kwargs.get('eta_cmp', 0.85)
        self.eta_tur = kwargs.get('eta_tur', 0.9)
        self.eta_me = kwargs.get('eta_me', 0.99)
        self.eta_mo = kwargs.get('eta_mo', 0.97)
        self.eta_ge = kwargs.get('eta_ge', 0.98)
        self.eta_comb = kwargs.get('eta_comb', 0.98)
        self.eta_ht = kwargs.get('eta_ht', 0.95)
        self.quality_cmp = kwargs.get('quality_cmp', 0.8)
        self.quality_tur = kwargs.get('quality_tur', 0.7)
        self.pi_cav_min = kwargs.get('pi_cav_min', 90)
        self.pi_cav_max = kwargs.get('pi_cav_max', 130)
        self.pi_cav_O_0 = kwargs.get('pi_cav_O_0', (self.pi_cav_max-self.pi_cav_min)/2)
        self.pi_loss = kwargs.get('pi_loss', 0.001)
        
    
    def fit_cmp(self, grid_num):
        """Return numerical fit of the compression mass flow."""
        power = np.linspace(self.P_cmp*0.5, self.P_cmp, grid_num)
        pressure = np.linspace(self.pi_cav_min, self.pi_cav_max, grid_num)
        power_mod = np.repeat(power, grid_num)
        pressure_mod = np.tile(pressure, grid_num)
        mass_flow = self.massflow_cmp(power_mod, pressure_mod)
        data = np.column_stack((power_mod, pressure_mod, mass_flow))
        X, Y = np.meshgrid(power, pressure)
        A = np.c_[data[:, 0], data[:, 1], np.ones(data.shape[0])]
        C_cmp, _, _, _ = scipy.linalg.lstsq(A, data[:, 2])
        Z = C_cmp[0]*X + C_cmp[1]*Y + C_cmp[2]
        results = {}
        results['C_cmp'] = C_cmp
        results['original_data'] = data
        results['X'] = X
        results['Y'] = Y
        results['Z'] = Z

        return results
    
    
    def fit_exp(self, grid_num):
        """Return numerical mass flow fit for the expansion."""
        power = np.linspace(self.P_exp*0.5, self.P_exp, grid_num)
        pressure = np.linspace(self.pi_cav_min, self.pi_cav_max, grid_num)
        power_mod = np.repeat(power, grid_num)
        pressure_mod = np.tile(pressure, grid_num)
        mass_flow = self.massflow_exp(power_mod, pressure_mod)
        data = np.column_stack((power_mod, pressure_mod, mass_flow))
        X, Y = np.meshgrid(power, pressure)
        A = np.c_[data[:, 0], data[:, 1], np.ones(data.shape[0])]
        C_exp, _, _, _ = scipy.linalg.lstsq(A, data[:, 2])
        Z = C_exp[0]*X + C_exp[1]*Y + C_exp[2]
        results = {}
        results['C_exp'] = C_exp
        results['original_data'] = data
        results['X'] = X
        results['Y'] = Y
        results['Z'] = Z

        return results
    
    
    def fit_q_comb(self, grid_num):
        """Return numerical fit for combustion."""
        power = np.linspace(self.P_exp*0.5, self.P_exp, grid_num)
        pressure = np.linspace(self.pi_cav_min, self.pi_cav_max, grid_num)
        power_mod = np.repeat(power, grid_num)
        pressure_mod = np.tile(pressure, grid_num)
        q_comb = self.q_comb(power_mod, pressure_mod)
        data = np.column_stack((power_mod, pressure_mod, q_comb))
        X, Y = np.meshgrid(power, pressure)
        A = np.c_[data[:, 0], data[:, 1], np.ones(data.shape[0])]
        C_comb, _, _, _ = scipy.linalg.lstsq(A, data[:, 2])
        Z = C_comb[0]*X + C_comb[1]*Y + C_comb[2]
        results = {}
        results['C_comb'] = C_comb
        results['original_data'] = data
        results['X'] = X
        results['Y'] = Y
        results['Z'] = Z

        return results
    
    def coefficents_linear_model(self, grid_num):
        """Return coefficients for linear approximation."""
        fit_cmp = self.fit_cmp(grid_num)
        fit_exp = self.fit_exp(grid_num)
        fit_comb = self.fit_q_comb(grid_num)
        dc = {}
        dc['cmp_a'] = fit_cmp['C_cmp'][2]
        dc['cmp_b'] = fit_cmp['C_cmp'][0]
        dc['cmp_c'] = fit_cmp['C_cmp'][1]
        dc['cmp_d'] = self.eta_me*self.eta_mo*self.eta_ht
        dc['cmp_e'] = self.q_cmp(self.P_cmp, 90)[1]
        dc['exp_a'] = fit_exp['C_exp'][2]
        dc['exp_b'] = fit_exp['C_exp'][0]
        dc['exp_c'] = fit_exp['C_exp'][1]
        dc['exp_d'] = fit_comb['C_comb'][2]
        dc['exp_e'] = fit_comb['C_comb'][0]
        dc['exp_f'] = fit_comb['C_comb'][1]
      
        return dc
    

        
    
    def optimize(self, el_price, cost, grid_num):
        
        coeff = self.coefficents_linear_model(grid_num)
        costs = pd.read_csv(cost+'.csv', sep=",",index_col=0).astype(np.float64)
        seq = pd.read_csv(el_price+'.csv', index_col=0).astype(np.float64)
       
        m = po.ConcreteModel()
        m.T = po.Set(initialize=seq.index.values)
        
        m.cmp_P_max = po.Param(initialize=self.P_cmp)
        m.cmp_P_min = po.Param(initialize=self.P_cmp*0.5)
        m.cmp_a = po.Param(initialize=coeff['cmp_a'])
        m.cmp_b = po.Param(initialize=coeff['cmp_b'])
        m.cmp_c = po.Param(initialize=coeff['cmp_c'])
        m.cmp_d = po.Param(initialize=coeff['cmp_d'])
        m.cmp_e = po.Param(initialize=coeff['cmp_e'])
        
        m.C_var_cmp = costs.loc[('c_var_cmp'), 'value']
        m.C_var_exp = costs.loc[('c_var_exp'), 'value']
        m.C_fuel = costs.loc[('c_fuel'), 'value']
        m.C_emi = costs.loc[('c_emi'), 'value']
        m.C_charge = costs.loc[('c_charges'), 'value']
        
        m.exp_P_max = po.Param(initialize=self.P_exp)
        m.exp_P_min = po.Param(initialize=self.P_exp*0.5)
        m.exp_a = po.Param(initialize=coeff['exp_a'])
        m.exp_b = po.Param(initialize=coeff['exp_b'])
        m.exp_c = po.Param(initialize=coeff['exp_c'])
        m.exp_d = po.Param(initialize=coeff['exp_d'])
        m.exp_e = po.Param(initialize=coeff['exp_e'])
        m.exp_f = po.Param(initialize=coeff['exp_f'])
        
        m.cas_m0 = po.Param(initialize=self.m_cav_0)
        m.cas_Pi_o_0 = po.Param(initialize=(self.pi_cav_max-self.pi_cav_min)/2)
        m.cas_Pi_min = po.Param(initialize=self.pi_cav_min)
        m.cas_Pi_o_max = po.Param(initialize= self.pi_cav_max - self.pi_cav_min)
        m.eta = po.Param(initialize=self.pi_loss)
        
        m.mkt_C_el = po.Param(m.T, initialize=dict(zip(seq.index.values,
                                               seq['2030NEPC'].values)))
        
       
        m.cmp_P = po.Var(m.T, domain=po.NonNegativeReals,
                 bounds=(0, self.P_cmp))
        m.cmp_m = po.Var(m.T, domain=po.NonNegativeReals)
        m.cmp_Q = po.Var(m.T, domain=po.NonNegativeReals)
        m.cmp_y = po.Var(m.T, domain=po.Binary)
        m.cmp_z = po.Var(m.T, domain=po.NonNegativeReals)
        m.exp_P = po.Var(m.T, domain=po.NonNegativeReals,
                 bounds=(0, self.P_exp))
        m.exp_m = po.Var(m.T, domain=po.NonNegativeReals)
        m.exp_Q = po.Var(m.T, domain=po.NonNegativeReals)
        m.exp_y = po.Var(m.T, domain=po.Binary)
        m.exp_z = po.Var(m.T, domain=po.NonNegativeReals)
        m.cas_Pi_o = po.Var(m.T, domain=po.NonNegativeReals,
                    bounds=(0, self.pi_cav_max - self.pi_cav_min))


        m.profit_test = po.Objective(sense=po.minimize, rule=ru.obj)

        m.boundary = po.Constraint(m.T, rule=ru.boundary)
        m.cas_pi = po.Constraint(m.T, rule=ru.cas_pi)

        m.cmp_p_range_min = po.Constraint(m.T, rule=ru.cmp_p_range_min)
        m.cmp_p_range_max = po.Constraint(m.T, rule=ru.cmp_p_range_max)

        m.massflow_cmp = po.Constraint(m.T, rule=ru.massflow_cmp)
        m.q_cmp = po.Constraint(m.T, rule=ru.q_cmp)
        m.cmp_z1 = po.Constraint(m.T, rule=ru.cmp_z1)
        m.cmp_z2 = po.Constraint(m.T, rule=ru.cmp_z2)
        m.cmp_z3 = po.Constraint(m.T, rule=ru.cmp_z3)
        m.cmp_z4 = po.Constraint(m.T, rule=ru.cmp_z4)

        m.exp_p_range_min = po.Constraint(m.T, rule=ru.exp_p_range_min)
        m.exp_p_range_max = po.Constraint(m.T, rule=ru.exp_p_range_max)

        m.massflow_exp = po.Constraint(m.T, rule=ru.massflow_exp)
        m.q_exp = po.Constraint(m.T, rule=ru.q_exp)
        m.exp_z1 = po.Constraint(m.T, rule=ru.exp_z1)
        m.exp_z2 = po.Constraint(m.T, rule=ru.exp_z2)
        m.exp_z3 = po.Constraint(m.T, rule=ru.exp_z3)
        m.exp_z4 = po.Constraint(m.T, rule=ru.exp_z4)

        m.cmp_exp_excl = po.Constraint(m.T, rule=ru.cmp_exp_excl)

        opt = SolverFactory('gurobi')
        opt.options["mipgap"]=0.02
        results = opt.solve(m, tee=True )
        m.solutions.load_from(results)

        data = {
        'cmp_P': [m.cmp_P[t].value for t in m.T],
        'cmp_Q': [m.cmp_Q[t].value for t in m.T],
        'cmp_y': [m.cmp_y[t].value for t in m.T],
        'exp_P': [m.exp_P[t].value for t in m.T],
        'exp_y': [m.exp_y[t].value for t in m.T],
        'exp_Q': [m.exp_Q[t].value for t in m.T],
        'cas_Pi_o': [m.cas_Pi_o[t].value for t in m.T],
        'cmp_m': [m.cmp_m[t].value for t in m.T],
        'exp_m': [m.exp_m[t].value for t in m.T]       
        }

        df = pd.DataFrame.from_dict(data)
        df.to_csv('results.csv', sep=',')
        
class Diabatic(CAES):

    def __init__(self, V_cas=None, P_cmp=None, P_exp=None, recuperation=True,
                 **kwargs):
        """Contruct class instance of the parent class."""
        CAES.__init__(self)
        """Contruct class instance."""
        self.V_cas = V_cas  # m^3
        self.P_cmp = P_cmp  # MW
        self.P_exp = P_exp  # MW
        self.Hu = kwargs.get('Hu', 52013000)
        self.m_cav_0 = (self.p0*self.V_cas)/(self.R0*self.T0)
        self.recuperation = recuperation
        self.pi_1_cmp = kwargs.get('pi_1_cmp', 6)
        self.pi_1_exp = kwargs.get('pi_1_exp', 11.63)
        self.pi_ex = kwargs.get('pi_ex', 1.1)
        self.T1_s_cmp = kwargs.get('T1_s_cmp', 338.15)
        self.T2_s_cmp = kwargs.get('T2_s_cmp', 338.15)
        self.T2_exp = kwargs.get('T2_exp', 763.15)
        self.T2_s_exp = kwargs.get('T2_s_exp', 338.15)
        self.T1_exp = kwargs.get('T1_exp', 1218.15)
        self.T_out = kwargs.get('T_rec_out', 372.15)

    def temperature_cmp(self, P_cmp, pi_cav):
        """Return the fluid temperature after each compression stage."""
        s0 = self.s0
        h0 = self.h0
        h1_is = CP.PropsSI('HMASS', 'S', s0, 'P', self.pi_1_cmp*self.p0, 'air')
        delta_h1_is = h1_is - h0
        s1_s = CP.PropsSI('SMASS', 'T', self.T1_s_cmp,
                          'P', self.pi_1_cmp*self.p0, 'air')
        h1_s = CP.PropsSI('HMASS', 'S', s1_s, 'P', self.pi_1_cmp*self.p0,
                          'air')
        h2_is = CP.PropsSI('HMASS', 'S', s1_s, 'P', pi_cav*self.p0, 'air')
        delta_h2_is = h2_is - h1_s
        T1_cmp = CP.PropsSI('T', 'HMASS', h0 + (delta_h1_is/((P_cmp/self.P_cmp*(1 - self.quality_cmp) + self.quality_cmp)*self.eta_cmp)), 'P', self.pi_1_cmp*self.p0, 'air')
        T2_cmp = CP.PropsSI('T', 'HMASS', h1_s + (delta_h2_is/((P_cmp/self.P_cmp*(1-self.quality_cmp)+self.quality_cmp)*self.eta_cmp)), 'P', pi_cav*self.p0, 'air')

        return T1_cmp, T2_cmp

    def massflow_cmp(self, P_cmp, pi_cav):
        """Return compression mass flow."""
        s0 = self.s0
        h0 = self.h0
        h1_is = CP.PropsSI('HMASS', 'S', s0,
                           'P', self.pi_1_cmp*self.p0, 'air')
        s1_s = CP.PropsSI('SMASS', 'T', self.T1_s_cmp,
                          'P', self.pi_1_cmp*self.p0, 'air')
        h1_s = CP.PropsSI('HMASS', 'S', s1_s,
                          'P', self.pi_1_cmp*self.p0, 'air')
        h2_is = CP.PropsSI('HMASS', 'S', s1_s, 'P', pi_cav*self.p0, 'air')
        delta_h1_is = h1_is - h0
        delta_h2_is = h2_is - h1_s
        mass_flow = (((P_cmp*10e5)*((P_cmp/self.P_cmp)*(1-self.quality_cmp)+self.quality_cmp)*self.eta_cmp*self.eta_me*self.eta_mo)/(delta_h1_is + delta_h2_is))

        return mass_flow


    def q_cmp(self, P_cmp, pi_cav):
        """Return excess heat from compression."""
        h = CP.PropsSI('HMASS', 'T', self.T2_s_exp, 'P', pi_cav*self.p0, 'air')
        h_0 = CP.PropsSI('HMASS', 'T', self.T0, 'P', self.p0, 'air')
        q_cmp = (P_cmp*self.eta_me*self.eta_mo*self.eta_ht-(self.massflow_cmp(P_cmp, pi_cav)*(h-h_0))/10e5)
        e_cmp = (h-h_0)/10e5

        return q_cmp, e_cmp

    def eta_exergy_cmp(self, P_cmp, pi_cav):

        h = CP.PropsSI('HMASS', 'T', self.T2_s_cmp, 'P', pi_cav*self.p0, 'air')
        s = CP.PropsSI('SMASS', 'T', self.T2_s_cmp, 'P', pi_cav*self.p0, 'air')
        e_st = (h-self.h0)-self.T0*(s-self.s0)-self.g*self.z
        E_st = self.massflow_cmp(P_cmp, pi_cav)*e_st
        eta_ex = E_st/(P_cmp*10e5)

        return eta_ex

    def temperature_exp(self, P_exp, pi_cav):
        """Return the fluid temperature after each expansion stage."""
        s2 = CP.PropsSI('SMASS', 'T', self.T2_exp, 'P', pi_cav*self.p0, 'air')
        h2 = CP.PropsSI('HMASS', 'S', s2, 'P', pi_cav*self.p0, 'air')
        h1_s_is = CP.PropsSI('HMASS', 'S', s2,
                             'P', self.pi_1_exp*self.pi_ex*self.p0, 'air')
        s1 = CP.PropsSI('SMASS', 'T', self.T1_exp,
                        'P', self.pi_1_exp*self.pi_ex*self.p0, 'air')
        h1 = CP.PropsSI('HMASS', 'S', s1,
                        'P',  self.pi_1_exp*self.pi_ex*self.p0, 'air')
        h_ex_is = CP.PropsSI('HMASS', 'S', s1, 'P', self.pi_ex*self.p0, 'air')
        delta_h2_is = h2 - h1_s_is
        delta_h1_is = h1 - h_ex_is
        T1_s_exp = CP.PropsSI('T', 'HMASS', h2-(delta_h2_is*(np.sqrt(P_exp/self.P_exp)*(1-self.quality_tur)+self.quality_tur)*self.eta_tur), 'P', self.pi_1_exp*self.pi_ex*self.p0, 'air')
        Tex_exp = CP.PropsSI('T', 'HMASS', h1-(delta_h1_is*(np.sqrt(P_exp/self.P_exp)*(1-self.quality_tur)+self.quality_tur)*self.eta_tur), 'P', self.pi_ex*self.p0, 'air')

        return T1_s_exp, Tex_exp

    def fuel_ratio(self, P_exp, pi_cav, T_start=537):
        """Return m1 for expansion."""
        h1 = CP.PropsSI('HMASS', 'T', self.T1_exp,
                        'P', self.pi_1_exp*self.pi_ex*self.p0, 'air')
        h2 = CP.PropsSI('HMASS', 'T', self.T2_exp, 'P', pi_cav*self.p0, 'air')
        s2 = CP.PropsSI('SMASS', 'T', self.T2_exp, 'P', pi_cav*self.p0, 'air')
        h1_s_is = CP.PropsSI('HMASS', 'S', s2,
                             'P', self.pi_1_exp*self.pi_ex*self.p0, 'air')
        delta_h1_is = h1 - h1_s_is
        h1_s = h2-(delta_h1_is*(np.sqrt(P_exp/self.P_exp)*(1-self.quality_tur)+self.quality_tur)*self.eta_tur)
        delta_h1 = h1-h1_s
        if self.recuperation is True:
            n = 0
            while n < 3:
                h2_s = CP.PropsSI('HMASS', 'T', self.T2_s_exp,
                                  'P', pi_cav*self.p0, 'air')
                s1 = CP.PropsSI('SMASS', 'T', self.T1_exp,
                                'P', self.pi_1_exp*self.pi_ex*self.p0, 'air')
                hex_is = CP.PropsSI('HMASS', 'S', s1,
                                    'P', self.pi_ex*self.p0, 'air')
                delta_hex_is = h1-hex_is
                h_ex = h1-(delta_hex_is*(np.sqrt(P_exp/self.P_exp)*(1-self.quality_tur)+self.quality_tur)*self.eta_tur)
                h_out = CP.PropsSI('HMASS', 'T', self.T_out,
                                   'P', self.p0, 'air')
                h2_rec = CP.PropsSI('HMASS', 'T', T_start,
                                    'P', pi_cav*self.p0, 'air')
                delta_h_rec = h2 - h2_rec
                m2 = 1+(delta_h_rec/(self.eta_comb*self.Hu))
                m1 = (m2+(delta_h1/(self.eta_comb*self.Hu)))
                cp_air = (h2_rec - h2_s)/(T_start-self.T2_s_exp)
                T_rec = (self.eta_ht*m1*(h_ex-h_out)/(cp_air))+self.T2_s_exp
                T_start = T_rec
                n = n+1
        else:
            h2_s = CP.PropsSI('HMASS', 'T', self.T2_s_exp,
                              'P', pi_cav*self.p0, 'air')
            delta_h2 = h2 - h2_s
            m2 = 1 + (delta_h2)/(self.eta_comb*self.Hu)
            m1 = (m2 + (delta_h1/(self.eta_comb*self.Hu)))
            T_rec = self.T2_s_exp

        return m1, m2, T_rec

    def massflow_exp(self, P_exp, pi_cav):
        """Return the mass flow rate for the expansion."""
        s2 = CP.PropsSI('SMASS', 'T', self.T2_exp, 'P', pi_cav*self.p0, 'air')
        h2 = CP.PropsSI('HMASS', 'S', s2, 'P', pi_cav*self.p0, 'air')
        h1_s_is = CP.PropsSI('HMASS', 'S', s2, 'P', self.pi_1_exp*self.pi_ex*self.p0, 'air')
        s1 = CP.PropsSI('SMASS', 'T', self.T1_exp, 'P', self.pi_1_exp*self.pi_ex*self.p0, 'air')
        h1 = CP.PropsSI('HMASS', 'S', s1, 'P',  self.pi_1_exp*self.pi_ex*self.p0, 'air')
        h_ex_is = CP.PropsSI('HMASS', 'S', s1, 'P', self.pi_ex*self.p0, 'air')
        delta_h2_is = h2 - h1_s_is
        delta_h1_is = h1 - h_ex_is
        massflow = (P_exp*1e6)/((np.sqrt(P_exp/self.P_exp)*(1-self.quality_tur)+self.quality_tur)*(self.eta_tur*self.eta_ge*self.eta_me*(self.fuel_ratio(P_exp, pi_cav)[1]*delta_h2_is+self.fuel_ratio(P_exp, pi_cav)[0]*delta_h1_is)))

        return massflow


    def q_comb(self, P_exp, pi_cav):
        """Return combustion heat flow within expansion."""
        h2 = CP.PropsSI('HMASS', 'T', self.T2_exp, 'P', pi_cav*self.p0, 'air')
        s2 = CP.PropsSI('SMASS', 'T', self.T2_exp, 'P', pi_cav*self.p0, 'air')
        h1_s_is = CP.PropsSI('HMASS', 'S', s2, 'P', self.pi_1_exp*self.pi_ex*self.p0, 'air')
        h1 = CP.PropsSI('HMASS', 'T', self.T1_exp, 'P', self.pi_1_exp*self.pi_ex*self.p0, 'air')
        h1_s = h2 - ((h2-h1_s_is)*(np.sqrt(P_exp/self.P_exp)*(1-self.quality_tur)+self.quality_tur)*self.eta_tur)
        delta_h1 = h1 - h1_s

        if self.recuperation is True:
            h2_rec = CP.PropsSI('HMASS', 'T', self.fuel_ratio(P_exp, pi_cav)[2], 'P', pi_cav*self.p0, 'air')
            delta_h_rec = h2 - h2_rec
            q_comb = self.massflow_exp(P_exp, pi_cav)*((delta_h_rec + delta_h1)/self.eta_comb)/(10e5)
        else:   
            h2_s = CP.PropsSI('HMASS', 'T', self.T2_s_exp, 'P', pi_cav*self.p0, 'air')
            delta_h2 = h2 - h2_s
            q_comb = self.massflow_exp(P_exp, pi_cav)*((delta_h2 + delta_h1)/self.eta_comb)/(10e5)

        return q_comb


    def power_heat_ratio(self, P_exp, pi_cav):

        if self.recuperation is True:
            h2_rec = CP.PropsSI('HMASS', 'T', self.fuel_ratio(P_exp, pi_cav)[2], 'P', pi_cav*self.p0, 'air')
            h2_s = CP.PropsSI('HMASS', 'T', self.T2_s_exp, 'P', pi_cav*self.p0, 'air')
            delta_h_rec = h2_rec-h2_s
            phr = (self.massflow_exp(P_exp, pi_cav)*delta_h_rec)/(self.P_exp*10e5)

        else:
            phr = 0

        return phr




    
    
# if __name__ == '__main__':
#     import doctest
#     doctest.testmod()

