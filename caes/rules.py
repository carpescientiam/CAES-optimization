"""Objective and constraint rules for Compressed Air Energy Storages (CAES)."""
import pyomo.environ as po


"objective function"

def obj(m):
    expr = (sum((m.mkt_C_el[t]+ m.C_var_cmp + m.C_charge) * m.cmp_P[t] +
                (m.C_fuel+ m.C_emi) * m.exp_Q[t] + (m.exp_P_max + m.cmp_P_max)*0.57 -
                (m.mkt_C_el[t]-m.C_var_exp) * m.exp_P[t]    
            for t in m.T))
    return expr


"compression rules"


def cmp_p_range_min(m, t):
    return(m.cmp_P[t] >= m.cmp_y[t] * m.cmp_P_min)


def cmp_p_range_max(m, t):
    return(m.cmp_P[t] <= m.cmp_y[t] * m.cmp_P_max)


def massflow_cmp(m, t):
    return(m.cmp_m[t] == m.cmp_a * m.cmp_y[t] + m.cmp_b * m.cmp_P[t] + m.cmp_c * m.cmp_z[t] + m.cmp_c * m.cas_Pi_min * m.cmp_y[t])


def cmp_z1(m, t):
    return(m.cmp_z[t] <= m.cas_Pi_o_max * m.cmp_y[t])


def cmp_z2(m, t):
    return(m.cmp_z[t] <= m.cas_Pi_o[t])


def cmp_z3(m, t):
    return(m.cmp_z[t] >= m.cas_Pi_o[t] - (1 - m.cmp_y[t]) * m.cas_Pi_o_max)


def cmp_z4(m, t):
    return(m.cmp_z[t] >= 0)


def q_cmp (m,t):
    return(m.cmp_Q[t] == m.cmp_d * m.cmp_P[t] - m.cmp_e * m.cmp_m[t])


"expansion rules"


def exp_p_range_min(m, t):
    return(m.exp_P[t] >= m.exp_y[t] * m.exp_P_min)


def exp_p_range_max(m, t):
    return(m.exp_P[t] <= m.exp_y[t] * m.exp_P_max)
    
def massflow_exp(m, t):
    return(m.exp_m[t] == m.exp_a * m.exp_y[t]+ m.exp_b *m.exp_P[t]+ m.exp_c * m.exp_z[t]  + m.exp_c* m.cas_Pi_min * m.exp_y[t])

    
def exp_z1(m, t):
    return(m.exp_z[t] <= m.cas_Pi_o_max * m.exp_y[t])


def exp_z2(m, t):
    return(m.exp_z[t] <= m.cas_Pi_o[t])


def exp_z3(m, t):
    return(m.exp_z[t] >= m.cas_Pi_o[t] - (1 - m.exp_y[t]) * m.cas_Pi_o_max)


def exp_z4(m, t):
    return(m.exp_z[t] >= 0)
      

def q_exp (m,t):
    return(m.exp_Q[t] == m.exp_d * m.exp_y[t] + m.exp_e * m.exp_P[t] + m.exp_f * m.exp_z[t] + m.exp_f * m.cas_Pi_min * m.exp_y[t])


def cmp_exp_excl(m, t):
    return(m.cmp_y[t] + m.exp_y[t] <= 1)

# -----------------------------------------------------------------------------
#                          CAVERN RULES
# -----------------------------------------------------------------------------

def cas_pi(m, t):
    if t == min(m.T):
        return(m.cas_Pi_o[t] == m.cas_Pi_o_0)
    elif t >= 0:
        return(m.cas_Pi_o[t] == (1-m.eta)*m.cas_Pi_o[t-1] +
               (3600/m.cas_m0)*(m.cmp_m[t] - m.exp_m[t]))
    else:
        return po.Constraint.Skip
    
def boundary(m,t):
    if t == min(m.T):
        return(m.cas_Pi_o[t] == m.cas_Pi_o_0)
    elif t == max(m.T):
        return(m.cas_Pi_o[t]== m.cas_Pi_o_0)
    else:
        return po.Constraint.Skip
    

