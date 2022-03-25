import numpy as np
from sklearn.linear_model import LinearRegression
import random
import math
import statistics
import scipy.stats as st
from scipy.stats import binom


lowerb = 1
upperbd = 30
x1 = 50
x2 = 2
R = 1
sig1 = []
sig2 = []
sig3 = []
for i in range(400):
    t1 = random.uniform(lowerb,upperbd)
    t2 = random.uniform(lowerb,upperbd)
    d = t2/t1
    rho = R/d
    c1 = (x1-x2*rho)/(math.sqrt(x1 +x2*rho**2))
    c2 = (x1 -x2*rho)/(math.sqrt((x1+x2)*rho))
    c3 = (np.log(x1/x2) - np.log(rho))/(math.sqrt(1/x1 + 1/x2))
    p1 = 1-st.norm.cdf(c1)
    p2 = 1-st.norm.cdf(c2)
    p3 = 1-st.norm.cdf(c3)
    sig1.append(p1)
    sig2.append(p2)
    sig3.append(p3)
cpval = []
for i in range(400): 
    t1 = random.uniform(lowerb,upperbd)
    t2 = random.uniform(lowerb,upperbd)
    d = t1/t2
    rho = R/d
    q = rho/(1+rho)
    c_val = binom.cdf(x1, x1+x2, q)
    p_val = 1-c_val
    actual_p_val = 2*min(c_val, p_val)
    cpval.append(actual_p_val)
print(statistics.mean(sig1))
print(statistics.mean(sig2))
print(statistics.mean(sig3))
print(statistics.mean(cpval))
print(statistics.variance(cpval))
