import numpy as np
import sympy as sym
from sympy.utilities.lambdify import lambdify
import matplotlib.pyplot as plt
from SFGcal import *
from nonzerobeta2 import nonzerob2

# parameters
inc1=60.*np.pi/180.0
inc2=55*np.pi/180.0
wl1=760.0
wl2=3500.0
n1SF=1.0
n2SF=1.36
n1VIS=1.0
n2VIS=1.36
n1IR=1.0
n2IR=1.36
incSF = incSFG(wl1,wl2,inc1,inc2)
print('incSF=',incSF*180.0/np.pi)
L = Lt(n1SF,n2SF,incSF)
L1 = Lt(n1VIS,n2VIS,inc1)
L2 = Lt(n1IR,n2IR,inc2)
ns = sym.symbols('ns')
chi2 = nonzerob2('chi2','civ')

# methyl group treated as C3v symmetry
beta2 = persym(nonzerob2('b2','c3v'))

# chi2eff of different polarization configuration
chi2ssp = chi2essp(L,L1,L2,incSF,inc1,inc2,chi2,beta2)
chi2sps = chi2esps(L,L1,L2,incSF,inc1,inc2,chi2,beta2)
chi2pss = chi2epss(L,L1,L2,incSF,inc1,inc2,chi2,beta2)
chi2ppp = chi2eppp(L,L1,L2,incSF,inc1,inc2,chi2,beta2)

# separate as and ss mode in CH3
asmodessp = chi2ssp.as_independent(beta2[0,0,2],beta2[2,2,2])[0]
asmodesps = chi2sps.as_independent(beta2[0,0,2],beta2[2,2,2])[0]
asmodepss = chi2pss.as_independent(beta2[0,0,2],beta2[2,2,2])[0]
asmodeppp = chi2ppp.as_independent(beta2[0,0,2],beta2[2,2,2])[0]

ssmodessp = chi2ssp.as_independent(beta2[0,0,2],beta2[2,2,2])[1]
ssmodesps = chi2sps.as_independent(beta2[0,0,2],beta2[2,2,2])[1]
ssmodepss = chi2pss.as_independent(beta2[0,0,2],beta2[2,2,2])[1]
ssmodeppp = chi2ppp.as_independent(beta2[0,0,2],beta2[2,2,2])[1]

# set the hyperpolarizability ratio beta2aac/beta2ccc = 3.4
expr0a = ((asmodessp.subs(beta2[0,2,0],1))**2).subs(theta,theta*sym.pi/180.0)
expr1a = ((asmodesps.subs(beta2[0,2,0],1))**2).subs(theta,theta*sym.pi/180.0)
expr2a = ((asmodepss.subs(beta2[0,2,0],1))**2).subs(theta,theta*sym.pi/180.0)
expr3a = ((asmodeppp.subs(beta2[0,2,0],1))**2).subs(theta,theta*sym.pi/180.0)

expr0s = ((ssmodessp.subs(beta2[0,0,2],3.4).subs(beta2[2,2,2],1))**2).subs(theta,theta*sym.pi/180.0)
expr1s = ((ssmodesps.subs(beta2[0,0,2],3.4).subs(beta2[2,2,2],1))**2).subs(theta,theta*sym.pi/180.0)
expr2s = ((ssmodepss.subs(beta2[0,0,2],3.4).subs(beta2[2,2,2],1))**2).subs(theta,theta*sym.pi/180.0)
expr3s = ((ssmodeppp.subs(beta2[0,0,2],3.4).subs(beta2[2,2,2],1))**2).subs(theta,theta*sym.pi/180.0)

# plot
xl=[theta]
exprla=[expr0a,expr1a,expr2a,expr3a]
exprls=[expr0s,expr1s,expr2s,expr3s]
xliml=[[0,90]]
plist=['ssp','sps','pss','ppp']
overlapexpr(xl, exprla, xliml, title=r'$CH_{3}-as$', labels=plist, xlabel=r'$\theta$', ylabel=r'$d^{2}R(\theta)$')
overlapexpr(xl, exprls, xliml, title=r'$CH_{3}-ss$', labels=plist, xlabel=r'$\theta$', ylabel=r'$d^{2}R(\theta)$')
plt.show()