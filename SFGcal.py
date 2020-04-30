import numpy as np
import sympy as sym
from sympy.utilities.lambdify import lambdify
import matplotlib.pyplot as plt
from nonzerobeta2 import nonzerob2

#list of optical field / 0 SF 1 VIS 2 IR

flst = ['SF', 'VIS', 'IR']
theta = sym.symbols('theta')
phi = sym.symbols('phi')
psi = sym.symbols('psi')


# incident angle of SF
def incSFG(wlVIS,wlIR,incVIS,incIR):
    kVIS = 2.*np.pi/wlVIS
    kIR = 2.*np.pi/wlIR
    return np.arcsin((kVIS*np.sin(incVIS)+kIR*np.sin(incIR))/(kVIS+kIR))

# refractive index of the interfacial layer (calculated with the Lorentz model)
def nplm(n2):
    return np.sqrt((n2**2*(n2**2+5))/(4*n2**2+2))

# refracted angle ref
def refa(n1,n2,inc):
    return np.arcsin(n1*np.sin(inc)/n2)

# wavelength (nm) to frequency (Hz)
def wl2f(wl):
    return 299792458/(wl*10**-9)

# coefficient A
def coeffA(fSF,incSF,n1SF,n1VIS,n1IR):
    return (8*np.pi**3*fSF**2*(1/np.cos(incSF)**2))/(299792458**3*n1SF*n1VIS*n1IR)

# Fresnel cofficient
def Lt(n1,n2,inc):
    L = sym.MutableDenseNDimArray(sym.zeros(3),(3,3))
    ref = refa(n1,n2,inc)
    np = nplm(n2)
    L[0,0] = (2*n1*sym.cos(ref))/(n1*sym.cos(ref)+n2*sym.cos(inc))
    L[1,1] = (2*n1*sym.cos(inc))/(n1*sym.cos(inc)+n2*sym.cos(ref))
    L[2,2] = (2*n2*sym.cos(inc))/(n1*sym.cos(ref)+n2*sym.cos(inc))*(n1/np)**2
    return L

# input field unit polarization vector ei | pa for polarization angle a for incidence angle
def ei(pa,inc):
    e = sym.MutableDenseNDimArray([0,0,0])
    e[0] = sym.cos(pa)*sym.cos(inc)
    e[1] = sym.sin(pa)
    e[2] = sym.cos(pa)*sym.sin(inc)
    return e

# output field eo polarization vector 
def eo(pa,inc):
    e = sym.MutableDenseNDimArray([0,0,0])
    e[0] = -sym.cos(pa)*sym.cos(inc)
    e[1] = sym.sin(pa)
    e[2] = sym.cos(pa)*sym.sin(inc)
    return e

# Permutation Symmetry of chi2/beta2 tensor (sym tensor) elements in SFG with both the SF and VIS frequencies off resonance
def persym(tensor):
    for k in range(3):
        for i in range(3):
            for j in range(3):
                if i > j:
                    tensor[i,j,k] = tensor[j,i,k]
    return tensor

# Euler rotational transformation matrix
def rtm(theta,phi,psi):
    m = sym.MutableDenseNDimArray(np.zeros(9),(3,3))
    m[0,0] = sym.cos(psi)*sym.cos(phi)-sym.cos(theta)*sym.sin(phi)*sym.sin(psi)
    m[0,1] = -sym.sin(psi)*sym.cos(phi)-sym.cos(theta)*sym.sin(phi)*sym.cos(psi)
    m[0,2] = sym.sin(theta)*sym.sin(phi)
    m[1,0] = sym.cos(psi)*sym.sin(phi)+sym.cos(theta)*sym.cos(phi)*sym.sin(psi)
    m[1,1] = -sym.sin(psi)*sym.sin(phi)+sym.cos(theta)*sym.cos(phi)*sym.cos(psi)
    m[1,2] = -sym.sin(theta)*sym.cos(phi)
    m[2,0] = sym.sin(theta)*sym.sin(psi)
    m[2,1] = sym.sin(theta)*sym.cos(psi)
    m[2,2] = sym.cos(theta)
    return m

# the average over different molecular orientations a is the function of (theta, phi, psi) dis is the distribution function
def oriave(a,dis):
    def f(theta, phi, psi):
        return dis*sym.sin(theta)
    def g(theta, phi, psi):
        return a*dis*sym.sin(theta)
    d = sym.integrate(g, (theta, 0, sym.pi), (phi, 0, 2*sym.pi), (psi, 0, 2*sym.pi))
    n = sym.integrate(f, (theta, 0, sym.pi), (phi, 0, 2*sym.pi), (psi, 0, 2*sym.pi))
    ave = n / d
    return ave

'''
# molecular hyperpolarizability beta2 lab and molecular coordinates conversion integrate over phi and psi (uniform distribution) theta (delta distribution)
def mh2chi2(beta2):
    # ns = sym.symbols('ns')
    chi2 = sym.MutableDenseNDimArray(sym.zeros(27),(3,3,3))
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for l in range(3):
                    for m in range(3):
                        for n in range(3):                        
                            dist = rtm(theta,phi,psi)[i,l]*rtm(theta,phi,psi)[j,m]*rtm(theta,phi,psi)[k,n]*1/(sym.pi*2)*1/(sym.pi*2)
                            chi2[i,j,k] += sym.integrate(dist, (phi, 0, 2*sym.pi),(psi, 0, 2*sym.pi))*beta2[l,m,n]
                chi2[i,j,k].simplify()
    # chi2 *= ns
    return chi2
'''

# calculate specific chi2 tensor element chi2eff*1/Ns
def mh2chi2e(beta2t,x,y,z):
    #ns = sym.symbols('ns')
    chi2xyz = 0
    for l in range(3):
        for m in range(3):
            for n in range(3):                        
                if beta2t[l,m,n] == 0:
                    continue
                else:
                    dist = rtm(theta,phi,psi)[x,l]*rtm(theta,phi,psi)[y,m]*rtm(theta,phi,psi)[z,n]*1/(sym.pi*2)*1/(sym.pi*2)
                    chi2xyz += (sym.integrate(dist, (phi, 0, 2*sym.pi),(psi, 0, 2*sym.pi))*beta2t[l,m,n]).simplify()
    #chi2xyz *= ns
    return chi2xyz

# the effective nonlinear susceptibility chi2eff interfacial unit field = L*ei(o) tensordot
def chi2eff(L,L1,L2,inc,inc1,inc2,pa,pa1,pa2,chi2):
    chi2e = 0
    a = sym.tensorcontraction(sym.tensorproduct(L,eo(pa,inc)),(1,2))
    b = sym.tensorcontraction(sym.tensorproduct(L1,ei(pa1,inc1)),(1,2))
    c = sym.tensorcontraction(sym.tensorproduct(L2,ei(pa2,inc2)),(1,2))
    for i in range(3):
        for j in range(3):
            for k in range(3):
                if (a[i] and b[j] and c[k] and chi2[i,j,k]) == 0:
                    continue
                else:
                    chi2e += a[i]*b[j]*c[k]*chi2[i,j,k]
    return chi2e

# chi2eff ssp polarization configuration
def chi2essp(L,L1,L2,inc,inc1,inc2,chi2,beta2):
    chi2e = 0
    a = sym.tensorcontraction(sym.tensorproduct(L,eo(sym.pi/2,inc)),(1,2))
    b = sym.tensorcontraction(sym.tensorproduct(L1,ei(sym.pi/2,inc1)),(1,2))
    c = sym.tensorcontraction(sym.tensorproduct(L2,ei(0,inc2)),(1,2))
    for i in range(3):
        for j in range(3):
            for k in range(3):
                if (a[i] and b[j] and c[k] and chi2[i,j,k]) == 0:
                    continue
                else:
                    chi2e += a[i]*b[j]*c[k]*chi2[i,j,k]
    chi2xxz = mh2chi2e(beta2,0,0,2)
    chi2e = chi2e.subs(chi2[0,0,2],chi2xxz)
    return chi2e

# chi2eff sps polarization configuration
def chi2esps(L,L1,L2,inc,inc1,inc2,chi2,beta2):
    chi2e = 0
    a = sym.tensorcontraction(sym.tensorproduct(L,eo(sym.pi/2,inc)),(1,2))
    b = sym.tensorcontraction(sym.tensorproduct(L1,ei(0,inc1)),(1,2))
    c = sym.tensorcontraction(sym.tensorproduct(L2,ei(sym.pi/2,inc2)),(1,2))
    for i in range(3):
        for j in range(3):
            for k in range(3):
                if (a[i] and b[j] and c[k] and chi2[i,j,k]) == 0:
                    continue
                else:
                    chi2e += a[i]*b[j]*c[k]*chi2[i,j,k]
    chi2xzx = mh2chi2e(beta2,0,2,0)
    chi2e = chi2e.subs(chi2[0,2,0],chi2xzx)    
    return chi2e

# chi2eff pss polarization configuration
def chi2epss(L,L1,L2,inc,inc1,inc2,chi2,beta2):
    chi2e = 0
    a = sym.tensorcontraction(sym.tensorproduct(L,eo(0,inc)),(1,2))
    b = sym.tensorcontraction(sym.tensorproduct(L1,ei(sym.pi/2,inc1)),(1,2))
    c = sym.tensorcontraction(sym.tensorproduct(L2,ei(sym.pi/2,inc2)),(1,2))
    for i in range(3):
        for j in range(3):
            for k in range(3):
                if (a[i] and b[j] and c[k] and chi2[i,j,k]) == 0:
                    continue
                else:
                    chi2e += a[i]*b[j]*c[k]*chi2[i,j,k]
    chi2zxx = mh2chi2e(beta2,2,0,0)
    chi2e = chi2e.subs(chi2[2,0,0],chi2zxx)
    return chi2e

# chi2eff ppp polarization configuration
def chi2eppp(L,L1,L2,inc,inc1,inc2,chi2,beta2):
    chi2e = 0
    a = sym.tensorcontraction(sym.tensorproduct(L,eo(0,inc)),(1,2))
    b = sym.tensorcontraction(sym.tensorproduct(L1,ei(0,inc1)),(1,2))
    c = sym.tensorcontraction(sym.tensorproduct(L2,ei(0,inc2)),(1,2))
    for i in range(3):
        for j in range(3):
            for k in range(3):
                if (a[i] and b[j] and c[k] and chi2[i,j,k]) == 0:
                    continue
                else:
                    chi2e += a[i]*b[j]*c[k]*chi2[i,j,k]
    chi2xxz = mh2chi2e(beta2,0,0,2)
    chi2xzx = mh2chi2e(beta2,0,2,0)
    chi2zxx = mh2chi2e(beta2,2,0,0)
    chi2zzz = mh2chi2e(beta2,2,2,2)
    chi2e = chi2e.subs(chi2[0,0,2],chi2xxz).subs(chi2[0,2,0],chi2xzx).subs(chi2[2,0,0],chi2zxx).subs(chi2[2,2,2],chi2zzz)
    return chi2e

# plot a sympy expression (function of x(sym symbols)) with matplotlib
def plotsexpr(x, expr, xlim, title=None, label=None, xlabel=None, ylabel=None):
    func = lambdify(x,expr,'numpy')
    a = np.linspace(xlim[0],xlim[1],500)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(a,func(a),label=label)
    ax.set_xlim(xlim)
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)    

# overlap curves of several sympy expressions
def overlapexpr(xlist, exprlist, xlimlist, title=None, labels=None, xlabel=None, ylabel=None):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    xmin = []
    xmax = []
    if len(xlist) == 1:
        xlist = xlist*len(exprlist)
    if len(xlimlist) == 1:
        xlimlist = xlimlist*len(exprlist)
    if labels is None:
        labels = []
        for expr in exprlist:
            labels.append(str(expr))
    for i,j,k,l in zip(xlist,exprlist,xlimlist,labels):
        xmin.append(k[0])
        xmax.append(k[1])
        func = lambdify(i,j,'numpy')
        a = np.linspace(k[0],k[1],500)
        ax.plot(a,func(a),label=l)        
    ax.set_xlim(min(xmin),max(xmax))
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.legend(loc='best')
