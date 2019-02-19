import sympy as sp
from sympy import I as i,E as e
from IPython.display import display
from matplotlib import pyplot as plt

x,y,z,t,r,theta,phi,k,omega = sp.symbols("x y z t r theta phi k omega") # Variables

V,f,g,psi,Psi = sp.symbols("V f g psi Psi", cls = sp.Function) # Functions

c,hbar,m,A,B,C,E = sp.symbols("c hbar m A B C E", constant = True) # Constants

m = 9.1093e-31
hbar = 1.0545718e-34
potential = V(x,t) #Time-Dependent Potential
wavefunc = Psi(x,t) #Time-Dependent Wave Function
tipotential = V(x) #Time-Independent Potential
tiwavefunc = psi(x) #Time-Independent Wave Function

potential = 1#sp.Piecewise((1000, sp.Abs(x) > C),(0, sp.Abs(x) <= C))
tipotential = 1#sp.Piecewise((1000, sp.Abs(x) > C),(0, sp.Abs(x) <= C))

SE = sp.Eq(i*hbar*wavefunc.diff(t),hbar**2/(2*m) * wavefunc.diff(x,2) + potential*wavefunc).simplify() # Time-Dependent Schrodinger Equation
tiSE = sp.Eq(E*tiwavefunc, hbar**2/(2*m) * tiwavefunc.diff(x,2)+ tipotential*tiwavefunc).simplify() # Time-Independent Schrodinger Equation


sp.init_printing(use_latex = True)
display(SE)
display(tiSE)
try:
    display(sp.pdsolve(SE,wavefunc))
    display(sp.dsolve(tiSE,tiwavefunc))
except:
    eq = sp.dsolve(tiSE,tiwavefunc)
    eq2 = eq.subs({"C2":0,"C1":1})
    c1 = 1/sp.integrate(eq2,(x,-sp.oo,sp.oo))
    display(eq.subs({"C2":0,"C1":c1}))