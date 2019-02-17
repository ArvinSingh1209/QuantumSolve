import sympy as sp
from sympy import I as i,E as e
from IPython.display import display

x,y,z,t,r,theta,phi,k,omega = sp.symbols("x y z t r theta phi k omega") # Variables

V,f,g,psi,Psi,E,p = sp.symbols("V f g psi Psi E p", cls = sp.Function) # Functions

c,hbar,m,A,B,C = sp.symbols("c hbar m A B C", constant = True) # Constants

E = hbar*omega # Energy
p = hbar*k # Momentum

potential = V(x,t) #Time-Dependent Potential
wavefunc = Psi(x,t) #Time-Dependent Wave Function

SE = sp.Eq(i*hbar*wavefunc.diff(t),hbar**2/(2*m) * wavefunc.diff(x,2) + potential*wavefunc).simplify() # Time-Dependent Schrodinger Equation

sp.init_printing(use_latex = True)
display(SE)

generalsol = sp.Eq(wavefunc,A*e**(-i*(p*x-E*t)/hbar)).simplify()
display(generalsol)

potential = None