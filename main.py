import sympy as sp
from IPython.display import display

x,y,z,t = sp.symbols("x y z t") # Variables

V,f,g,psi,Psi = sp.symbols("V f g psi Psi", cls = sp.Function) # Functions

k,c,hbar,m = sp.symbols("k c hbar m", constant = True) # Constants
# hbar = h/(2*sp.pi)

potential = V(x,t) #Time-Dependent Potential
wavefunc = Psi(x,t) #Time-Dependent Wave Function

SE = sp.Eq(sp.I*hbar*wavefunc.diff(t),hbar**2/(2*m) * wavefunc.diff(x,2) + potential*wavefunc) # Time-Dependent Schrodinger Equation

sp.init_printing(use_latex = True)

display(SE)