import math
import scipy
from scipy import linalg
from scipy import optimize
#import scipy.linalg
#import scipy.optimize
#import numpy

#class GATPFitPy:
    
def CpFit(atoms, rotors, linearity, t, cp):
    #the "main" function; eventually I will need to add routines for S, Hf, but these should be simple post-processing steps, as each should involve two equations and two unknowns
    #input: number of atoms, number of rotors, and linearity (0 for non-linear, 1 for linear), vector of temperatures, t (in Kelvin), corresponding vector of Cp (cal/mol-K)
    #output: NASA parameters for Cp/R, b1, b2, b3, b4, b5 (low temp parameters) and b6, b7, b8, b9, b10 (high temp parameters), and Tint
    R = 1.9872 #cal/mol-K
    B = 500. #K
    Tmin = 300. #K
    Tmax = 6000. #K
    Tintg = 1000. #K

    #convert from K to kK
    t = [x/1000. for x in t] 
    B = B/1000. 
    Tmin = Tmin/1000. 
    Tmax = Tmax/1000. 
    Tintg = Tintg/1000.

    cp = [x/R for x in cp] #convert to Cp/R
    (cp0, cpInf)=CpLimits(atoms, rotors, linearity)
    (a0, a1, a2, a3) = WilhoitFit(t, cp, cp0, cpInf, B)
    print [a0, a1, a2, a3]
    (b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, tint) = Wilhoit2NASA_TintOpt(cp0, cpInf, B, a0, a1, a2, a3, Tmin, Tmax, Tintg)
    #restore to conventional units of K for Tint and units based on K rather than kK in NASA polynomial coefficients
    tint=tint*1000.
    b2 = b2/1000.
    b7 = b7/1000.
    b3 = b3/1000000.
    b8 = b8/1000000.
    b4 = b4/1000000000.
    b9 = b9/1000000000.
    b5 = b5/1000000000000.
    b10= b10/1000000000000.

    return b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, tint

def CpFitFixedTint(atoms, rotors, linearity, t, cp):
    #the "main" function with Tint fixed; eventually I will need to add routines for S, Hf, but these should be simple post-processing steps, as each should involve two equations and two unknowns
    #input: number of atoms, number of rotors, and linearity (0 for non-linear, 1 for linear), vector of temperatures, t (in Kelvin), corresponding vector of Cp (cal/mol-K)
    #output: NASA parameters for Cp/R, b1, b2, b3, b4, b5 (low temp parameters) and b6, b7, b8, b9, b10 (high temp parameters), and Tint
    R = 1.9872 #cal/mol-K
    B = 500. #K
    Tmin = 300. #K
    Tmax = 6000. #K
    Tintg = 1000. #K

    #convert from K to kK
    t = [x/1000. for x in t] 
    B = B/1000. 
    Tmin = Tmin/1000. 
    Tmax = Tmax/1000. 
    Tintg = Tintg/1000.

    cp = [x/R for x in cp] #convert to Cp/R
    (cp0, cpInf)=CpLimits(atoms, rotors, linearity)
    (a0, a1, a2, a3) = WilhoitFit(t, cp, cp0, cpInf, B)
    print [a0, a1, a2, a3]
    #(b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, tint) = Wilhoit2NASA_TintOpt(cp0, cpInf, B, a0, a1, a2, a3, Tmin, Tmax, Tintg)
    (b1, b2, b3, b4, b5, b6, b7, b8, b9, b10) = Wilhoit2NASA(cp0, cpInf, B, a0, a1, a2, a3, Tmin, Tmax, Tintg)
    #restore to conventional units of K for Tint and units based on K rather than kK in NASA polynomial coefficients
    #tint=tint*1000.
    b2 = b2/1000.
    b7 = b7/1000.
    b3 = b3/1000000.
    b8 = b8/1000000.
    b4 = b4/1000000000.
    b9 = b9/1000000000.
    b5 = b5/1000000000000.
    b10= b10/1000000000000.

    #return b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, tint
    return b1, b2, b3, b4, b5, b6, b7, b8, b9, b10

def CpLimits(atoms, rotors, linearity):
    #(based off of lsfp_wilh1.f in GATPFit)
    #input: number of atoms, number of rotors, and linearity (0 for non-linear, 1 for linear)
    #output: Cp0/R, CpInf/R
    #monatomic case
    if(atoms == 1):
        cp0 = 2.5
        cpInf = 2.5
    #non-linear case
    elif(linearity == 0):
        cp0  = 4.0
        cpInf=3.*atoms-(2.+0.5*rotors)
        #linear case
    else:
        cp0  = 3.5
        cpInf=3.*atoms-1.5

    return cp0, cpInf



def WilhoitFit(t, cp, cp0, cpInf, B):
    #input: vector of temperatures, t (in kiloKelvin), corresponding vector of Cp/R, Cp0/R, CpInf/R, and B parameter in y=t/t+B (in kiloKelvin; will typically be 0.5)
    #output: fitted Wilhoit parameters, a0, a1, a2, a3

    m = len(t)

    #A = mx4
    #b = mx1
    #x = 4x1

    #from scipy import zeros
    A = scipy.zeros([m,4])
    b = scipy.zeros([m,1])
    for i in range(m):
        y = t[i]/(t[i]+B)
        A[i,0] = y*y*(y-1)
        A[i,1] = A[i,0]*y
        A[i,2] = A[i,1]*y
        A[i,3] = A[i,2]*y
        b[i] = (cp[i]-cp0)/(cpInf-cp0) - y*y
    #solve least squares problem A*x = b
    #http://docs.scipy.org/doc/scipy/reference/tutorial/linalg.html#solving-linear-least-squares-problems-and-pseudo-inverses
   # from linalg import lstsq
    x,resid,rank,sigma = linalg.lstsq(A,b)
    a0 = x[0]
    a1 = x[1]
    a2 = x[2]
    a3 = x[3]

    return a0, a1, a2, a3

def Wilhoit2NASA(cp0, cpInf, B, a0, a1, a2, a3, tmin, tmax, tint):
    #input: Wilhoit parameters, Cp0/R, CpInf/R, and B (kK), a0, a1, a2, a3, Tmin (minimum temperature (in kiloKelvin), Tmax (maximum temperature (in kiloKelvin), Tint (intermediate temperature, in kiloKelvin)
    #output: NASA parameters for Cp/R, b1, b2, b3, b4, b5 (low temp parameters) and b6, b7, b8, b9, b10 (high temp parameters)

    #construct 13*13 symmetric A matrix (in A*x = b); other elements will be zero
    #from scipy import zeros
    A = scipy.zeros([13,13])
    b = scipy.zeros([13,1])
    A[0,0] = 2*(tint - tmin)
    A[0,1] = tint*tint - tmin*tmin
    A[0,2] = 2.*(tint*tint*tint - tmin*tmin*tmin)/3
    A[0,3] = (tint*tint*tint*tint - tmin*tmin*tmin*tmin)/2
    A[0,4] = 2.*(tint*tint*tint*tint*tint - tmin*tmin*tmin*tmin*tmin)/5
    A[1,1] = A[0,2]
    A[1,2] = A[0,3]
    A[1,3] = A[0,4]
    A[1,4] = (tint*tint*tint*tint*tint*tint - tmin*tmin*tmin*tmin*tmin*tmin)/3
    A[2,2] = A[0,4]
    A[2,3] = A[1,4]
    A[2,4] = 2.*(tint*tint*tint*tint*tint*tint*tint - tmin*tmin*tmin*tmin*tmin*tmin*tmin)/7
    A[3,3] = A[2,4]
    A[3,4] = (tint*tint*tint*tint*tint*tint*tint*tint - tmin*tmin*tmin*tmin*tmin*tmin*tmin*tmin)/4
    A[4,4] = 2.*(tint*tint*tint*tint*tint*tint*tint*tint*tint - tmin*tmin*tmin*tmin*tmin*tmin*tmin*tmin*tmin)/9

    A[5,5] = 2*(tmax - tint)
    A[5,6] = tmax*tmax - tint*tint
    A[5,7] = 2.*(tmax*tmax*tmax - tint*tint*tint)/3
    A[5,8] = (tmax*tmax*tmax*tmax - tint*tint*tint*tint)/2
    A[5,9] = 2.*(tmax*tmax*tmax*tmax*tmax - tint*tint*tint*tint*tint)/5
    A[6,6] = A[5,7]
    A[6,7] = A[5,8]
    A[6,8] = A[5,9]
    A[6,9] = (tmax*tmax*tmax*tmax*tmax*tmax - tint*tint*tint*tint*tint*tint)/3
    A[7,7] = A[5,9]
    A[7,8] = A[6,9]
    A[7,9] = 2.*(tmax*tmax*tmax*tmax*tmax*tmax*tmax - tint*tint*tint*tint*tint*tint*tint)/7
    A[8,8] = A[7,9]
    A[8,9] = (tmax*tmax*tmax*tmax*tmax*tmax*tmax*tmax - tint*tint*tint*tint*tint*tint*tint*tint)/4
    A[9,9] = 2.*(tmax*tmax*tmax*tmax*tmax*tmax*tmax*tmax*tmax - tint*tint*tint*tint*tint*tint*tint*tint*tint)/9

    A[0,10] = 1.
    A[1,10] = tint
    A[1,11] = 1.
    A[2,10] = tint*tint
    A[2,11] = 2*tint
    A[2,12] = 2.
    A[3,10] = A[2,10]*tint
    A[3,11] = 3*A[2,10]
    A[3,12] = 6*tint
    A[4,10] = A[3,10]*tint
    A[4,11] = 4*A[3,10]
    A[4,12] = 12*A[2,10]

    A[5,10] = -A[0,10]
    A[6,10] = -A[1,10]
    A[6,11] = -A[1,11]
    A[7,10] = -A[2,10]
    A[7,11] = -A[2,11]
    A[7,12] = -A[2,12]
    A[8,10] = -A[3,10]
    A[8,11] = -A[3,11]
    A[8,12] = -A[3,12]
    A[9,10] = -A[4,10]
    A[9,11] = -A[4,11]
    A[9,12] = -A[4,12]

    # make the matrix symmetric
    for i in range(1,13):
        for j in range(0, i):
            A[i,j] = A[j,i]

    #construct b vector
    #store values at tint (this will avoid evaluating them twice)
    w0int = WilhoitInt0(cp0,cpInf,B,a0,a1,a2,a3,tint)
    w1int = WilhoitInt1(cp0,cpInf,B,a0,a1,a2,a3,tint)
    w2int = WilhoitInt2(cp0,cpInf,B,a0,a1,a2,a3,tint)
    w3int = WilhoitInt3(cp0,cpInf,B,a0,a1,a2,a3,tint)
    w4int = WilhoitInt4(cp0,cpInf,B,a0,a1,a2,a3,tint)

    b[0] = 2*(w0int - WilhoitInt0(cp0,cpInf,B,a0,a1,a2,a3,tmin))
    b[1] = 2*(w1int - WilhoitInt1(cp0,cpInf,B,a0,a1,a2,a3,tmin))
    b[2] = 2*(w2int - WilhoitInt2(cp0,cpInf,B,a0,a1,a2,a3,tmin))
    b[3] = 2*(w3int - WilhoitInt3(cp0,cpInf,B,a0,a1,a2,a3,tmin))
    b[4] = 2*(w4int - WilhoitInt4(cp0,cpInf,B,a0,a1,a2,a3,tmin))
    b[5] = 2*(WilhoitInt0(cp0,cpInf,B,a0,a1,a2,a3,tmax) - w0int)
    b[6] = 2*(WilhoitInt1(cp0,cpInf,B,a0,a1,a2,a3,tmax) - w1int)
    b[7] = 2*(WilhoitInt2(cp0,cpInf,B,a0,a1,a2,a3,tmax) - w2int)
    b[8] = 2*(WilhoitInt3(cp0,cpInf,B,a0,a1,a2,a3,tmax) - w3int)
    b[9] = 2*(WilhoitInt4(cp0,cpInf,B,a0,a1,a2,a3,tmax) - w4int)
   # b[10] = 0
   # b[11] = 0
   # b[12] = 0

    #solve A*x=b for x (note that factor of 2 in b vector and 10*10 submatrix of A matrix is not required; not including it should give same result, except Lagrange multipliers will differ by a factor of two)
    #from linalg import solve
    #print A
    x = linalg.solve(A,b,overwrite_a=1,overwrite_b=1)
    b1 = x[0]
    b2 = x[1]
    b3 = x[2]
    b4 = x[3]
    b5 = x[4]
    b6 = x[5]
    b7 = x[6]
    b8 = x[7]
    b9 = x[8]
    b10 = x[9]
    #return x[0:10]
    return b1, b2, b3, b4, b5, b6, b7, b8, b9, b10
    
def Wilhoit2NASA_TintOpt(cp0, cpInf, B, a0, a1, a2, a3, tmin, tmax, tintg):
    #input: Wilhoit parameters, Cp0/R, CpInf/R, and B (kK), a0, a1, a2, a3, Tmin (minimum temperature (in kiloKelvin), Tmax (maximum temperature (in kiloKelvin), Tintg (guess intermediate temperature, in kiloKelvin)
    #output: NASA parameters for Cp/R, b1, b2, b3, b4, b5 (low temp parameters) and b6, b7, b8, b9, b10 (high temp parameters), and Tint
    #1. vary Tint, using Tintg as a starting guess, to minimize TintOpt_objFun
    #from optimize import fminbound
    #tint = optimize.fminbound(TintOpt_objFun, tmin+0.2, tmax-0.5, args=(cp0, cpInf,B,a0,a1,a2,a3,tmin,tmax))
    tint = optimize.fminbound(TintOpt_objFun, tmin, tmax, args=(cp0, cpInf,B,a0,a1,a2,a3,tmin,tmax))
    #note that we have not used the specified guess, tintg when using this minimization routine
    #2. determine the bi parameters based on the optimized Tint (alternatively, maybe we could have TintOpt_objFun also return these parameters, along with the objective function, which would avoid an extra calculation)
    (b1, b2, b3, b4, b5, b6, b7, b8, b9, b10) = Wilhoit2NASA(cp0,cpInf,B,a0,a1,a2,a3,tmin,tmax,tint)
    #bvec = Wilhoit2NASA(cp0,cpInf,B,a0,a1,a2,a3,tmin,tmax,tint)
    return b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, tint
                                                                                                    
def TintOpt_objFun(tint, cp0, cpInf, B, a0, a1, a2, a3, tmin, tmax):
	#input: Tint (intermediate temperature, in kiloKelvin); Wilhoit parameters, Cp0/R, CpInf/R, and B (kK), a0, a1, a2, a3, Tmin (minimum temperature (in kiloKelvin), Tmax (maximum temperature (in kiloKelvin)
	#output: the quantity Integrate[(Cp(Wilhoit)/R-Cp(NASA)/R)^2, {t, tmin, tmax}]
	(b1, b2, b3, b4, b5, b6, b7, b8, b9, b10) = Wilhoit2NASA(cp0,cpInf,B,a0,a1,a2,a3,tmin,tmax,tint)
	result = (Wilhoit2Int(cp0,cpInf,B,a0,a1,a2,a3,tmax) - Wilhoit2Int(cp0,cpInf,B,a0,a1,a2,a3,tmin) +
                 NASA2Int(b1,b2,b3,b4,b5,tint)-NASA2Int(b1,b2,b3,b4,b5,tmin) + NASA2Int(b6,b7,b8,b9,b10,tmax) - NASA2Int(b6,b7,b8,b9,b10,tint)
                   - 2* (b6*WilhoitInt0(cp0,cpInf,B,a0,a1,a2,a3,tmax)+(b1-b6)*WilhoitInt0(cp0,cpInf,B,a0,a1,a2,a3,tint) - b1*WilhoitInt0(cp0,cpInf,B,a0,a1,a2,a3,tmin)
                 +b7*WilhoitInt1(cp0,cpInf,B,a0,a1,a2,a3,tmax)+(b2-b7)*WilhoitInt1(cp0,cpInf,B,a0,a1,a2,a3,tint) - b2*WilhoitInt1(cp0,cpInf,B,a0,a1,a2,a3,tmin)
                 +b8*WilhoitInt2(cp0,cpInf,B,a0,a1,a2,a3,tmax)+(b3-b8)*WilhoitInt2(cp0,cpInf,B,a0,a1,a2,a3,tint) - b3*WilhoitInt2(cp0,cpInf,B,a0,a1,a2,a3,tmin)
                 +b9*WilhoitInt3(cp0,cpInf,B,a0,a1,a2,a3,tmax)+(b4-b9)*WilhoitInt3(cp0,cpInf,B,a0,a1,a2,a3,tint) - b4*WilhoitInt3(cp0,cpInf,B,a0,a1,a2,a3,tmin)
                 +b10*WilhoitInt4(cp0,cpInf,B,a0,a1,a2,a3,tmax)+(b5-b10)*WilhoitInt4(cp0,cpInf,B,a0,a1,a2,a3,tint) - b5*WilhoitInt4(cp0,cpInf,B,a0,a1,a2,a3,tmin)))
	print tint
	return result


#analytical integrals:

#input (for all functions WilhoitXIntN): Wilhoit parameters: Cp0(/R), CpInf(/R), B, a0, a1, a2, a3 and t (in kiloKelvin)

def WilhoitInt0(cp0, cpInf, B, a0, a1, a2, a3, t):
    #output: the quantity Integrate[Cp(Wilhoit)/R, t'] evaluated at t'=t
    result =        ( cpInf*t + (a3*B**6*(cp0 - cpInf))/(5.*(B + t)**5) - ((a2 + 5*a3)*B**5*(cp0 - cpInf))/(4.*(B + t)**4) + ((a1 + 4*a2 + 10*a3)*B**4*(cp0 - cpInf))/(3.*(B + t)**3) - 
           ((a0 + 3*a1 + 6*a2 + 10*a3)*B**3*(cp0 - cpInf))/(2.*(B + t)**2) + ((1 + 2*a0 + 3*a1 + 4*a2 + 5*a3)*B**2*(cp0 - cpInf))/(B + t) + (2 + a0 + a1 + a2 + a3)*B*(cp0 - cpInf)*math.log(B + t))
    return result

def WilhoitInt1(cp0, cpInf, B, a0, a1, a2, a3, t):
    #output: the quantity Integrate[Cp(Wilhoit)/R*t, t'] evaluated at t'=t
    result =       ( (2 + a0 + a1 + a2 + a3)*B*(cp0 - cpInf)*t + (cpInf*t**2)/2. + (a3*B**7*(-cp0 + cpInf))/(5.*(B + t)**5) + ((a2 + 6*a3)*B**6*(cp0 - cpInf))/(4.*(B + t)**4) - 
           ((a1 + 5*(a2 + 3*a3))*B**5*(cp0 - cpInf))/(3.*(B + t)**3) + ((a0 + 4*a1 + 10*(a2 + 2*a3))*B**4*(cp0 - cpInf))/(2.*(B + t)**2) - 
           ((1 + 3*a0 + 6*a1 + 10*a2 + 15*a3)*B**3*(cp0 - cpInf))/(B + t) - (3 + 3*a0 + 4*a1 + 5*a2 + 6*a3)*B**2*(cp0 - cpInf)*math.log(B + t))
    return result

def WilhoitInt2(cp0, cpInf, B, a0, a1, a2, a3, t):
    #output: the quantity Integrate[Cp(Wilhoit)/R*t^2, t'] evaluated at t'=t
    result =       (  -((3 + 3*a0 + 4*a1 + 5*a2 + 6*a3)*B**2*(cp0 - cpInf)*t) + ((2 + a0 + a1 + a2 + a3)*B*(cp0 - cpInf)*t**2)/2. + (cpInf*t**3)/3. + (a3*B**8*(cp0 - cpInf))/(5.*(B + t)**5) - 
           ((a2 + 7*a3)*B**7*(cp0 - cpInf))/(4.*(B + t)**4) + ((a1 + 6*a2 + 21*a3)*B**6*(cp0 - cpInf))/(3.*(B + t)**3) - ((a0 + 5*(a1 + 3*a2 + 7*a3))*B**5*(cp0 - cpInf))/(2.*(B + t)**2) + 
           ((1 + 4*a0 + 10*a1 + 20*a2 + 35*a3)*B**4*(cp0 - cpInf))/(B + t) + (4 + 6*a0 + 10*a1 + 15*a2 + 21*a3)*B**3*(cp0 - cpInf)*math.log(B + t))
    return result

def WilhoitInt3(cp0, cpInf, B, a0, a1, a2, a3, t):
    #output: the quantity Integrate[Cp(Wilhoit)/R*t^3, t'] evaluated at t'=t
    result =       ( (4 + 6*a0 + 10*a1 + 15*a2 + 21*a3)*B**3*(cp0 - cpInf)*t + ((3 + 3*a0 + 4*a1 + 5*a2 + 6*a3)*B**2*(-cp0 + cpInf)*t**2)/2. + ((2 + a0 + a1 + a2 + a3)*B*(cp0 - cpInf)*t**3)/3. + 
           (cpInf*t**4)/4. + (a3*B**9*(-cp0 + cpInf))/(5.*(B + t)**5) + ((a2 + 8*a3)*B**8*(cp0 - cpInf))/(4.*(B + t)**4) - ((a1 + 7*(a2 + 4*a3))*B**7*(cp0 - cpInf))/(3.*(B + t)**3) + 
           ((a0 + 6*a1 + 21*a2 + 56*a3)*B**6*(cp0 - cpInf))/(2.*(B + t)**2) - ((1 + 5*a0 + 15*a1 + 35*a2 + 70*a3)*B**5*(cp0 - cpInf))/(B + t) - 
           (5 + 10*a0 + 20*a1 + 35*a2 + 56*a3)*B**4*(cp0 - cpInf)*math.log(B + t))
    return result

def WilhoitInt4(cp0, cpInf, B, a0, a1, a2, a3, t):
    #output: the quantity Integrate[Cp(Wilhoit)/R*t^4, t'] evaluated at t'=t
    result=       ( -((5 + 10*a0 + 20*a1 + 35*a2 + 56*a3)*B**4*(cp0 - cpInf)*t) + ((4 + 6*a0 + 10*a1 + 15*a2 + 21*a3)*B**3*(cp0 - cpInf)*t**2)/2. + 
           ((3 + 3*a0 + 4*a1 + 5*a2 + 6*a3)*B**2*(-cp0 + cpInf)*t**3)/3. + ((2 + a0 + a1 + a2 + a3)*B*(cp0 - cpInf)*t**4)/4. + (cpInf*t**5)/5. + (a3*B**10*(cp0 - cpInf))/(5.*(B + t)**5) - 
           ((a2 + 9*a3)*B**9*(cp0 - cpInf))/(4.*(B + t)**4) + ((a1 + 8*a2 + 36*a3)*B**8*(cp0 - cpInf))/(3.*(B + t)**3) - ((a0 + 7*(a1 + 4*(a2 + 3*a3)))*B**7*(cp0 - cpInf))/(2.*(B + t)**2) + 
           ((1 + 6*a0 + 21*a1 + 56*a2 + 126*a3)*B**6*(cp0 - cpInf))/(B + t) + (6 + 15*a0 + 35*a1 + 70*a2 + 126*a3)*B**5*(cp0 - cpInf)*math.log(B + t))
    return result

def Wilhoit2Int(cp0, cpInf, B, a0, a1, a2, a3, t):
    #output: the quantity Integrate[(Cp(Wilhoit)/R)^2, t'] evaluated at t'=t
    result =         (cpInf**2*t - (a3**2*B**12*(cp0 - cpInf)**2)/(11.*(B + t)**11) + (a3*(a2 + 5*a3)*B**11*(cp0 - cpInf)**2)/(5.*(B + t)**10) - 
           ((a2**2 + 18*a2*a3 + a3*(2*a1 + 45*a3))*B**10*(cp0 - cpInf)**2)/(9.*(B + t)**9) + ((4*a2**2 + 36*a2*a3 + a1*(a2 + 8*a3) + a3*(a0 + 60*a3))*B**9*(cp0 - cpInf)**2)/(4.*(B + t)**8) - 
           ((a1**2 + 14*a1*(a2 + 4*a3) + 2*(14*a2**2 + a3 + 84*a2*a3 + 105*a3**2 + a0*(a2 + 7*a3)))*B**8*(cp0 - cpInf)**2)/(7.*(B + t)**7) + 
           ((3*a1**2 + a2 + 28*a2**2 + 7*a3 + 126*a2*a3 + 126*a3**2 + 7*a1*(3*a2 + 8*a3) + a0*(a1 + 6*a2 + 21*a3))*B**7*(cp0 - cpInf)**2)/(3.*(B + t)**6) - 
           (B**6*(cp0 - cpInf)*(a0**2*(cp0 - cpInf) + 15*a1**2*(cp0 - cpInf) + 10*a0*(a1 + 3*a2 + 7*a3)*(cp0 - cpInf) + 2*a1*(1 + 35*a2 + 70*a3)*(cp0 - cpInf) + 
                2*(35*a2**2*(cp0 - cpInf) + 6*a2*(1 + 21*a3)*(cp0 - cpInf) + a3*(5*(4 + 21*a3)*cp0 - 21*(cpInf + 5*a3*cpInf)))))/(5.*(B + t)**5) + 
           (B**5*(cp0 - cpInf)*(14*a2*cp0 + 28*a2**2*cp0 + 30*a3*cp0 + 84*a2*a3*cp0 + 60*a3**2*cp0 + 2*a0**2*(cp0 - cpInf) + 10*a1**2*(cp0 - cpInf) + 
                a0*(1 + 10*a1 + 20*a2 + 35*a3)*(cp0 - cpInf) + a1*(5 + 35*a2 + 56*a3)*(cp0 - cpInf) - 15*a2*cpInf - 28*a2**2*cpInf - 35*a3*cpInf - 84*a2*a3*cpInf - 60*a3**2*cpInf))/
            (2.*(B + t)**4) - (B**4*(cp0 - cpInf)*((1 + 6*a0**2 + 15*a1**2 + 32*a2 + 28*a2**2 + 50*a3 + 72*a2*a3 + 45*a3**2 + 2*a1*(9 + 21*a2 + 28*a3) + a0*(8 + 20*a1 + 30*a2 + 42*a3))*cp0 - 
                (1 + 6*a0**2 + 15*a1**2 + 40*a2 + 28*a2**2 + 70*a3 + 72*a2*a3 + 45*a3**2 + a0*(8 + 20*a1 + 30*a2 + 42*a3) + a1*(20 + 42*a2 + 56*a3))*cpInf))/(3.*(B + t)**3) + 
           (B**3*(cp0 - cpInf)*((2 + 2*a0**2 + 3*a1**2 + 9*a2 + 4*a2**2 + 11*a3 + 9*a2*a3 + 5*a3**2 + a0*(5 + 5*a1 + 6*a2 + 7*a3) + a1*(7 + 7*a2 + 8*a3))*cp0 - 
                (2 + 2*a0**2 + 3*a1**2 + 15*a2 + 4*a2**2 + 21*a3 + 9*a2*a3 + 5*a3**2 + a0*(6 + 5*a1 + 6*a2 + 7*a3) + a1*(10 + 7*a2 + 8*a3))*cpInf))/(B + t)**2 - 
           (B**2*((2 + a0 + a1 + a2 + a3)**2*cp0**2 - 2*(5 + a0**2 + a1**2 + 8*a2 + a2**2 + 9*a3 + 2*a2*a3 + a3**2 + 2*a0*(3 + a1 + a2 + a3) + a1*(7 + 2*a2 + 2*a3))*cp0*cpInf + 
                (6 + a0**2 + a1**2 + 12*a2 + a2**2 + 14*a3 + 2*a2*a3 + a3**2 + 2*a1*(5 + a2 + a3) + 2*a0*(4 + a1 + a2 + a3))*cpInf**2))/(B + t) + 
           2*(2 + a0 + a1 + a2 + a3)*B*(cp0 - cpInf)*cpInf*math.log(B + t))
    return result

def NASA2Int(c1,c2,c3,c4,c5,t) :
    #input: NASA parameters for Cp/R, c1, c2, c3, c4, c5 (either low or high temp parameters), temperature t (in kiloKelvin; an endpoint of the low or high temp range
    #output: the quantity Integrate[(Cp(NASA)/R)^2, t'] evaluated at t'=t 
    #can speed further by precomputing and storing e.g. thigh^2, tlow^2, etc.
    result = c1*c1*t + c1*c2*t*t + (2*c1*c3+c2*c2)/3*t*t*t + (c1*c4+c2*c3)/2*t*t*t*t + (2*c1*c5 + 2*c2*c4 + c3*c3)/5*t*t*t*t*t + (c2*c5 + c3*c4)/3*t*t*t*t*t*t + (2*c3*c5 + c4*c4)/7*t*t*t*t*t*t*t + c4*c5/4*t*t*t*t*t*t*t*t + c5*c5/9*t*t*t*t*t*t*t*t*t
    return result


