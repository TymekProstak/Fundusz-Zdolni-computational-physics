import numpy as np
import math
import cmath
import matplotlib.pyplot as plt
import copy
from scipy import sparse
from scipy import linalg
from scipy.sparse.linalg import spsolve
from numpy import save
Ny=401   #ilość punktów na siatce
Nx=401

xmin=-12
xmax=12
ymin=-12
ymax=12

X=np.linspace(xmin,xmax,Nx)
Y=np.linspace(ymin,ymax,Ny)

dx=abs(X[0]-X[1])
dy=abs(Y[0]-Y[1])
print(dx)
print(dy)
e=10**(-8)
dt=1j*0.01

n=1000
g=0.015
K=10

beta=500

omega=3

# R=(3*g*(n-1)/2/K)**(1/3)
# print(R)
# miu=K/2*(R**2)

def abs2(u):
    unew=np.array(u)
    unew=unew.reshape(-1)
    return(np.multiply(abs(unew),abs(unew)))
def norm(u):
    norm_val = dx*dy*(np.sum(abs2(u)))
    return norm_val

def normalize(u):
    norm_val = norm(u)
    u_normalized =u*1/math.sqrt(norm(u))
    return(u_normalized)
#stan początkowy
def psi_0(p,k):
    psi=np.zeros((Ny,Nx),dtype=complex)
    for i in range(0,Ny):
        for j in range(0,Nx):
            psi[i][j]=cmath.exp(-p*(X[j]**2+Y[i]**2)+1j*k*math.sqrt(X[j]**2+Y[i]**2))

    psi=psi.flatten()
    psi=normalize(psi)
    psi=psi.reshape(-1, 1)
    return(psi)

#kinetic operator
def kin():
    
    kin_x=-1/(2*dx**2)*sparse.diags([np.ones(Nx*Ny)*(-2),np.ones(Nx*Ny-1),np.ones(Nx*Ny-1)],[0,1,-1])
    
    kin_y=-1/(2*dy**2)*sparse.diags([np.ones(Nx*Ny)*(-2),np.ones(Nx*Ny-Ny),np.ones(Nx*Ny-Ny)],[0,Ny,-Ny])
    
    kin=kin_y+kin_x
    return( kin)
def macierz_gestości(f):
    macierz_gestości=abs2(f)
    return(macierz_gestości)

def H_Nl(u):
    unew1=np.array(u)
    unew1=unew1.reshape(-1)
    H_Nldiag=beta*macierz_gestości(u)
    H_Nlsparese=sparse.diags([H_Nldiag],[0],dtype=complex)
    return(H_Nlsparese)
V=np.zeros((Ny,Nx))
V0=15
kappa1=(math.pi)/3
kappa2=math.pi/3

#genarcaj loswego

random_numbers = np.random.normal(1, 0, Nx*Ny)

for i in range(0,Ny):
    for j in range(0,Nx):
            V[i][j]= V[i][j]+(((X[j])**2+(Y[i])**2))
            
iterator=0        
for i in range(0,Ny):
    for j in range(0,Nx):
        V[i][j]= V[i][j]+ V0*((math.sin(kappa1*X[j]))**2+(math.sin(kappa2*Y[i]))**2)* random_numbers[iterator]
        iterator+=1

V=V.flatten()
V[0:Nx]=10**15
V[-Nx:]=10**15


for i in range(0,len(V)):
    if((i+1)%Nx==0 or i%Nx==0):
        V[i]=10**15

V=sparse.diags([V],[0],dtype=complex)


#identity
def iden(a):
    return (a*(sparse.eye(Ny*Nx)))
#x-y

def L1(): 
    oy=(1j/dy/2)*(sparse.diags([np.ones(Nx*Ny-Ny),np.ones(Nx*Ny-Ny)*(-1)],[Ny,-Ny]))
    przekatna=[]#liczy w którym x jestesmy
    for i in range(0,Nx*Ny):
        przekatna.append(X[(i%Nx)])
  #  print(ox.toarray()) #MS: odkomentuj ten wiersz jak bedziesz chcial sprawdzać czy dostajesz prawidlowa macierz
    L=np.dot(sparse.diags([przekatna],[0],dtype=complex),oy)
    
    return(omega*L)
def L2():
    ox=(1j/dx/2)*(sparse.diags([np.ones(Nx*Ny-1),np.ones(Nx*Ny-1)*(-1)],[1,-1]))
    przekatna=[]#liczy w którym y jestesmy
    for i in range(0,Nx*Ny):
        przekatna.append(Y[(i//Nx)])
    L=np.dot(sparse.diags([przekatna],[0],dtype=complex),ox)
  #  print(ox.toarray())

    return(-omega*L)   




def H(u):
     return(kin()+V+H_Nl(u)+L1()+L2())
    
    
def Hprim(u):
    return(kin()+V+H_Nl(u)+L1()+L2()+iden(1j/dt))

#MS: Gdzieś w dole odwoływałeś się do Hprim ale nigdzie nie był zdefiniowany. Hprim ma dodatkowy element iden 
# w porównaniu z H. Jeśli zostawisz dt urojone to odpowiednio też musisz zmienić 1/dt na 1j/dt
def Moment_pedu(u):
    u.reshape(Nx*Ny,1)
    hermit=np.conj(u.T)
    hermit=hermit.reshape(-1)
    #print("conj u:",hermit.shape)
    Lmatrix=(L1()+L2()   )
    #print("Hmatrix:",Hmatrix.shape)
    moment_tab1= Lmatrix*u
    #print("Hmatrix*u:",enrgia_tab1.shape)
   
    moment_tab2=(np.dot(hermit,moment_tab1))
    #print("u* H u:",enrgia_tab2.shape)
    moment_tab2=moment_tab2.reshape(-1)
    #print("result", enrgia_tab2.shape)
    return dx*dx*(moment_tab2.real)


# def Energia(u):
#     u.reshape(Nx*Ny,1)
#     hermit=np.conj(u.T)
#     #hermit=hermit.reshape(-1)
#     print("conj u:",hermit.shape)
#     Hmatrix=(H(u))
#     print("Hmatrix:",Hmatrix.shape)

#     enrgia_tab1= Hmatrix*u
#     print("Hmatrix*u:",enrgia_tab1.shape)
   
#     enrgia_tab2=(np.dot(hermit,enrgia_tab1))
#     print("u* H u:",enrgia_tab2.shape)
#    # enrgia_tab2=enrgia_tab2.reshape(-1)
#     print("result", enrgia_tab2.shape)
#     return dx*dx*(enrgia_tab2.real)

def delta_E(u1,u2):
    diff= (np.sum(np.abs(abs2(u1)-abs2(u2))))
    return diff
def solve_num(psi):
    psiw=[]
    dE=1000
    iterator=0
    print("epsilon dE = ", e)
    while dE>e: # MS: dobrze dodać tutaj warunek na maksymalną ilość iteracji (np 100 000) 

        
        psi_prev=psi
        psi = spsolve(Hprim(psi_prev),psi_prev*(1j/dt)) 
        #MS: tutaj zabraklo podzielenia przez dt lub przez dt/j jeśli dt jest urojone. Zmienilem też żeby funkcja
        #przyjmowała tylko jeden argument dla przejrzystości.
        
        psi=normalize(psi)
        
        dE=delta_E(psi,psi_prev)
        iterator+=1
        if(iterator%200==0):
            print(iterator)
            print(dE)
        if(iterator<200 and iterator%10==0):
            psiw.append(psi)
	    psi_L.append(Moment_pedu(psi))
        if(iterator>100000):
            return(psiw)
            break
    return(psiw)
psi_t=psi_0(0.02,0.02)
psi_t = np.reshape(abs2(psi_t),(Nx,Ny))
im = plt.imshow( psi_t )
psi_t = np.reshape(abs2(psi_t),Nx*Ny)
psi_L=[]
psiwynik=solve_num(psi_t)

np.save("C=0,R=401,12,O=3,B==500,V0=15,kpina3,ewolcuaja_czasowa_psi",, psiwynik)
np.save("C=0,R=401,12,O=3,B==500,V0=15,kpina3,ewolcuaja_czasowa_L",psi_L)             