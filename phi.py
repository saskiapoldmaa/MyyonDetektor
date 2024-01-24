import matplotlib.pyplot as plt
import numpy as np
import math

E = np.logspace(0, 3, 500) #in GeV/c^2
E_mev=[1000*q for q in E]
p=[x for x in E]

#Bugaev-Reyna model from table 1
bugaev=[]
ad=[[0.3061, 1.2743, -0.2630, 0.0252], [1.7910, 0.3040, 0, 0], [3.6720, 0, 0, 0], [4, 0, 0, 0]]
Ad=[2.950e-3,1.781e-2,14.35,1000]
for i in p:
    if 0<=i<9.2765e2:
        a=ad[0]
        A=Ad[0]
        bugaev.append(A*i**-(a[3]*np.log10(i)**3+a[2]*np.log10(i)**2+a[1]*np.log10(i)+a[0]))

    elif 9.2765e2<=i<1.5878e3:
        a=ad[1]
        A=Ad[1]
        bugaev.append(A*i**-(a[3]*np.log10(i)**3+a[2]*np.log10(i)**2+a[1]*np.log10(i)+a[0]))

    elif 1.5878e3<=i<4.1625e5:
        a=ad[2]
        A=Ad[2]
        bugaev.append(A*i**-(a[3]*np.log10(i)**3+a[2]*np.log10(i)**2+a[1]*np.log10(i)+a[0]))

    elif i>=4.1625e5:
        a=ad[3]
        A=Ad[3]
        bugaev.append(A*i**-(a[3]*np.log10(i)**3+a[2]*np.log10(i)**2+a[1]*np.log10(i)+a[0]))

#Smith and Duller-Chatzidaki model from table 2
phi=[]
A,r,a,y0,gamma=2.382e-3,0.76,2.5,1000,-8/3
bmu,mmu,taumu,rho0,c=0.8,105.659,2.2e-6,1.23e-3,3e8
lampi,b,tau0,mpi=120,0.771,2.61e-8,139.58
for y in E_mev:
    Bmu=bmu*mmu*3e10*1/(taumu*3.321e28)
    Epi=1/r*(y+0.9*a*y0)
    jpi=mpi*y0*3e10/(tau0*rho0)
    Pmu=(0.1*(1-a*(y0-100)/(r*Epi)))**(Bmu/(r*Epi+100*a))
    phi.append(A*Epi**gamma*Pmu*lampi*b*jpi/(Epi*y+b*jpi))
#new chatzidaki
chatzidaki=[]
for q in E:
    #Epi=1/r*(q+0.9*a*y0) #võetud kokkuvõrvast artiklist ja on probleemi juur
    Epi=q
    y0,R,Te,M,g,rho0,c=1000,8.314,220,28.966,981.3,0.00123,3e10
    alpha,mmu,mpi,mk,taumu,tau0,tauk,b,bmu,lammu,lampi,lamk,jpi,jk=2.5,105.659,139.58,493.8,2.2e-6,2.61e-8,1.24e-8,0.771,0.8,120,120,120,148.16,1105.80
    A,gamma=0.002382,2.645
    jpi=mpi*y0*c/tau0/rho0
    Bmu=bmu*mmu*y0*c/(taumu*rho0)
    #Pmu=(y0*q/(y0*q))**(Bmu/(q*1))
    Pmu=0.5
    chatzidaki.append(A*Epi**(-gamma)*Pmu*lampi*b*jpi/(Epi*1+b*jpi))
#Tang model
p1,p2,p3,p4,p5=0.102573,-0.068287,0.958633,0.0407253,0.817285
cos=math.sqrt((1+0.102573**2-0.068287+0.0407253)/(1+0.102573**2-0.068287+0.0407253))
epspi,epskap,gam,BG=115/1.1,819/1.1,-2.7,0.054
tang=[]
for z in E:
    if z>100/cos:
        AT,rc=1,0
        Eroof=z
    elif 1/cos<z<=100/cos:
        Eroof=z+2.06e-3*(950/cos-90)
        AT,rc=1.1*(90*math.sqrt(1.001)/1030)**(4.5/(z*cos)),1e-4
    elif z<=1/cos:
        AT=0.001
        Eroof=0.3*z+0.7
        rc=0
    tang.append(AT*0.14*z**gam*(1/(1+Eroof*cos/epspi)+BG/(1+Eroof*cos/epskap)+rc))
#Gaisser model
gaisser=[]
for x in E:
    gaisser.append(0.14*x**(-2.7)*(1/(1+1.1*x/115)+0.054*(1+1.1*x/810)))
plt.plot(E,bugaev) #sirgem
plt.plot(E, chatzidaki)
plt.plot(E,tang) #jõnksuga
plt.plot(E,gaisser) #samuti sirge
plt.legend(["Bugaev", "Chatzidaki", "Tang", "Gaisser"])
plt.yscale('log')
plt.xscale('log')
plt.title('Müüonvoog maapinnal')
plt.ylabel('Voog [$GeV^{-1}cm^{-2}s^{-1}sr^{-1}]$')
plt.xlabel('Energia[GeV]')
plt.grid(True, which="both", ls="-", color='0.65')
plt.show()


