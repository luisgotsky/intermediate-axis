"""
Código: Luis Lucas García
Grado en Física - Prácticas de mecánica newtoniana
Práctica 3 - Teorema de la raqueta de ténis
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import matplotlib.animation as anim

#Constantes de Python

ql = 1200
size=(8, 6)

#Constantes físicas

I1 = 4
I2 = 2
I3 = 1

#Resolución del problema con odeint

def omega(z, t):
    
    o1, o2, o3 = z
    
    return [(I3 - I2)*o2*o3/I1, (I1 - I3)*o1*o3/I2, (I2 - I1)*o1*o2/I3]

o1_0 = 2
o2_0 = 0.1
o3_0 = 0.1

z_0 = [o1_0, o2_0, o3_0]
nt = np.linspace(0, 30, 25000)

z = odeint(omega, z_0, nt)

o1 = z[:,0]
o2 = z[:,1]
o3 = z[:,2]

#Buscamos la frecuencia de las oscilaciones y comparamos con la predicha

expected = np.sqrt((I1 - I2)*(I1 - I3)/I2*I3)*o1_0
print("El valor esperado es: ", expected)

def velAng(o):

    cont = -1
    PrimeraVez = True
    for i in range(1, len(o)-1):
        
        if o[i-1]<o[i] and o[i]>o[i+1]:
            
            cont += 1
            
            if PrimeraVez:
                
                t0 = nt[i]
                PrimeraVez = False
                
            tf = nt[i]
            
    t = tf-t0
    return 2*np.pi*cont/t

real2, real3 = velAng(o2), velAng(o3)
print("El valor real es: ", real2, real3)

#Reducción del error

data = open("data", "r")
Prim = True
iteraciones, val1, val2 = [], [], []

for linea in data:
    
    if not Prim:
        
        cad = linea.split()
        iteraciones.append(float(cad[0]))
        val1.append(float(cad[1]))
        val2.append(float(cad[2]))
        
    else: Prim = False

data.close()

pol1 = np.polyfit(iteraciones, val1, 3)
pol2 = np.polyfit(iteraciones, val2, 2)

#Graficamos

plt.figure()
plt.scatter(iteraciones, val1, label="$val_1$")
plt.scatter(iteraciones, val2, label="$val_2$")
plt.plot(iteraciones, np.polyval(pol1, iteraciones), label="aprox. 1")
plt.plot(iteraciones, np.polyval(pol2, iteraciones), label="aprox. 2")
plt.plot(iteraciones, [expected for i in range(len(iteraciones))], "--")
plt.legend(loc="upper right")
plt.xlabel("Iteraciones")
plt.ylabel("Frecuencia")
plt.savefig("tendencia.png", dpi=ql)

#Estudiamos las energías y el momento angular

E = 0.5*np.sqrt(I1*o1**2 + I2*o2**2 + I3*o3**2)
Lx, Ly, Lz = I1*o1, I2*o2, I3*o3
L = np.sqrt(Lx**2 + Ly**2 + Lz**2)

#Deducimos el periodo de la oscilación

PerOx = 2*np.pi/velAng(o1)
PerOy = 2*np.pi/velAng(o2)
PerOz = 2*np.pi/velAng(o3)
PerLx = 2*np.pi/velAng(Lx)
PerLy = 2*np.pi/velAng(Ly)
PerLz = 2*np.pi/velAng(Lz)

print("El periodo de la velocidad angular es: ", PerOx, PerOy, PerOz)
print("El periodo del momento angular es: ", PerLx, PerLy, PerLz)

#Graficamos las velocidades angulares

plt.figure(figsize=size)
plt.plot(nt, o1, label="$\\omega_1$")
plt.plot(nt, o2, label="$\\omega_2$")
plt.plot(nt, o3, label="$\\omega_3$")
plt.legend(loc="lower right")
#plt.savefig("omegaInt.png", dpi=ql)

#Animamos la velocidad angular

fig2 = plt.figure()
ax2 = plt.axes(projection="3d")
ax2.set_xlim([min(o1)-0.5, max(o1)+0.5])
ax2.set_ylim([min(o2)-0.5, max(o2)+0.5])
ax2.set_zlim([min(o3)-0.5, max(o3)+0.5])
q2 = ax2.quiver(0, 0, 0, o1[0], o2[0], o3[0])
line2, = ax2.plot(o1[:1], o2[:1], o3[:1])

def vectAnim(i):
    
    ax2.cla()
    ax2.set_xlim([min(o1)-0.5, max(o1)+0.5])
    ax2.set_ylim([min(o2)-0.5, max(o2)+0.5])
    ax2.set_zlim([min(o3)-0.5, max(o3)+0.5])
    q2 = ax2.quiver(0, 0, 0, o1[25*i], o2[25*i], o3[25*i], length=0.5)
    line2, = ax2.plot(o1[:25*i+1], o2[:25*i+1], o3[:25*i+1])
    ims = [[q2, line2]]
    return ims

ani2 = anim.FuncAnimation(fig2, vectAnim, frames=1000, interval=10)
#ani2.save("oInt.gif")

#Graficamos las energías y el momento angular

plt.figure(figsize=size)
plt.plot(nt, E, label="E (J)")
plt.plot(nt, L, "--", label="L ($kg\\frac{m^2}{s}$)")
plt.legend(loc="lower right")
#plt.savefig("EyLInt", dpi=ql)

fig = plt.figure(figsize=size)
ax = plt.axes(projection="3d")
q = plt.quiver(0, 0, 0, Lx[0], Ly[0], Lz[0])
ax.set_xlim([0, 15])
ax.set_ylim([0, 3])
ax.set_zlim([0, 3])
fig.add_axes(ax)
#plt.savefig("Vector.png", dpi=ql)

#Animamos el vector momento angular

fig = plt.figure()
ax = plt.axes(projection="3d")
ax.set_xlim([min(Lx)-0.5, max(Lx)+0.5])
ax.set_ylim([min(Ly)-0.5, max(Ly)+0.5])
ax.set_zlim([min(Lz)-0.5, max(Lz)+0.5])
q = ax.quiver(0, 0, 0, Lx[0], Ly[0], Lz[0])
line, = ax.plot(Lx[:1], Ly[:1], Lz[:1])
aL = 25

def vectAnim(i):
    
    ax.cla()
    ax.set_xlim([min(Lx)-0.5, max(Lx)+0.5])
    ax.set_ylim([min(Ly)-0.5, max(Ly)+0.5])
    ax.set_zlim([min(Lz)-0.5, max(Lz)+0.5])
    q = ax.quiver(0, 0, 0, Lx[aL*i], Ly[aL*i], Lz[aL*i], length=0.4)
    line, = ax.plot(Lx[:aL*i+1], Ly[:aL*i+1], Lz[:aL*i+1])
    ims = [[q, line]]
    return ims

ani = anim.FuncAnimation(fig, vectAnim, frames=len(Lx)//aL, interval=10)
#ani.save("LInter.gif")