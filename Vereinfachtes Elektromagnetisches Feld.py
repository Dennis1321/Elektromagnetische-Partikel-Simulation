import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.widgets import Slider


#Positionen
x=1e-10
x2=0.0
y=0.0
y2=0.0
#Masse
m=1.672621923e-27
m2=9.1093837139e-31
#Geschwindigkeit
vx=0.0
vx2=0
vy2=2e6
vy=-vy2 * (m2 / m)

dt=2e-19#Zeitschritt
#Elektrische Ladung
Q=-1.602176634e-19
q=1.602176634e-19
Epsilon_0=8.854187817e-12 #Feldkonstante

Bz = 30 # Magnetfeld

#Elektrisches Feld
def field(x,y): 
    #Abstand
    dx=x-x2
    dy=y-y2
    r=(dx**2+dy**2+1e-40)**0.5#Radius
    
    vorfaktor= 1 / (4 * np.pi * Epsilon_0) * (q * Q) / (r ** 3)#Kraftbetrag
    #Kraft
    f_x=dx*vorfaktor
    f_y=dy*vorfaktor

    return f_x,f_y

def energy():
    dx=x-x2
    dy=y-y2
    r=(dx**2+dy**2+1e-40)**0.5
    #Kinetische Energie
    kin_energy_pro=1/2*m*(vx**2+vy**2)
    kin_energy_ele=1/2*m2*(vx2**2+vy2**2)

    pot_energy= 1 / (4 * np.pi * Epsilon_0) * (q * Q) /r#Potenzielle Energie

    return pot_energy+kin_energy_ele+kin_energy_pro

#Grafik
fig, ax = plt.subplots()

plt.subplots_adjust(bottom=0.18)

electron, = ax.plot([0], [0], 'go', markersize=3,label="Elektron")
trail_ele, = ax.plot([],[], 'g-',alpha=0.5)
proton, = ax.plot([x], [y], 'ro', markersize=6,label="Proton")
trail_pro, = ax.plot([], [], 'r-',alpha=0.5)

energy_text = ax.text(-4.5e-10, 4.2e-10, '', fontsize=10, fontfamily='monospace')

#Koordinatensystem
ax.set_xlim(-6e-10,6e-10)
ax.set_ylim(-6e-10,6e-10)
ax.set_aspect('equal')
ax.legend()
#Trail
xs, ys = [x], [y]
xs2,ys2 = [x2], [y2]

ax_slider = plt.axes([0.2, 0.05, 0.6, 0.03])
slider_B = Slider(
    ax=ax_slider,
    label='Magnetfeld in T',
    valmin=-1000,
    valmax=1000,
    valinit=Bz
)

def update_B(val):
    global Bz
    Bz = slider_B.val

slider_B.on_changed(update_B)

def update(frame):
    global x, y, vx, vy,x2,y2,vx2,vy2
    
    for _ in range(100):
        feldx, feldy = field(x,y)#Elektrisches Feld
        
        # Halb-Schritt (E-Feld)
        vx_minus = vx + (feldx / m) * (dt / 2)
        vy_minus = vy + (feldy / m) * (dt / 2)
        vx2_minus = vx2 + (-feldx / m2) * (dt / 2)
        vy2_minus = vy2 + (-feldy / m2) * (dt / 2)

        # Boris Rotation (Magnetfeld)
        t_pro = (q * Bz / m) * (dt / 2)
        t_ele = (Q * Bz / m2) * (dt / 2)
        
        s_pro = (2 * t_pro) / (1 + t_pro**2)
        s_ele = (2 * t_ele) / (1 + t_ele**2)
        
        v_prime_x_pro = vx_minus + vy_minus * t_pro
        v_prime_y_pro = vy_minus - vx_minus * t_pro
        v_prime_x_ele = vx2_minus + vy2_minus * t_ele
        v_prime_y_ele = vy2_minus - vx2_minus * t_ele
        
        vx_plus = vx_minus + v_prime_y_pro * s_pro
        vy_plus = vy_minus - v_prime_x_pro * s_pro
        vx2_plus = vx2_minus + v_prime_y_ele * s_ele
        vy2_plus = vy2_minus - v_prime_x_ele * s_ele

        # Zweiter Halb-Schritt (E-Feld)
        vx = vx_plus + (feldx / m) * (dt / 2)
        vy = vy_plus + (feldy / m) * (dt / 2)
        vx2 = vx2_plus + (-feldx / m2) * (dt / 2)
        vy2 = vy2_plus + (-feldy / m2) * (dt / 2)

        # Position
        x += vx * dt
        y += vy * dt
        x2 += vx2 * dt
        y2 += vy2 * dt
    #Trail übernehmen
    xs.append(x)
    ys.append(y)
    xs2.append(x2)
    ys2.append(y2)
    #Werte übernehmen
    proton.set_data([x], [y])
    trail_pro.set_data(xs, ys)
    electron.set_data([x2], [y2])
    trail_ele.set_data(xs2, ys2)

    energy_text.set_text(f"Gesamtenergie: {energy():.6e} J")#Text
    
    return proton, trail_pro, electron, trail_ele, energy_text

ani=FuncAnimation(fig,update,frames=1000000,interval=1,blit=False)#Animationseinstellungen

plt.show()#Animation anzeigen