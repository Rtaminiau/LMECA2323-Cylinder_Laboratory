import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
import scipy.linalg as solver


def find_x(pol, y):
    """

    :param pol: polynome from np.poly1D du second degré
    :param y:  ordonnée du polynome tel que pol(x)= y
    :return: x1 & x2 les racines du polynôme
    """
    a = pol.c[0]
    b = pol.c[1]
    c = pol.c[2] - y

    delta = b * b - 4 * a * c

    x1 = (-b + np.sqrt(delta)) / (2 * a)
    x2 = (-b - np.sqrt(delta)) / (2 * a)

    return x1, x2
txt = np.loadtxt("labo1.txt")

angle = txt[:,0]
pression = txt[:,1]
Md = np.array([0.0, 0.0682, 0.182, 0.295, 0.408,0.522, 0.748, 0.975])
Ud = np.array([0.38, 0.7, 1.23, 1.76, 2.3,2.82, 3.89, 4.9]) 

model = np.poly1d(np.polyfit(angle[0:5], pression[0:5], 2))

axe = np.linspace(-10,10,40)
angle_max = axe[model(axe).argmax()]
pressure_max = max(model(axe))
pressure_min = min(model(axe))

#Question 2 stagnation point
def plot_Q2():
    fig, ax = plt.subplots()
    ax.scatter(angle[0:5], pression[0:5])
    ax.plot(axe, model(axe))
    ax.vlines(angle_max,pressure_min-5,pressure_max,linestyles="dashed",colors='red')
    ax.hlines(pressure_max,-10,angle_max,linestyles="dashed",colors='red')

    ax.set_xlabel("Angle [deg]")
    ax.set_ylabel(" p(θ) - p$\infty$ [Pa]")
    #ax.set_title("Evolution de la pression dynamique  \n autour du point de stagnation du cylindre")
    ax.set_xlim(-10,11)
    ax.set_ylim(pressure_min-5,pressure_max+5)
    ax.tick_params(axis='both', which='minor', colors='red')
    ax.yaxis.set_major_locator(ticker.FixedLocator([115,120,125,130,135,140,145]))
    ax.yaxis.set_minor_locator(ticker.FixedLocator([pressure_max]))

    plt.savefig("Stagnation.pdf")
    plt.show()

#Question 1 U_inf et Re
rho = 1.125
nu = 15.6e-6
D = 0.05
b = 0.5

U = np.sqrt(2*pressure_max/rho)
Re = D*U/nu

#Question 4 Cp avec la balance

Tension_mesure = 2.62
C_De = 0.09
Md = Md*9.81
model_force = np.poly1d(np.polyfit(Ud, Md, 1))

tension = np.linspace(Ud[0],5,100)
Force = model_force(tension)
Force_mesure = model_force(Tension_mesure)
def plot_force_balance():
    plt.plot(tension,Force)
    plt.scatter(Ud,Md)
    plt.plot(Tension_mesure,Force_mesure,'.r')
    plt.vlines(Tension_mesure,min(Md),Force_mesure,linestyles="dashed",colors="red")
    plt.hlines(Force_mesure,0,Tension_mesure,linestyles="dashed",colors="red")
    plt.xlim(Ud[0],5)
    plt.ylim(0,10)

    plt.xlabel("Tension [V]")
    plt.ylabel("Force [N]")
    #plt.title("Calibration de la balance, Force/tension")

    plt.savefig("Balance_drag.pdf")
    plt.show()

Total_drag = model_force(Tension_mesure)
Total_CD = Total_drag/(1/2 * rho * U**2 * b*D)
Drag_cylinder = Total_CD - C_De

#Question 3 CP coeff et drag


x = np.linspace(-10,180,100)
y = 1 - 4 * np.sin(x * np.pi/180) **2


plt.scatter(angle,pression/(1/2 * rho * U**2))
plt.plot(angle,pression/(1/2 * rho * U**2), label="expérimental")
plt.plot(x,y,'r', label="théorie")
plt.xlabel("Angle [deg]")
plt.ylabel("Cp coefficient [-]")
plt.legend()
plt.savefig("Comp.pdf")
plt.show()

CD = 0
for i in range(2,len(pression)-1):
    CDi1 = pression[i]/((1/2) * rho * U**2) * np.cos(np.pi/180 * angle[i])
    CDi2 = pression[i+1] / ((1 / 2) * rho * U ** 2) * np.cos(np.pi/180 * angle[i+1])
    CD += (CDi2+CDi1)/2 * (angle[i+1] - angle[i]) * np.pi/180

print("Here is CD : ",CD)

#Prints==============================

print("max angle pression",angle_max)
print("Max pressure",pressure_max)

print("Upstream velocity",U)
print("Reynolds number ",Re)


print("Force balance aero",Force_mesure)
print("Total drag coeff",Total_CD)
print("Drag cylinder coeff ",Drag_cylinder)
plot_Q2()
plot_force_balance()



