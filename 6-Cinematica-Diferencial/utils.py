import sympy as sp
import numpy as np
import matplotlib.pyplot as plt

def srotx(ang):
    """ Matriz de rotación alrededor de x
    """
    Rx = sp.Matrix([[1, 0, 0],
                    [0, sp.cos(ang), -sp.sin(ang)],
                    [0, sp.sin(ang), sp.cos(ang)]])
    return Rx

def sroty(ang):
    """ Matriz de rotación alrededor de y
    """
    Ry = sp.Matrix([[sp.cos(ang), 0, sp.sin(ang)],
                    [0, 1, 0],
                    [-sp.sin(ang), 0, sp.cos(ang)]])
    return Ry

def srotz(ang):
    """ Matriz de rotación alrededor de z
    """
    Rz = sp.Matrix([[sp.cos(ang),-sp.sin(ang),0],
                    [sp.sin(ang), sp.cos(ang),0],
                    [0,0,1] ])
    return Rz

def sTdh(d, th, a, alpha):
    """ Matriz de transformación homogénea de Denavit-Hartenberg
    """
    cth = sp.cos(th); sth = sp.sin(th)
    ca = sp.cos(alpha); sa = sp.sin(alpha)
    Tdh = sp.Matrix([[cth, -ca*sth,  sa*sth, a*cth],
                     [sth,  ca*cth, -sa*cth, a*sth],
                     [0,        sa,     ca,      d],
                     [0,         0,      0,      1]])
    return Tdh

def sVectorFromSkew(S):
    return sp.Matrix([S[2,1],S[0,2],S[1,0]])



def Trotz(ang):
    """ Transformación homogénea que representa rotación en z
    """
    T = np.array([[np.cos(ang), -np.sin(ang), 0., 0.],
                  [np.sin(ang),  np.cos(ang), 0., 0.],
                  [0., 0., 1., 0.],
                  [0., 0., 0., 1.]])
    return T

def Trasl(x, y, z):
    """ Transformación homogénea que representa traslación pura
    """
    T = np.array([[1,0,0,x],
                  [0,1,0,y],
                  [0,0,1,z],
                  [0,0,0,1]])
    return T

def cdirecta_robot2d(q, l1, l2):
    """ Retorna los sistemas de referencia de cada eslabón con respecto a la base
    """
    # Sistemas con respecto al anterior
    T01 = Trotz(q[0]).dot(Trasl(l1,0,0))
    T1e = Trotz(q[1]).dot(Trasl(l2,0,0))
    # Sistemas con respecto a la base
    T0e = T01.dot(T1e)
    return T01, T0e

def plot_robot2d(q, l1, l2, ax=None, show_axis=True, k=0.25):
    """ Grafica el robot RR (2D) según la configuración articular
    """
    T1, Te = cdirecta_robot2d(q, l1, l2)
    if (ax==None):
        plt.clf()
        ax = plt
    # Cuerpo del robot
    ax.plot([0, T1[0,3]],[0, T1[1,3]], linewidth=3, color='b')
    ax.plot([T1[0,3], Te[0,3]],[T1[1,3], Te[1,3]], linewidth=3, color='b')
    # Puntos en las articulaciones
    ax.plot(0, 0, color='r', marker='o', markersize=8)
    ax.plot(T1[0,3], T1[1,3], color='r', marker='o', markersize=8)
    # Efector final (definido por 4 puntos)
    p1 = np.array([0, 0.1, 0, 1]); p2 = np.array([0.2, 0.1, 0, 1])
    p3 = np.array([0, -0.1, 0, 1]); p4 = np.array([0.2, -0.1, 0, 1])
    p1 = Te.dot(p1); p2 = Te.dot(p2); p3 = Te.dot(p3); p4 = Te.dot(p4)
    # Gráfico del efector final
    ax.plot([p1[0],p2[0]], [p1[1],p2[1]], color='b', linewidth=3)
    ax.plot([p3[0],p4[0]], [p3[1],p4[1]], color='b', linewidth=3)
    ax.plot([p1[0],p3[0]], [p1[1],p3[1]], color='b', linewidth=3)
    # Sistema de referencia del efector final y de la base
    if (show_axis):
        ax.plot([Te[0,3],Te[0,3]+k*Te[0,0]], [Te[1,3],Te[1,3]+k*Te[1,0]],color='r')
        ax.plot([Te[0,3],Te[0,3]+k*Te[0,1]], [Te[1,3],Te[1,3]+k*Te[1,1]],color='g')
        ax.plot([0,k],[0,0],color='r');plt.plot([0,0],[0,k],color='g')
    # Plot
    ax.axis('equal')
    ax.grid('on')