import numpy as np

cos = np.cos
sin = np.sin
pi = np.pi


def Trotx(ang):
    """ Transformación homogénea que representa rotación en x
    """
    T = np.array([[1., 0., 0.,0.],
                  [0., cos(ang), -sin(ang), 0.],
                  [0., sin(ang), cos(ang), 0.],
                  [0., 0., 0., 1.]])
    return T

def Troty(ang):
    """" Transformación homogénea que representa rotación en y
    """
    T = np.array([[cos(ang), 0., sin(ang), 0.],
                  [0., 1., 0., 0.],
                  [-sin(ang), 0., cos(ang), 0.],
                  [0., 0., 0., 1.]])	
    return T

def Trotz(ang):
    """ Transformación homogénea que representa rotación en z
    """   
    T = np.array([[cos(ang), -sin(ang), 0., 0.],
                  [sin(ang), cos(ang), 0., 0.],
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

def skew(w):
    """ Matriz antisimétrica (skew) a partir de un vector
    """
    S = np.array([[0., -w[2], w[1]],
                  [w[2], 0., -w[0]],
                  [-w[1], w[0], 0.]])
    return S    
