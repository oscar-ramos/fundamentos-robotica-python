"""
Forward Kinematics using Denavit-Hartenberg convention

Oscar Ramos Ponce
Universidad de Ingenieria y Tecnologia - UTEC

"""

import numpy as np
import sympy as sp

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from matplotlib.patches import Circle
import mpl_toolkits.mplot3d.art3d as art3d


def plot_cylinder(ax, p0, p1, radius, resolution=6, alpha=1, color='g'):
    """
	Plot a cylinder in 3D
	
	    ax - An axis for the plot
		p0, p1 - 3D points that specify the initial and final point of the axis
		radius - radius of the cylinder
        resolution - number of polygons used to approximate the circular surface
		alpha - level of transparency
		color - color of the cylinder
		
	"""
    # Vector in the direction of the axis
    v = np.array([p1[0]-p0[0], p1[1]-p0[1], p1[2]-p0[2]])
    # Magnitude of the vector (length of the cylinder)
    length = np.linalg.norm(v)
    # Unit vector in the direction of the axis
    v = v/length
    # A vector not in the same direction as v
    not_v = np.array([1, 0, 0])
    if (v == not_v).all():
        not_v = np.array([0, 1, 0])
    # Make the vector perpendicular to v
    n1 = np.cross(v, not_v)
    # Normalize n1
    n1 /= np.linalg.norm(n1)
    # Make a unit vector perpendicular to v and n1
    n2 = np.cross(v, n1)

    # The surface ranges over t from 0 to the length of the axis and 0 to 2*pi
    t = np.linspace(0, length, 2)
    th = np.linspace(0, 2*np.pi, resolution+1)
    rsample = np.linspace(0, radius, 2)
    # Use meshgrid to make 2d arrays
    t, th2 = np.meshgrid(t, th)
    rsample, th = np.meshgrid(rsample, th)

    # Coordinates for the body of the cylinder
    X1, Y1, Z1 = [p0[i]+v[i]*t + radius*np.sin(th2)*n1[i] + radius*np.cos(th2)*n2[i] for i in [0, 1, 2]]
    # Coordinates for the bottom of the cylinder
    X2, Y2, Z2 = [p0[i] + rsample[i]*np.sin(th)*n1[i] + rsample[i]*np.cos(th)*n2[i] for i in [0, 1, 2]]
    # Coordinates for the top of the cylinder
    X3, Y3, Z3 = [p0[i] + v[i]*length + rsample[i]*np.sin(th)*n1[i] + rsample[i]*np.cos(th)*n2[i] for i in [0, 1, 2]]
    # Plot
    ax.plot_surface(X1, Y1, Z1, alpha=alpha, color=color)
    ax.plot_surface(X2, Y2, Z2, alpha=alpha, color=color)
    ax.plot_surface(X3, Y3, Z3, alpha=alpha, color=color)
	

class SerialRobot(object):
    def __init__(self, L, name='robot'):
        """
        Constructor for a generic serial robot
		
          L - List containing the DH parameters in the following format
              [[d1, th1, a1, alpha1, rp1], [d2, th2, a2, alpha2, rp2],...]
              where rp can be 'r' for revolute joint, or 'p' for prismatic joint.
			  
          name - String containing the robot name
          
        """
        self.name = name
        self.ndof = len(L)
        self.type = [i[4] for i in L]
        # DH parameters
        self.d = [i[0] for i in L] 
        self.th = [i[1] for i in L]
        self.a = [i[2] for i in L]
        self.alpha = [i[3] for i in L]
        # Store the transformation matrices
        self.Ts = self.ndof*[0,]
        # Store indication that d is zero
        self.isdzero = self.ndof*[True, ]
        for k in range(self.ndof):
            if self.type[k] == 'p':
                self.isdzero[k] = False
            if np.abs(self.d[k])>1e-6:
                self.isdzero[k] = False
        # Store points for plot
        self.p = (self.isdzero.count(False)+self.ndof)*[0,]
        
    def _Tdh(self, k, d, th):
        cth = np.cos(self.th[k]+th); sth = np.sin(self.th[k]+th)
        ca = np.cos(self.alpha[k]); sa = np.sin(self.alpha[k])
        return np.array([[cth, -ca*sth, sa*sth, self.a[k]*cth],
                         [sth, ca*cth, -sa*cth, self.a[k]*sth],
                         [0., sa, ca, self.d[k]+d],
                         [0., 0., 0., 1.]])
    
    def _sym_Tdh(self, k):
        cth = sp.cos(self.th[k]); sth = sp.sin(self.th[k])
        ca = sp.cos(self.alpha[k]); sa = sp.sin(self.alpha[k])
        return Matrix([[cth, -ca*sth, sa*sth, self.a[k]*cth],
                       [sth, ca*cth, -sa*cth, self.a[k]*sth],
                       [0, sa, ca, self.d[k]],
                       [0, 0, 0, 1]])

    def get_name(self):
        """Returns the robot name
        """
        return self.name
        
    def sym_fkine(self, verbose=False):
        """
        Compute the symbolic forward kinematics of the robot.
           verbose - shows intermediate matrices
        
        """
        Tf = sp.eye(4)
        for k in range(self.ndof):
            if self.type[k]=='r':
                T = self._sym_Tdh(k)
                if(verbose):
                    print('\nT'+str(k)+str(k+1)+':')
                    try:
                        display(T)
                    except NameError:
                        sp.pprint(T)
            elif self.type[k]=='p':
                T = self._sym_Tdh(k)
                if(verbose):
                    print('\nT'+str(k)+str(k+1)+':')
                    try:
                        display(T)
                    except NameError:
                        sp.pprint(T)
            else:
                print('not supported joint type')
            Tf = sp.simplify(Tf*T)
        return Tf
    
    def fkine(self, q, verbose=False):
        """
        Compute the forward kinematics of the robot given a joint configuration
           q - Joint configuration
           verbose - shows intermediate matrices
        
        """
        if len(q)!=self.ndof:
            print("Error: incorrect number of joints")
            return 0
        Tf = np.eye(4)
        for k in range(self.ndof):
            if self.type[k]=='r':
                T = self._Tdh(k, 0., q[k])
                if(verbose):
                    print('\nT'+str(k)+str(k+1)+':')
                    print(np.round(T,3));
            elif self.type[k]=='p':
                T = self._Tdh(k, q[k], 0.)
                if(verbose):
                    print('\nT'+str(k)+str(k+1)+':')
                    print(np.round(T,3));
            else:
                print('not supported joint type')
            self.Ts[k] = T
            Tf = Tf.dot(T)
        return Tf

    def plot(self, ax, q, ee=True, axlimits=None, elev=25, azim=45, ascale=0.4, cscale=0.1, radius=0.05, color='g', lw=3, colab=False):
        """
        Basic plot of the robot using lines and showing the joints with green "cylinders"
        
        Arguments
           q - Joint configuration vector (1D)
           ee - If true (default), an end-effector is plotted
           axlimits - Limits for the figure axes [[xmin,xmax],[ymin,ymax],[zmin,zmax]]. If not
                      specified, the default axes are used
           elev - elevation [in degrees] for the plot view
           azim - azimuth [in degrees] for the plot view
           ascale - scale for the axis at the base and end effector
           cscale - scale for the "cylinders"
		   radius - radius of the cylinders
		   color - color for the cylinders
		   lw - width of the lines joining the joints
           colab - True when using Google Colab
           
        """
        # Compute intermediate DH homogeneous transformations
        pcnt = 0
        # Initial intermediate point
        if (not self.isdzero[0]):
            if self.type[0]=='r':
                self.p[pcnt] = np.array([0.,0.,self.d[0],1.])
                pcnt+=1
            elif self.type[0]=='p':
                self.p[pcnt] = np.array([0.,0.,self.d[0]+q[0],1.])
                pcnt+=1
        # Loop around all the ndofs
        for k in range(self.ndof):
            if self.type[k]=='r':
                T = self._Tdh(k, 0., q[k])
                if (not self.isdzero[k] and k!=0):
                    self.p[pcnt] = np.array([0.,0.,self.d[k],1.])
            elif self.type[k]=='p':
                T = self._Tdh(k, q[k], 0.)
                if (not self.isdzero[k] and k!=0):
                    self.p[pcnt] = np.array([0.,0.,self.d[k]+q[k],1.])
            else:
                print('wrong joint type')
            if k==0:
                self.Ts[k] = T
                self.p[pcnt] = self.Ts[k][:,3]
                pcnt+=1
            else:
                self.Ts[k] = self.Ts[k-1].dot(T)
                if (not self.isdzero[k]):
                    self.p[pcnt] = self.Ts[k-1].dot(self.p[pcnt])
                    pcnt+=1
                self.p[pcnt] = self.Ts[k][:,3]
                pcnt+=1
        # Clear the figure
        # plt.clf(); ax = plt.axes(projection='3d')
        if (not colab):
            ax.cla()
        # Names for the axes
        ax.set_xlabel('x'); ax.set_ylabel('y'); ax.set_zlabel('z')
        # Points (base of the robot)
        ax.scatter(0, 0, 0, color='g', s=50)
        # Body of the robot
        ax.plot([0, self.p[0][0]], [0, self.p[0][1]], [0, self.p[0][2]], linewidth=lw, color='k')
        for k in range(1, len(self.p)):
            ax.plot([self.p[k-1][0],self.p[k][0]], [self.p[k-1][1],self.p[k][1]],
                    [self.p[k-1][2], self.p[k][2]], linewidth=lw, color='k')
        e = 0
        # Ensure cylinders for end-effector are not so large
        if ee:
            e = 1
            # scale for the cylinders of the end effector
            escale = cscale*0.9
            for k in np.flip(range(self.ndof)):
                if np.linalg.norm(self.Ts[k][:,3]-self.Ts[k-1][:,3])<1e-6:
                    e=e+1
                else:
                    break
                    
        # Cylinders that represent the joint directions (except for the last joint)
        resolution = 10
        plot_cylinder(ax, [0,0,0-cscale], [0,0,0+cscale], radius, resolution, alpha=1, color=color)
        for k in range(self.ndof-e):
            plot_cylinder(ax, 
                          [self.Ts[k][0,3]-cscale*self.Ts[k][0,2], 
                           self.Ts[k][1,3]-cscale*self.Ts[k][1,2],
                           self.Ts[k][2,3]-cscale*self.Ts[k][2,2]], 
                          [self.Ts[k][0,3]+cscale*self.Ts[k][0,2], 
                           self.Ts[k][1,3]+cscale*self.Ts[k][1,2], 
                           self.Ts[k][2,3]+cscale*self.Ts[k][2,2]], radius, resolution, alpha=1, color=color)
        for k in range(self.ndof-e, self.ndof):
            plot_cylinder(ax, 
                          [self.Ts[k][0,3]-escale*self.Ts[k][0,2], 
                           self.Ts[k][1,3]-escale*self.Ts[k][1,2],
                           self.Ts[k][2,3]-escale*self.Ts[k][2,2]], 
                          [self.Ts[k][0,3]+escale*self.Ts[k][0,2], 
                           self.Ts[k][1,3]+escale*self.Ts[k][1,2], 
                           self.Ts[k][2,3]+escale*self.Ts[k][2,2]], radius*0.7, resolution, alpha=1, color=color)
        # Plot an end effector
        if ee:
            # End effector (defined by 4 points)
            p1 = np.array([0, 0.1, 0, 1]); p2 = np.array([0, 0.1, 0.2, 1])
            p3 = np.array([0, -0.1, 0, 1]); p4 = np.array([0, -0.1, 0.2, 1])
            p1 = self.Ts[-1].dot(p1); p2 = self.Ts[-1].dot(p2); p3 = self.Ts[-1].dot(p3); p4 = self.Ts[-1].dot(p4)
            # Plot an end effector
            ax.plot([p1[0],p2[0]], [p1[1],p2[1]], [p1[2],p2[2]], color='k', linewidth=3)
            ax.plot([p3[0],p4[0]], [p3[1],p4[1]], [p3[2],p4[2]], color='k', linewidth=3)
            ax.plot([p1[0],p3[0]], [p1[1],p3[1]], [p1[2],p3[2]], color='k', linewidth=3)
        # Reference frame for the end effector (with respect to frame 0)
        ax.plot([self.Ts[-1][0,3],self.Ts[-1][0,3]+ascale*self.Ts[-1][0,0]], 
                [self.Ts[-1][1,3],self.Ts[-1][1,3]+ascale*self.Ts[-1][1,0]], 
                [self.Ts[-1][2,3],self.Ts[-1][2,3]+ascale*self.Ts[-1][2,0]], color='r')
        ax.plot([self.Ts[-1][0,3],self.Ts[-1][0,3]+ascale*self.Ts[-1][0,1]], 
                [self.Ts[-1][1,3],self.Ts[-1][1,3]+ascale*self.Ts[-1][1,1]], 
                [self.Ts[-1][2,3],self.Ts[-1][2,3]+ascale*self.Ts[-1][2,1]], color='g')
        ax.plot([self.Ts[-1][0,3],self.Ts[-1][0,3]+ascale*self.Ts[-1][0,2]], 
                [self.Ts[-1][1,3],self.Ts[-1][1,3]+ascale*self.Ts[-1][1,2]], 
                [self.Ts[-1][2,3],self.Ts[-1][2,3]+ascale*self.Ts[-1][2,2]], color='b')
        # Reference frame for the base (0)
        ax.plot([0,ascale], [0,0], [0,0], color='r')
        ax.plot([0,0], [0,ascale], [0,0], color='g')
        ax.plot([0,0], [0,0], [0,ascale], color='b')
        # Point of view
        ax.view_init(elev=elev, azim=azim)
        # Limits fot the figure axes
        if axlimits!=None:
            ax.set_xlim3d(axlimits[0][0], axlimits[0][1])
            ax.set_ylim3d(axlimits[1][0], axlimits[1][1])
            ax.set_zlim3d(axlimits[2][0], axlimits[2][1])
