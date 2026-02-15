import numpy as np 
import tkinter as tk    #tkinter needed for GUI
from tkinter import ttk
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation 
from mpl_toolkits.mplot3d import Axes3D         #for 3D plotting 

#physical constants
g = 9.81  #acceleration due to gravity (m/s^2)
e = 0.8   #coefficient of restitution (= energy loss on each bounce)

#BouncingBall class
#defining the physics and trajectory calculation of the ball
class BouncingBall:
    def __init__(self, H, vx, vy, vz):
        """
        Balls initial conditions here.

        Parameters:
        H: initial height (m)
        vx: velocity in x-direction (m/s)
        vy: velocity in y-direction (m/s)
        vz: initial vertical velocity (m/s)
        """
        self.H = H
        self.vx = vx
        self.vy = vy
        self.vz = vz

    def trajectory(self, dt=0.01, tmax=10):
        """
        Computing the 3D trajectory of the ball over time.

        Uses time-stepping (Euler integration) with
        gravity and inelastic bounces (as defined above) off the ground.

        Parameters:
        dt: time step (s)
        tmax: maximum simulation time (s), so extreme input doesn't run forever

        Returns:
        Arrays of x, y, z positions over time
        """
        #initial time and position
        t = 0.0
        x, y, z = 0.0, 0.0, self.H
        vz = self.vz

        #lists to store the trajectory
        xs, ys, zs = [], [], []

        #time integration loop
        while t < tmax:
            #updates vertical position and velocity using kinematics
            z_new = z + vz * dt - 0.5 * g * dt**2
            vz_new = vz - g * dt

            #handles collision with the ground
            if z_new <= 0 and vz < 0:
                z_new = 0
                vz_new = -e * vz_new  #reverse vertical velocity

            #update horizontal positions (with constant velocity)
            x += self.vx * dt
            y += self.vy * dt

            #commit updated values
            z = z_new
            vz = vz_new

            #store trajectory points
            xs.append(x)
            ys.append(y)
            zs.append(z)

            #stop simulation once motion is very small
            if abs(vz) < 0.05 and z < 0.01:
                break

            t += dt

        return np.array(xs), np.array(ys), np.array(zs)



#simulation + visualization
def run_simulation():
    """
    Reads user input from the GUI, runs the physics simulation, and displays a 3D animated drop. End values are also printed in the GUI.
    """
    #read initial conditions from GUI input
    H = float(entry_H.get())
    vx = float(entry_vx.get())
    vy = float(entry_vy.get())
    vz = float(entry_vz.get())

    #create ball and compute trajectory
    ball = BouncingBall(H, vx, vy, vz)
    x, y, z = ball.trajectory()

    #display final resting position
    landing_label.config(
        text=f"Final rest position: x = {x[-1]:.2f} m, y = {y[-1]:.2f} m"
    )

    #create 3D figure
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    #set axis limits and labels
    ax.set_xlim(min(x) - 1, max(x) + 1)
    ax.set_ylim(min(y) - 1, max(y) + 1)
    ax.set_zlim(0, max(z) + 1)
    ax.set_xlabel("X (m)")
    ax.set_ylabel("Y (m)")
    ax.set_zlabel("Z (m)")

    #initialize plot objects:
    #point: current ball position
    #line: full trajectory so far
    point, = ax.plot([], [], [], 'ro')
    line, = ax.plot([], [], [], 'b-', lw=2)

    def update(frame):
        """
        Animation update-function.
        Advances the ball position frame-by-frame.
        """
        #update ball position
        point.set_data([x[frame]], [y[frame]])
        point.set_3d_properties([z[frame]])

        #update trajectory line
        line.set_data(x[:frame], y[:frame])
        line.set_3d_properties(z[:frame])

        return point, line

    #create actual animation
    ani = FuncAnimation(fig, update, frames=len(x), interval=20)

    #display the plot window
    plt.show()



#GUI setup
root = tk.Tk()
root.title("3D Bouncing Ball Simulator")

#main container frame
frame = ttk.Frame(root, padding=10)
frame.grid()

#input: initial height
ttk.Label(frame, text="Initial height H (m):").grid(row=0, column=0)
entry_H = ttk.Entry(frame)
entry_H.insert(0, "10")
entry_H.grid(row=0, column=1)

#input: x velocity
ttk.Label(frame, text="vx (m/s):").grid(row=1, column=0)
entry_vx = ttk.Entry(frame)
entry_vx.insert(0, "2")
entry_vx.grid(row=1, column=1)

#input: y velocity
ttk.Label(frame, text="vy (m/s):").grid(row=2, column=0)
entry_vy = ttk.Entry(frame)
entry_vy.insert(0, "1")
entry_vy.grid(row=2, column=1)

#iput: vertical velocity
ttk.Label(frame, text="vz (m/s):").grid(row=3, column=0)
entry_vz = ttk.Entry(frame)
entry_vz.insert(0, "0")
entry_vz.grid(row=3, column=1)

#run button/command
ttk.Button(
    frame,
    text="Run Simulation",
    command=run_simulation
).grid(row=4, column=0, columnspan=2)

#output label for final position
landing_label = ttk.Label(frame, text="Final rest position: ")
landing_label.grid(row=5, column=0, columnspan=2)

#start GUI (our actual program)
root.mainloop()
