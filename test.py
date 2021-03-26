import matplotlib.pyplot as plt
import numpy as np

fig1, ax1 = plt.subplots()
fig2, ax2 = plt.subplots()

x=np.linspace(-10,10)
y1=x**2
y2 = x**3

ax1.plot(x, y1)


ax2.plot(x,y2)
plt.show()
