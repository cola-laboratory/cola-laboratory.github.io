import matplotlib.pyplot as plt
import numpy as np
import sys
from matplotlib.pyplot import show, plot

n   = int(input())
fig = plt.figure()
ax  = fig.add_subplot(111)

while (n != -1):
    gen = int(input())
    ax.cla()
    ax.set_title("generation:" + str(gen))
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_xlim(0, 3)
    ax.set_ylim(0, 3)
    while (n>0):
        xs = float(input())
        ys = float(input())
        n = n - 1
        ax.scatter(xs, ys)
    plt.draw()
    plt.show(block = False)   # this creates an empty frozen window.
    n = int(input())
    if (n == -1):
        plt.show(block = True)
plt.close()


