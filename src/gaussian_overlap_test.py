#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

def gaussian(x, sigma):
    return np.exp(-x**2/(2*sigma**2))/(sigma*np.sqrt(2*(np.pi)))

# def overlap(f, g):

def main():
    sigma1 = 2
    sigma2 = sigma1
    del_x = 1
    x = np.linspace(-10, 10, 1000000); dL = (20/1000000.)
    
    f = gaussian(x, sigma1)
    g = gaussian(x-del_x, sigma2)
    overlap = np.minimum(f, g)

    print(sum(f)*dL)
    print(sum(g)*dL)
    print(sum(overlap)*dL)
    print(gaussian(del_x, sigma1))
    plt.plot(x, f)
    plt.plot(x, g)
    plt.plot(x, overlap)
    
    plt.show()

if __name__ == '__main__':
    main()