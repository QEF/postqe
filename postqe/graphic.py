import numpy as np
import matplotlib.pyplot as plt

#variables for plots and movie

frps = 25


def plot_sp(x, s, rgb):
    plt.plot(x, s, '-', color = rgb)
    plt.gcf().subplots_adjust(bottom=0.15)
    plt.fill_between(x, s, 0.3, color = rgb)
    plt.xlabel('wavelength [nm]', fontsize=45, color='black')
    plt.ylabel('absorption', fontsize=45, color='black')
    plt.axis([380, 780, 0.30, 25])
    plt.tick_params(axis='x', which = 'major', length=50, width=10, labelsize='40', top = 'off', bottom = 'off')
    plt.xticks([450,600,750])
    plt.yticks([])


def plot_movie(x, s, rgb):
    plt.plot(x, s, '-', color = rgb)
    plt.gcf().subplots_adjust(bottom=0.20)
    plt.fill_between(x, s, 0.3, color = rgb)
    plt.xlabel('wavelength [nm]', fontsize=40, color='black')
    plt.ylabel('absorption', fontsize=40, color='black')
    plt.axis([380, 780, 0.30, 25])
    plt.tick_params(axis='x', which = 'major', length=50, width=10, labelsize='35', top = 'off', bottom = 'off')
    plt.xticks([450,600,750])
    plt.yticks([])


def npoints(x, y, m, n):
    """
    Interpolation function for movie.

    :param x:
    :param y:
    :param m:
    :param n:
    :return:
    """
    m_out = list()
    nDim = n * (x-1)
    for i in np.arange(1,x):
        m_out.append(np.transpose(np.array([np.linspace(a, b, n) for a, b in zip(m[i-1,:], m[i,:])])))
    m_out = np.asarray(m_out).reshape(nDim,y)
    return m_out






