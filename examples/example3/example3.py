#!/usr/bin/env python3
#encoding: UTF-8

"""
This is an example showing how to compute the electronic density of states (DOS or edos).

First it computes the DOS for silicon (input file Si.xml), then for magnetic Ni (lsda) showing how to
obtain the dos for spin up and down and how to plot it.
"""
    
if __name__ == "__main__":

    from postqe import compute_dos, simple_plot_xy, multiple_plot_xy

    e, dos_up, dos_down = compute_dos('Si.xml', filedos='filedosSi', e_min=-10, e_max=20,
                                      e_step=0.01, degauss=0.02, ngauss=0)

    # plot the DOS
    fig1 = simple_plot_xy(e,dos_up,xlabel="E (eV/cell)",ylabel="DOS (states/eV/cell)")
    fig1.savefig("figure_Sidos.png")

    e2, dos_up2, dos_down2 = compute_dos('Ni.xml', filedos='filedosNi', e_min=0, e_max=50,
                                      e_step=0.01, degauss=0.02, ngauss=0)

    # plot the DOS
    fig2 = simple_plot_xy(e2,dos_up2,xlabel="E (eV/cell)",ylabel="DOS (states/eV/cell)")
    fig2.savefig("figure_Nidosup.png")
    fig3 = simple_plot_xy(e2,dos_down2,xlabel="E (eV/cell)",ylabel="DOS (states/eV/cell)")
    fig3.savefig("figure_Nidosdown.png")

    # create a numpy matrix *doss* for plotting both spin up and down on the same plot
    import numpy as np
    doss = np.zeros((len(e2),2))
    doss[:,0] = dos_up2
    doss[:,1] = -dos_down2
    fig4 = multiple_plot_xy(e2,doss,xlabel="E (eV/cell)",ylabel="DOS (states/eV/cell)")
    fig4.savefig("figure_Nidosupanddown.png")
