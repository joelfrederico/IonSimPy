from . import classes as _cl

import numpy as _np
import scisalt.matplotlib as _sm
import matplotlib.pyplot as _plt
import scipy.special as _spp
import scipy.constants as _spc
from matplotlib.gridspec import GridSpec as _GridSpec
import matplotlib.widgets as _wid


def readfield(filename):
    # filename = 'large.h5'
    # filename = 'output.h5'
    
    sim    = _cl.Sim(filename)
    step1 = sim.Step_0000
    n_field = sim.n_field_z

    q = sim.q_tot * _spc.elementary_charge / ( _spc.epsilon_0 * sim.sz )

    # ================================
    # Setup plot
    # ================================
    # fig, ax = _sm.setup_axes(rows=4, cols=3, figsize=(4, 12), expand=False)
    rows = 4
    cols = 3
    ax = _np.empty((rows, cols), dtype=object)
    fig = _plt.figure(figsize=(8, 11))
    print(fig.get_size_inches())
    gs = _GridSpec(rows, cols, top=0.97, bottom=0.13, left=0.05, right=0.95, hspace=0.25, wspace=0.25)
    # gs = _GridSpec(rows, cols)
    
    for row in range(rows):
        for col in range(cols):
            ax[row, col] = fig.add_subplot(gs[row, col])
    
    gs_gui = _GridSpec(1, 1, top=0.08, bottom=0.02)
    ax_gui = fig.add_subplot(gs_gui[0, 0])
    smin = 0
    smax = n_field-1
    savg = _np.round((smin+smax)/2)
    slider = _wid.Slider(ax_gui, "E-field\nSlice", smin, smax, savg, dragging=False, valfmt='%d')

    def update(index):
        index = _np.round(index)
        if index < smin:
            index = smin
            slider.set_val(index)
        elif index > smax:
            index = smax
            slider.set_val(index)

        index = _np.int(index)
    
        update_plot(ax, q, step1, index)

    def press(event):
        if event.key == 'right':
            slider.set_val(slider.val+1)
        elif event.key == 'left':
            slider.set_val(slider.val-1)
    
    def scroll(event):
        slider.set_val(slider.val + event.step/60)
    
    slider.on_changed(update)
    fig.canvas.mpl_connect('key_press_event', press)
    fig.canvas.mpl_connect('scroll_event', scroll)

    _plt.show()

    return ax


def update_plot(ax, q, step1, index):
    field = step1.field
    Dat_Ex = field.Ex(index)
    Dat_Ey = field.Ey(index)
    # ================================
    # Plot Data
    # ================================
    vmag = 0.1
    extent = _np.array([field.x_grid[0], field.x_grid[-1], field.y_grid[0], field.y_grid[-1]]) * 1e6

    # rbkwargs = {"vmin": -vmag, "vmax": vmag, "cmap": "RdBu"}
    rbkwargs = {"cmap": "RdBu", "add_cbar": True, "extent": extent}
    genkwargs = {"add_cbar": False, "extent": extent}
    
    Dat_Em = _np.sqrt(Dat_Ex**2+Dat_Ey**2)
    # Dat_Em = Dat_Ex
    axc = ax[0, 0]
    axc.clear()
    
    _sm.imshow(Dat_Em, ax=axc, **genkwargs)
    _sm.addlabel(ax=axc, toplabel="Data: Em")
    
    axc = ax[0, 1]
    axc.clear()
    
    _sm.imshow(Dat_Ex, ax=axc, **rbkwargs)
    # _sm.contour(Dat_Ex, ax=axc, **rbkwargs)
    _sm.addlabel(ax=axc, toplabel="Data: Ex")
    
    axc = ax[0, 2]
    axc.clear()
    
    _sm.imshow(Dat_Ey, ax=axc, **rbkwargs)
    _sm.addlabel(ax=axc, toplabel="Data: Ey")
    
    # ================================
    # Basetti-Erskine
    # ================================
    
    # sr = _np.mean([2.42769e-6, 2.44567e-6])
    sr = 2.4276628847185805e-06

    sx = 2.44568e-6
    sy = 2.42769e-6
    # var_x = sx**2
    # var_y = sy**2
    # var_x_minus_var_y = var_x-var_y
    
    n_pts = Dat_Ex.shape[0]
    
    xvals = field.x_grid
    yvals = field.y_grid
    
    BE_Emag = _np.empty((n_pts, n_pts))
    BE_Ex_plot = _np.empty((n_pts, n_pts))
    BE_Ey_plot = _np.empty((n_pts, n_pts))
    
    for i, x in enumerate(xvals):
        for j, y in enumerate(yvals):
            # x_in = _np.abs(x)
            fval = _np.sqrt(2*(sx**2-sy**2))
            x_in = _np.abs(x)
            y_in = _np.abs(y)
            r = sy/sx
            a = x_in / fval
            b = y_in / fval
            # aib = a + 1j*b
            aribr = a*r + 1j*b/r
            E = (q/(2*_np.sqrt(2*_np.pi*(sx**2-sy**2)))) * (_spp.wofz((x_in + 1j*y_in)/fval) - _np.exp(-x_in**2/(2*sx**2)-y_in**2/(2*sy**2)) * _spp.wofz(aribr))
            Ex = _np.imag(E) * _np.sign(x)
            Ey = _np.real(E) * _np.sign(y)
    
            BE_Ex_plot[i, j] = Ex
            BE_Ey_plot[i, j] = Ey
            BE_Emag[i, j]    = _np.absolute(E)
    
    temp       = BE_Ex_plot.transpose()
    BE_Ex_plot = BE_Ey_plot.transpose()
    BE_Ey_plot = temp
    BE_Emag    = BE_Emag.transpose()
    
    axc = ax[1, 0]
    axc.clear()
    _sm.imshow(BE_Emag, ax=axc, **genkwargs)
    _sm.addlabel(ax=axc, toplabel="B/E: E_mag")
    
    axc = ax[1, 1]
    axc.clear()
    _sm.imshow(BE_Ex_plot, ax=axc, **rbkwargs)
    _sm.addlabel(ax=axc, toplabel="B/E: Ex")
    
    axc = ax[1, 2]
    axc.clear()
    _sm.imshow(BE_Ey_plot, ax=axc, **rbkwargs)
    _sm.addlabel(ax=axc, toplabel="B/E: Ey")
    
    # ================================
    # Difference
    # ================================
    axc = ax[2, 0]
    axc.clear()
    toplot = (Dat_Em-BE_Emag)/Dat_Em
    toplot[_np.int((n_pts-1)/2), :] = 0
    vmag = _np.max(_np.abs(toplot))
    _sm.imshow(toplot, ax=axc, vmin=-vmag, vmax=vmag, **rbkwargs)
    _sm.addlabel(ax=axc, toplabel="% Diff: E_mag")
    
    toplot = (Dat_Ex-BE_Ex_plot)/Dat_Ex
    toplot[_np.int((n_pts-1)/2), :] = 0
    vmag = _np.max(_np.sqrt(toplot**2))
    axc = ax[2, 1]
    axc.clear()
    _sm.imshow(toplot, ax=axc, vmin=-vmag, vmax=vmag, **rbkwargs)
    _sm.addlabel(ax=axc, toplabel="% Diff: Ex")
    
    toplot = (Dat_Ey-BE_Ey_plot)/Dat_Ey
    toplot[:, _np.int((n_pts-1)/2)] = 0
    vmag = _np.max(_np.sqrt(toplot**2))
    axc = ax[2, 2]
    axc.clear()
    _sm.imshow(toplot, ax=axc, vmin=-vmag, vmax=vmag, **rbkwargs)
    _sm.addlabel(ax=axc, toplabel="% Diff: Ey")
    
    # ================================
    # Radial fields
    # ================================
    
    # sr = _np.mean([sx, sy])
    sr = _np.mean([2.42769e-6, 2.44567e-6])
    # sr = 1
    
    xvals = field.x_grid
    yvals = field.y_grid
    
    Emag = _np.empty((n_pts, n_pts))
    Ex_plot = _np.empty((n_pts, n_pts))
    Ey_plot = _np.empty((n_pts, n_pts))
    
    for i, x in enumerate(xvals):
        for j, y in enumerate(yvals):
            # E = 
            # Ex = _np.imag(E)
            # Ey = _np.real(E)
    
            rsq = x**2+y**2
            r = _np.sqrt(rsq)
            if rsq == 0:
                Em = 0
            else:
                Em = q * (1-_np.exp(-rsq/(2*sr**2)))/(2*_np.pi*r)
            Emag[i, j] = Em
            Ex_plot[i, j] = Em * x/r
            Ey_plot[i, j] = Em * y/r
    
    axc = ax[3, 0]
    axc.clear()
    _sm.imshow(Emag, ax=axc, **genkwargs)
    _sm.addlabel(ax=axc, toplabel="Radial Sym.: E_mag")
    
    axc = ax[3, 1]
    axc.clear()
    _sm.imshow(Ex_plot, ax=axc, **rbkwargs)
    _sm.addlabel(ax=axc, toplabel="Radial Sym.: Ex")
    
    axc = ax[3, 2]
    axc.clear()
    _sm.imshow(Ey_plot, ax=axc, **rbkwargs)
    _sm.addlabel(ax=axc, toplabel="Radial Sym.: Ey")
