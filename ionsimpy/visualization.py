from . import classes as _cl

import numpy as _np
import scisalt.matplotlib as _sm
import matplotlib.pyplot as _plt
import scipy.special as _spp
import scipy.constants as _spc
from matplotlib.gridspec import GridSpec as _GridSpec
import matplotlib.widgets as _wid


class readfield(object):
    def __init__(self, filename):
        # ================================
        # Setup important vars
        # ================================
        self._filename  = filename
        self.sim        = _cl.Sim(filename)
        self.step1      = self.sim.Step_0000
        self.n_field    = self.sim.n_field_z
        self.q          = self.sim.q_tot * _spc.elementary_charge / ( _spc.epsilon_0 * self.sim.sz )
        self._clear     = _np.array([], dtype=object)
        self._firstpass = True

        # ================================
        # Setup plot
        # ================================
        # fig, ax = _sm.setup_axes(rows=4, cols=3, figsize=(4, 12), expand=False)

        self.smin = 0
        self.smax = self.n_field-1

        savg = _np.round((self.smin+self.smax)/2)
    
        rows         = 4
        cols         = 3
        self.fig     = _plt.figure(figsize=(8, 11))
        self.ax      = _np.empty((rows, cols), dtype=object)
        self.imshows = _np.empty((rows, cols), dtype=object)

        print(self.fig.get_size_inches())
        gs = _GridSpec(rows, cols, top=0.97, bottom=0.13, left=0.05, right=0.95, hspace=0.25, wspace=0.25)
        # gs = _GridSpec(rows, cols)
        
        for row in range(rows):
            for col in range(cols):
                self.ax[row, col] = self.fig.add_subplot(gs[row, col])
        
        gs_gui = _GridSpec(1, 1, top=0.08, bottom=0.02)
        ax_gui = self.fig.add_subplot(gs_gui[0, 0])

        self.slider = _wid.Slider(ax_gui, "E-field\nSlice", self.smin, self.smax, savg, dragging=False, valfmt='%d')

        self.slider.on_changed(self.update)
        self.fig.canvas.mpl_connect('key_press_event', self.press)
        self.fig.canvas.mpl_connect('scroll_event', self.scroll)
        _plt.show()

    def update(self, index):
        index = _np.round(index)
        if index < self.smin:
            index = self.smin
            self.slider.set_val(index)
        elif index > self.smax:
            index = self.smax
            self.slider.set_val(index)

        index = _np.int(index)
    
        self.update_plot(index)

    def press(self, event):
        if event.key == 'right':
            self.slider.set_val(self.slider.val+1)
        elif event.key == 'left':
            self.slider.set_val(self.slider.val-1)
    
    def scroll(self, event):
        self.slider.set_val(self.slider.val + event.step/60)

    def update_plot(self, index):
        firstpass = self._firstpass
        if self._firstpass:
            self._firstpass = False

        field = self.step1.field
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
        
        if firstpass:
            row = 0
            col = 0
            axc = self.ax[row, col]
            self.imshows[row, col] = _sm.imshow(Dat_Em, ax=axc, **genkwargs)
            _sm.addlabel(ax=axc, toplabel="Data: Em")

            col = 1

            axc = self.ax[row, col]

            self.imshows[row, col] = _sm.imshow(Dat_Ex, ax=axc, **rbkwargs)
            # _sm.contour(Dat_Ex, ax=axc, **rbkwargs)
            _sm.addlabel(ax=axc, toplabel="Data: Ex")
            
            col = 2

            axc = self.ax[row, col]
            
            self.imshows[row, col] = _sm.imshow(Dat_Ey, ax=axc, **rbkwargs)
            _sm.addlabel(ax=axc, toplabel="Data: Ey")
        else:
            row = 0
            self.imshows[row, 0].im.set_data(_np.transpose(Dat_Em))
            self.imshows[row, 1].im.set_data(_np.transpose(Dat_Ex))
            self.imshows[row, 2].im.set_data(_np.transpose(Dat_Ey))
        
        # ================================
        # Basetti-Erskine
        # ================================
        
        if firstpass:
            # sr = _np.mean([2.42769e-6, 2.44567e-6])
            sr = 2.4276628847185805e-06
    
            sx = 2.44568e-6
            sy = 2.42769e-6
            # var_x = sx**2
            # var_y = sy**2
            # var_x_minus_var_y = var_x-var_y
            
            self.n_pts = Dat_Ex.shape[0]
            
            xvals = field.x_grid
            yvals = field.y_grid
            
            self.BE_Emag = _np.empty((self.n_pts, self.n_pts))
            self.BE_Ex_plot = _np.empty((self.n_pts, self.n_pts))
            self.BE_Ey_plot = _np.empty((self.n_pts, self.n_pts))
            
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
                    E = (self.q/(2*_np.sqrt(2*_np.pi*(sx**2-sy**2)))) * (_spp.wofz((x_in + 1j*y_in)/fval) - _np.exp(-x_in**2/(2*sx**2)-y_in**2/(2*sy**2)) * _spp.wofz(aribr))
                    Ex = _np.imag(E) * _np.sign(x)
                    Ey = _np.real(E) * _np.sign(y)
            
                    self.BE_Ex_plot[i, j] = Ex
                    self.BE_Ey_plot[i, j] = Ey
                    self.BE_Emag[i, j]    = _np.absolute(E)
            
            temp            = self.BE_Ex_plot.transpose()
            self.BE_Ey_plot = self.BE_Ey_plot.transpose()
            self.BE_Ey_plot      = temp
            self.BE_Emag         = self.BE_Emag.transpose()
        
            row = 1
            col = 0
            axc = self.ax[row, col]
            self.imshows[row, col] = _sm.imshow(self.BE_Emag, ax=axc, **genkwargs)
            _sm.addlabel(ax=axc, toplabel="B/E: E_mag")
            
            col = 1
            axc = self.ax[row, col]
            self.imshows[row, col] = _sm.imshow(self.BE_Ex_plot, ax=axc, **rbkwargs)
            _sm.addlabel(ax=axc, toplabel="B/E: Ex")
            
            col = 2
            axc = self.ax[row, col]
            self.imshows[row, col] = _sm.imshow(self.BE_Ey_plot, ax=axc, **rbkwargs)
            _sm.addlabel(ax=axc, toplabel="B/E: Ey")
        # else:
            # self.imshows[row, 0].im.set_data(_np.transpose(BE_Emag))
            # self.imshows[row, 1].im.set_data(_np.transpose(BE_Ex_plot))
            # self.imshows[row, 2].im.set_data(_np.transpose(BE_Ey_plot))
        
        # ================================
        # Difference
        # ================================
        row = 2
        col = 0
        axc = self.ax[row, col]
        toplot = (Dat_Em-self.BE_Emag)/Dat_Em
        toplot[_np.int((self.n_pts-1)/2), :] = 0
        vmag = _np.max(_np.abs(toplot))
        if firstpass:
            self.imshows[row, col] = _sm.imshow(toplot, ax=axc, vmin=-vmag, vmax=vmag, **rbkwargs)
            _sm.addlabel(ax=axc, toplabel="% Diff: E_mag")
        else:
            self.imshows[row, col].im.set_data(_np.transpose(toplot))
            self.imshows[row, col].im.set_clim(vmin=-vmag, vmax=vmag)
        
        col = 1
        toplot = (Dat_Ex-self.BE_Ex_plot)/Dat_Ex
        midpt = _np.int((self.n_pts-1)/2)
        toplot[midpt-2:midpt+3, :] = 0
        vmag = _np.max(_np.sqrt(toplot**2))
        axc = self.ax[row, col]
        if firstpass:
            self.imshows[row, col] = _sm.imshow(toplot, ax=axc, vmin=-vmag, vmax=vmag, **rbkwargs)
            _sm.addlabel(ax=axc, toplabel="% Diff: Ex")
        else:
            self.imshows[row, col].im.set_data(_np.transpose(toplot))
            self.imshows[row, col].im.set_clim(vmin=-vmag, vmax=vmag)
        
        col = 2
        toplot = (Dat_Ey-self.BE_Ey_plot)/Dat_Ey
        midpt = _np.int((self.n_pts-1)/2)
        toplot[:, midpt-2:midpt+3] = 0
        vmag = _np.max(_np.sqrt(toplot**2))
        axc = self.ax[row, col]
        if firstpass:
            self.imshows[row, col] = _sm.imshow(toplot, ax=axc, vmin=-vmag, vmax=vmag, **rbkwargs)
            _sm.addlabel(ax=axc, toplabel="% Diff: Ey")
        else:
            self.imshows[row, col].im.set_data(_np.transpose(toplot))
            self.imshows[row, col].im.set_clim(vmin=-vmag, vmax=vmag)
        
        # ================================
        # Radial fields
        # ================================
        
        if firstpass:
            # sr = _np.mean([sx, sy])
            sr = _np.mean([2.42769e-6, 2.44567e-6])
            # sr = 1
            
            xvals = field.x_grid
            yvals = field.y_grid
            
            Emag = _np.empty((self.n_pts, self.n_pts))
            Ex_plot = _np.empty((self.n_pts, self.n_pts))
            Ey_plot = _np.empty((self.n_pts, self.n_pts))
            
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
                        Em = self.q * (1-_np.exp(-rsq/(2*sr**2)))/(2*_np.pi*r)
                    Emag[i, j] = Em
                    Ex_plot[i, j] = Em * x/r
                    Ey_plot[i, j] = Em * y/r
            
            axc = self.ax[3, 0]
            _sm.imshow(Emag, ax=axc, **genkwargs)
            _sm.addlabel(ax=axc, toplabel="Radial Sym.: E_mag")
            
            axc = self.ax[3, 1]
            _sm.imshow(Ex_plot, ax=axc, **rbkwargs)
            _sm.addlabel(ax=axc, toplabel="Radial Sym.: Ex")
            
            axc = self.ax[3, 2]
            _sm.imshow(Ey_plot, ax=axc, **rbkwargs)
            _sm.addlabel(ax=axc, toplabel="Radial Sym.: Ey")
