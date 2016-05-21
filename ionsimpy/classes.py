import h5py as _h5
import numpy as _np


class Particles(object):
    """
    Fundamental class containing particles.
    """
    def __init__(self, particles):
        self._particles = particles
        self._x   = particles[:, 0]
        self._xp  = particles[:, 1]
        self._y   = particles[:, 2]
        self._yp  = particles[:, 3]
        self._z   = particles[:, 4]
        self._zp  = particles[:, 5]

    @property
    def x(self):
        """
        Beam particles coordinates :math:`x`.
        """
        return self._particles.value[:, 0]

    @property
    def xp(self):
        """
        Beam particles coordinates :math:`x'`.
        """
        return self._particles.value[:, 1]

    @property
    def y(self):
        """
        Beam particles coordinates :math:`y`.
        """
        return self._particles.value[:, 2]

    @property
    def yp(self):
        """
        Beam particles coordinates :math:`y'`.
        """
        return self._particles.value[:, 3]
        # return self._yp

    @property
    def z(self):
        """
        Beam particles coordinates :math:`z`.
        """
        return self._particles.value[:, 4]
        # return self._z

    @property
    def zp(self):
        """
        Beam particles coordinates :math:`z'`.
        """
        return self._particles.value[:, 5]
        # return self._zp


class Ions(object):
    """
    A class for loading ion particles.

    Parameters
    ----------

    filename : str
        Filename holding ion particles.
    """
    def __init__(self, step_obj):
        self._step_obj = step_obj

        for key in step_obj.keys():
            setattr(self, key, I_Step(step_obj[key]))


class Beam_Electrons(Particles):
    """
    A class for loading beam electrons.
    """
    def __init__(self, beam_obj):
        self._particles = beam_obj
        super().__init__(self._particles)

    @property
    def particles(self):
        return self._particles


class Field(object):
    """
    A class for handling fields.
    """
    def __init__(self, field_obj):
        self._Ex = _np.transpose(field_obj['Ex'].value)
        self._Ey = _np.transpose(field_obj['Ey'].value)
        self._Ez = _np.transpose(field_obj['Ez'].value)

        for key in field_obj.attrs.keys():
            setattr(self, key, field_obj.attrs[key])

    def Ex(self, index):
        """
        The ::math::`$E_x$` component of the electric field ::math::`$\vec{E}$`.
        """
        return self._Ex[index, :, :]

    def Ey(self, index):
        """
        The ::math::`$E_y$` component of the electric field ::math::`$\vec{E}$`.
        """
        return self._Ey[index, :, :]

    def Ez(self, index):
        """
        The ::math::`$E_z$` component of the electric field ::math::`$\vec{E}$`.
        """
        return self._Ez[index, :, :]

    # @property
    # def x_grid(self):
    #     """
    #     The coordinates for each grid point in the field ::math::`$E_x`.
    #     """
    #     return self._x_grid

    # @property
    # def y_grid(self):
    #     """
    #     The coordinates for each grid point in the field ::math::`$E_y`.
    #     """
    #     return self._y_grid

    # @property
    # def z_grid(self):
    #     """
    #     The coordinates for each grid point in the field ::math::`$E_z`.
    #     """
    #     return self._z_grid


class I_Step(Particles):
    """
    Class for each step in ion motion integration.
    """
    def __init__(self, step_obj):
        super().__init__(step_obj)
        self._step_obj = step_obj

        self._step = step_obj.attrs['Step']

    @property
    def particles(self):
        return self._particles

    @property
    def step(self):
        return self._step


class E_Step(object):
    """
    Class for each step in electron motion integration.
    """
    def __init__(self, step_obj):
        self._step_obj = step_obj

        for key in step_obj.keys():
            if key == 'field':
                setattr(self, 'field', Field(step_obj['field']))
            elif key == 'electrons':
                setattr(self, 'electrons', Beam_Electrons(step_obj['electrons']))
            elif key == 'ions_steps':
                setattr(self, 'ions_steps', Ions(step_obj['ions_steps']))


class Sim(object):
    """
    Container for simulation.
    """
    def __init__(self, filename):

        self.file = _h5.File(filename, 'r')

        for key in self.file.keys():
            setattr(self, key, E_Step(self.file[key]))

        for key in self.file.attrs.keys():
            setattr(self, key, self.file.attrs[key])
