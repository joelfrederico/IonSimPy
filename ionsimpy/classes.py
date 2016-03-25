import h5py as _h5
import numpy as _np


class Particles(object):
    def __init__(self, particles):
        self._particles = particles
        # self._x   = particles[:, 0]
        # self._xp  = particles[:, 1]
        # self._y   = particles[:, 2]
        # self._yp  = particles[:, 3]
        # self._z   = particles[:, 4]
        # self._zp  = particles[:, 5]

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


class Ions(Particles):
    """
    A class for loading ion particles.

    Parameters
    ----------

    filename : str
        Filename holding ion particles.
    """
    def __init__(self, file):
        self.file       = file
        self._ion_group = self.file['ions']
        n_keys          = len(self._ion_group)
        dataset_list = list(self._ion_group)
        n_parts = self._ion_group[dataset_list[0]].shape[0]
        self._data = _np.empty((n_keys, n_parts, 6), dtype=float)
        for i in range(n_keys):
            self._data[i, :, :] = self._ion_group['step_{:03d}'.format(i)].value

    @property
    def x(self):
        return self._data[:, :, 0]


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

    @property
    def Ex(self):
        return self._Ex

    @property
    def Ey(self):
        return self._Ey


class Step(object):
    def __init__(self, step_obj):
        self._step_obj = step_obj

        for key in step_obj.keys():
            if key == 'field':
                setattr(self, 'field', Field(step_obj['field']))
            elif key == 'ebeam':
                setattr(self, 'ebeam', Beam_Electrons(step_obj['ebeam']))
            elif key == 'ions':
                setattr(self, 'ions', Beam_Electrons(step_obj['ions']))


class Sim(object):
    def __init__(self, filename):

        self.file = _h5.File(filename)

        for key in self.file.keys():
            setattr(self, key, Step(self.file[key]))

        for key in self.file.attrs.keys():
            setattr(self, key, self.file.attrs[key][0])
