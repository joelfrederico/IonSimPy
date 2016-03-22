import h5py as _h5
import numpy as _np


class Particles(object):
    def __init__(self, particles):
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
        return self._x

    @property
    def xp(self):
        """
        Beam particles coordinates :math:`x'`.
        """
        return self._xp

    @property
    def y(self):
        """
        Beam particles coordinates :math:`y`.
        """
        return self._y

    @property
    def yp(self):
        """
        Beam particles coordinates :math:`y'`.
        """
        return self._yp

    @property
    def z(self):
        """
        Beam particles coordinates :math:`z`.
        """
        return self._z

    @property
    def zp(self):
        """
        Beam particles coordinates :math:`z'`.
        """
        return self._zp


class Ions(Particles):
    """
    A class for loading electron beam particles.

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


class Field(object):
    """
    A class for handling fields.
    """
    def __init__(self, field_obj):
        self._Ex = field_obj['Ex'].value

    @property
    def Ex(self):
        return self._Ex


class Step(object):
    def __init__(self, step_obj):
        for key in step_obj.keys():
            if key == 'field':
                setattr(self, 'field', Field(step_obj['field']))


class Sim(object):
    def __init__(self, filename):

        self.file = _h5.File(filename)

        for key in self.file.keys():
            setattr(self, key, Step(self.file[key]))

        for key in self.file.attrs.keys():
            setattr(self, key, self.file.attrs[key][0])
