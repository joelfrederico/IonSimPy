import h5py as _h5


class Ebeam(object):
    """
    A class for loading electron beam particles.

    Parameters
    ----------

    filename : str
        Filename holding particles.
    """
    def __init__(self, filename):
        self._filename = filename

        with _h5.File(filename) as f:
            particles = f['ebeam']['particles'].value
            
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
