#!/usr/bin/env python3
#encoding: UTF-8

class Charge:
    def __init__(self, *args, **kwargs):
        """Create charge object from """
        self.setvars(*args, **kwargs)

    def setvars(self, nr, charge, charge_diff=None):
        assert nr.shape == (3)
        self.nr = nr
        assert charge.shape == (nr[0], nr[1], nr[2])
        self.charge = charge
        if charge_diff != None:
            assert charge_diff.shape == (nr[0], nr[1], nr[2])
            self.charge_diff = charge_diff

        self.ax = None
        self.xcoords = None


if __name__ == "__main__":

    read_charge_file_hdf5(filename, nr)