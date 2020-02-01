class ExpOptions:

    def __init__(self, field: float, slC: float):
        """
        :param field: static magnetic field in mT, needed to find omega_e and _H
        :param slC: unit is M, spin label concentration for scaling relaxations to get "relaxivities"
        """
        self.field = field

        self.slC = slC

        self.smaxMod = 'tethered'  # either 'tethered' or 'free'

        self.ksig_bulk = 95.4  # unit is s^-1 M^-1
        # The only place I can find this is Franck, JM, et. al.; "Anomalously Rapid
        # Hydration Water Diffusion Dynamics Near DNA Surfaces" J. Am. Chem. Soc.
        # 2015, 137, 12013−12023. Figure 3 caption

        self.T100 = 2.5
        # this is the T1 without spin label and without mw power, unit is sec

        self.tcorr_bulk = 54  # (section 2.5), "corrected" bulk tcorr, unit is ps

        self.dH2O = 2.3e-9
        # (Eq. 19-20) bulk water diffusivity, unit is d^2/s where d is distance in
        # meters. *didnt use m to avoid confusion with mass

        self.dSL = 4.1e-10
        # (Eq. 19-20) spin label diffusivity, unit is d^2/s where d is distance in
        # meters. *didnt use m to avoid confusion with mass

        self.k_low_bulk = 366  # unit is s^-1 M^-1
        # The only place I can find this is Franck, JM, et. al.; "Anomalously Rapid
        # Hydration Water Diffusion Dynamics Near DNA Surfaces" J. Am. Chem. Soc.
        # 2015, 137, 12013−12023. Figure 3 caption

    @property
    def smaxMod(self):
        return self.smaxMod

    @smaxMod.setter
    def smaxMod(self, value: str):
        if value == 'tethered':
            self.smaxMod = value
        elif value == 'free':
            self.smaxMod = value
        else:
            raise ValueError('smaxMod should be either `tethered` or `free`')

    def print(self):
        # TODO: implement this
        raise NotImplementedError

