import astropy.units as u
from radis import Spectrum
import numpy as np
w = np.linspace(200, 300, 1000) * u.nm
I = np.random.rand(len(w)) * u.mW/u.cm**2/u.sr/u.nm
s = Spectrum.from_array(w, I, 'radiance_noslit')
s.apply_slit(0.5, "nm",shape="gaussian")
s.plot("radiance_noslit")
s.plot("radiance")