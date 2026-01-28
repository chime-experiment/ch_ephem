"""CHIME ephemeris routines.

Instrument Observers and General Ephemeris Routines
====================================================

Any ephemeris routine which needs to know the location of the
observer on the Earth are accessible via instrument Observer
instances importable from :py:meth:`ch_ephem.observers`:

    >>> from ch_ephem.observers import pathfinder
    >>> pathfinder.solar_transit(...)

The `Observer` instances provided by this module are subclassed from
`caput.astro.observer.Observer` and can be used as normal `caput` observers.
Location and geometry data for the instrument observers are defined
in the data file `instruments.yaml` provided as part of `ch_ephem`.

Celestial Intermediate Reference System
=======================================

The precession of the Earth's axis gives noticeable shifts in object
positions over the life time of CHIME. To minimise the effects of this we
need to be careful and consistent with our ephemeris calculations.
Historically Right Ascension has been given with respect to the Vernal
Equinox which has a significant (and unnecessary) precession in the origin of
the RA axis. To avoid this we use the new Celestial Intermediate Reference
System which does not suffer from this issue.

Practically this means that when calculating RA, DEC coordinates for a source
position at a *given time* you must be careful to obtain CIRS coordinates
(and not equinox based ones). Internally using `ch_ephem.coord.object_coords` does
exactly that for you, so for any lookup of coordinates you should use that on
your requested body.

Note that the actual coordinate positions of sources must be specified using
RA, DEC coordinates in ICRS (which is roughly equivalent to J2000). The
purpose of `object_coords` is to transform into new RA, DEC coordinates taking
into account the precession and nutation of the Earth's polar axis since
then.

These kind of coordinate issues are tricky, confusing and hard to debug years
later, so if you're unsure you are recommended to seek some advice.

Radio Source Catalogs
=====================

This package provides access to the standard radio source catalogue used
by CHIME, which is available as `ch_ephem.sources.source_dictionary`:

    >>> from ch_ephem.sources import source_dictionary
    >>> source_dictionary['CAS_A'].ra
    <Angle 23h 23m 27.94s>

The four standard CHIME cailbrators are also available by name:

    >>> from ch_ephem.sources import CasA, CygA, VirA, TauA
    >>> CasA.ra
    <Angle 23h 23m 27.94s>

Additional catalogues which are not part of the standard radio source catalogue
are also available through the `radiocosmology/fluxcat` package
<https://github.com/radiocosmology/fluxcat/>

Submodules
==========

Note: the submodules `observers` and `sources` read data from disk at import time.

.. autosummary::
    :toctree: _autosummary

    coord
    observers
    pointing
    sources
    time
"""

__all__ = ["__version__", "coord", "pointing", "time"]

# We deliberately do not import .observers and .sources here
# because they both perform reads from disk
from importlib.metadata import PackageNotFoundError, version

from . import coord, pointing, time

try:
    __version__ = version("ch_ephem")
except PackageNotFoundError:
    # package is not installed
    pass

del version, PackageNotFoundError
