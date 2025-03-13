"""CHIME co-ordinate transformations.

Functions
=========

- :py:meth:`cirs_radec`
- :py:meth:`star_cirs`
- :py:meth:`object_coords`
- :py:meth:`hadec_to_bmxy`
- :py:meth:`bmxy_to_hadec`
- :py:meth:`range_rate`
- :py:meth:`peak_ra`

"""

from __future__ import annotations

import warnings
from typing import TYPE_CHECKING

import numpy as np
import skyfield.starlib
from caput.time import skyfield_wrapper, unix_to_skyfield_time

if TYPE_CHECKING:
    import caput.time
    import skyfield.jpllib
    import skyfield.timelib
    import skyfield.units
    import skyfield.vectorlib

    SkyfieldSource = (
        skyfield.starlib.Star
        | skyfield.vectorlib.VectorSum
        | skyfield.jpllib.ChebyshevPosition
    )

del TYPE_CHECKING


def cirs_radec(
    body: skyfield.starlib.Star, obs: caput.time.Observer | None = None
) -> skyfield.starlib.Star:
    """Convert CIRS to ICRS.

    Converts a Skyfield body in CIRS coordinates at a given epoch to
    ICRS coordinates observed from `obs`.

    Parameters
    ----------
    body : skyfield.starlib.Star
        Skyfield Star object with positions in CIRS coordinates.
    obs : `caput.time.Observer` or None
        The observer instance to use. If not supplied use `chime`.

    Returns
    -------
    new_body : skyfield.starlib.Star
        Skyfield Star object with positions in ICRS coordinates
    """
    import skyfield.functions
    from skyfield.units import Angle

    if obs is None:
        from .observers import chime as obs

    ts = skyfield_wrapper.timescale

    epoch = ts.tt_jd(np.median(body.epoch))

    pos = obs.skyfield_obs().at(epoch).observe(body)

    # Matrix CT transforms from CIRS to ICRF (https://rhodesmill.org/skyfield/time.html)
    r_au, dec, ra = skyfield.functions.to_polar(
        np.einsum("ij...,j...->i...", epoch.CT, pos.position.au)
    )

    return skyfield.starlib.Star(
        ra=Angle(radians=ra, preference="hours"), dec=Angle(radians=dec), epoch=epoch
    )


def star_cirs(
    ra: skyfield.units.Angle,
    dec: skyfield.units.Angle,
    epoch: skyfield.timelib.Time,
    obs: caput.time.Observer | None = None,
) -> skyfield.starlib.Star:
    """Create wrapper for `skyfield.starlib.Star`.

    Creates a position given CIRS coordinates observed from `obs`.

    Parameters
    ----------
    ra, dec : skyfield.api.Angle
        RA and dec of the source in CIRS coordinates
    epoch : skyfield.api.Time
        Time of the observation
    obs : `caput.time.Observer` or None
        The observer instance to use. If not supplied use `chime`.

    Returns
    -------
    body : skyfield.starlib.Star
        Star object in ICRS coordinates
    """
    return cirs_radec(skyfield.starlib.Star(ra=ra, dec=dec, epoch=epoch), obs=obs)


def object_coords(
    body: SkyfieldSource,
    date: float | None = None,
    deg: bool = False,
    obs: caput.time.Observer | None = None,
) -> tuple[float, float]:
    """Calculate the RA and DEC of a source.

    Gives the ICRS coordinates if no date is given (=J2000), or if a date is
    specified gives the CIRS coordinates at that epoch.

    This also returns the *apparent* position, including abberation and
    deflection by gravitational lensing. This shifts the positions by up to
    20 arcseconds.

    Parameters
    ----------
    body : skyfield source
        skyfield.starlib.Star or skyfield.vectorlib.VectorSum or
        skyfield.jpllib.ChebyshevPosition body representing the source.
    date : float
        Unix time at which to determine ra of source If None, use Jan 01
        2000.
    deg : bool
        Return RA ascension in degrees if True, radians if false (default).
    obs : `caput.time.Observer` or None
        An observer instance to use. If not supplied use `chime`. For many
        calculations changing from this default will make little difference.

    Returns
    -------
    ra, dec: float
        Position of the source.
    """
    if obs is None:
        from .observers import chime as obs

    if date is None:  # No date, get ICRS coords
        if isinstance(body, skyfield.starlib.Star):
            ra, dec = body.ra.radians, body.dec.radians
        else:
            raise ValueError(
                "Body is not fixed, cannot calculate coordinates without a date."
            )

    else:  # Calculate CIRS position with all corrections
        date = unix_to_skyfield_time(date)
        radec = obs.skyfield_obs().at(date).observe(body).apparent().cirs_radec(date)

        ra, dec = radec[0].radians, radec[1].radians

    # If requested, convert to degrees
    if deg:
        ra = np.degrees(ra)
        dec = np.degrees(dec)

    # Return
    return ra, dec


def hadec_to_bmxy(ha_cirs, dec_cirs):
    """Convert CIRS hour angle and declination to CHIME/FRB beam-model XY coordinates.

    Parameters
    ----------
    ha_cirs : array_like
        The CIRS Hour Angle in degrees.
    dec_cirs : array_like
        The CIRS Declination in degrees.

    Returns
    -------
    bmx, bmy : ndarray
        The CHIME/FRB beam model X and Y coordinates in degrees as defined in
        the beam-model coordinate conventions:
        https://chime-frb-open-data.github.io/beam-model/#coordinate-conventions
    """
    from caput.interferometry import rotate_ypr, sph_to_ground

    from .observers import chime

    # Convert CIRS coordinates to CHIME "ground fixed" XYZ coordinates,
    # which constitute a unit vector pointing towards the point of interest,
    # i.e., telescope cartesian unit-sphere coordinates.
    # chx: The EW coordinate (increases to the East)
    # chy: The NS coordinate (increases to the North)
    # chz: The vertical coordinate (increases to the sky)
    chx, chy, chz = sph_to_ground(
        np.deg2rad(ha_cirs), np.deg2rad(chime.latitude), np.deg2rad(dec_cirs)
    )

    # Correct for CHIME telescope rotation with respect to North
    ypr = np.array([np.deg2rad(-chime.rotation), 0, 0])
    chx_rot, chy_rot, chz_rot = rotate_ypr(ypr, chx, chy, chz)

    # Convert rotated CHIME "ground fixed" XYZ coordinates to spherical polar
    # coordinates with the pole towards almost-North and using CHIME's meridian as the
    # prime meridian.
    # Note that the azimuthal angle theta in these spherical polar coordinates increases
    # to the West (to ensure that phi and theta here have the same meaning as the
    # variables with the same names in the beam_model package and DocLib #1203).
    # phi (polar angle): almost-North = 0 deg; zenith = 90 deg; almost-South = 180 deg
    # theta (azimuthal angle): almost-East = -90 deg; zenith = 0 deg;
    # almost-West = +90 deg
    phi = np.arccos(chy_rot)
    theta = np.arctan2(-chx_rot, +chz_rot)

    # Convert polar angle and azimuth to CHIME/FRB beam model XY position
    bmx = np.rad2deg(theta * np.sin(phi))
    bmy = np.rad2deg(np.pi / 2.0 - phi)

    return bmx, bmy


def bmxy_to_hadec(bmx, bmy):
    """Convert CHIME/FRB beam-model XY coordinates to CIRS hour angle and declination.

    Parameters
    ----------
    bmx, bmy : array_like
        The CHIME/FRB beam model X and Y coordinates in degrees as defined in
        the beam-model coordinate conventions:
        https://chime-frb-open-data.github.io/beam-model/#coordinate-conventions
        X is degrees west from the meridian
        Y is degrees north from zenith

    Returns
    -------
    ha_cirs : array_like
        The CIRS Hour Angle in degrees.
    dec_cirs : array_like
        The CIRS Declination in degrees.
    """
    from caput.interferometry import ground_to_sph, rotate_ypr

    from .observers import chime

    # Convert CHIME/FRB beam model XY position to spherical polar coordinates
    # with the pole towards almost-North and using CHIME's meridian as the prime
    # meridian. Note that the CHIME/FRB beam model X coordinate increases westward
    # and so does the azimuthal angle theta in these spherical polar coordinates
    # (to ensure that phi and theta here have the same meaning as the variables
    # with the same names in the beam_model package and DocLib #1203).
    # phi (polar angle): almost-North = 0 deg; zenith = 90 deg; almost-South = 180 deg
    # theta (azimuthal angle): almost-East = -90 deg; zenith = 0 deg;
    # almost-West = +90 deg
    phi = np.pi / 2.0 - np.deg2rad(bmy)
    theta = np.deg2rad(bmx) / np.sin(phi)

    # Warn for input beam-model XY positions below the horizon
    scalar_input = np.isscalar(theta)
    theta = np.atleast_1d(theta)
    if (theta < -1.0 * np.pi / 2.0).any() or (theta > np.pi / 2.0).any():
        warnings.warn("Input beam model XY coordinate(s) below horizon.")
    if scalar_input:
        theta = np.squeeze(theta)

    # Convert spherical polar coordinates to rotated CHIME "ground fixed" XYZ
    # coordinates (i.e., cartesian unit-sphere coordinates, rotated to correct
    # for the CHIME telescope's rotation with respect to North).
    # chx_rot: The almost-EW coordinate (increases to the almost-East)
    # chy_rot: The almost-NS coordinate (increases to the almost-North)
    # chz_rot: The vertical coordinate (increases to the sky)
    chx_rot = np.sin(phi) * np.sin(-theta)
    chy_rot = np.cos(phi)
    chz_rot = np.sin(phi) * np.cos(-theta)

    # Undo correction for CHIME telescope rotation with respect to North
    ypr = np.array([np.deg2rad(chime.rotation), 0, 0])
    chx, chy, chz = rotate_ypr(ypr, chx_rot, chy_rot, chz_rot)

    # Convert CHIME "ground fixed" XYZ coordinates to CIRS hour angle and declination
    ha_cirs, dec_cirs = ground_to_sph(chx, chy, np.deg2rad(chime.latitude))

    return np.rad2deg(ha_cirs), np.rad2deg(dec_cirs)


def get_range_rate(
    source: SkyfieldSource,
    date: float | list,
    obs: caput.time.Observer | None = None,
) -> float | np.array:
    """Calculate rate at which distance between observer and source changes.

    Parameters
    ----------
    source
        Position(s) on the sky.
    date
        Unix time(s) for which to calculate range rate.
    obs
        An Observer instance to use. If not supplied use `chime`. For many
        calculations changing from this default will make little difference.

    Returns
    -------
    range_rate
        Rate (in m/s) at which the distance between the observer and source
        changes (i.e., the velocity of observer in direction of source, but
        positive for observer and source moving appart). If either `source`
        or `date` contains multiple entries, `range_rate` will be an array.
        Otherwise, `range_rate` will be a float.

    Notes
    -----
    Only one of `source` and `date` can contain multiple entries.

    This routine uses an :class:`skyfield.positionlib.Apparent` object
    (rather than an :class:`skyfield.positionlib.Astrometric` object) to find
    the velocity of the observatory and the position of the source. This
    accounts for the gravitational deflection and the aberration of light.
    It is unclear if the latter should be taken into account for this Doppler
    shift calculation, but its effects are negligible.
    """
    if hasattr(source.ra._degrees, "__iter__") and hasattr(date, "__iter__"):
        raise ValueError(
            "Only one of `source` and `date` can contain multiple entries."
        )

    if obs is None:
        from .observers import chime as obs

    # Convert unix times to skyfield times
    date = unix_to_skyfield_time(date)

    # Create skyfield Apparent object of source position seen from observer
    position = obs.skyfield_obs().at(date).observe(source).apparent()

    # Observer velocity vector in ICRS xyz coordinates in units of m/s
    obs_vel_m_per_s = position.velocity.m_per_s

    # Normalized source position vector in ICRS xyz coordinates
    source_pos_m = position.position.m
    source_pos_norm = source_pos_m / np.linalg.norm(source_pos_m, axis=0)

    # Dot product of observer velocity and source position gives observer
    # velocity in direction of source; flip sign to get range rate (positive
    # for observer and source moving appart)
    return -np.sum(obs_vel_m_per_s.T * source_pos_norm.T, axis=-1)


def peak_ra(
    body: SkyfieldSource,
    date: float | None = None,
    deg: bool = False,
    obs: caput.time.Observer | None = None,
):
    """Calculate peak RA for a body.

    Calculates the right ascension where a source is expected to peak
    in the beam.  This is different than the transit RA for observers which
    are rotated with-respect-to north.

    Parameters
    ----------
    body : skyfield source
        skyfield.starlib.Star or skyfield.vectorlib.VectorSum or
        skyfield.jpllib.ChebyshevPosition representing the source.
    date : float
        Unix time at which to determine ra of source
        If None, 1 January, 2000 is used.
    deg : bool
        If true, return right ascension in degrees.  Otherwise,
        right ascension will be returned in radians.
    obs : caput.time.Observer or None
        The observer.  If not given, or None, the CHIME instrument will be used.

    Returns
    -------
    peak_ra : float
        RA when the transiting source peaks.
    """
    if obs is None:
        from .observers import chime as obs

    # Extract RA and dec of object
    ra, dec = object_coords(body, date=date)

    # Estimate the RA at which the transiting source peaks
    ra = ra + np.tan(obs.rotation) * (dec - obs.latitude) / np.cos(obs.latitude)

    # If requested, convert to degrees
    if deg:
        ra = np.degrees(ra)

    # Return
    return ra


if __name__ == "__main__":
    import doctest

    doctest.testmod()
