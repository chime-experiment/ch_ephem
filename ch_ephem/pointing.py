"""CHIME ephemeris pointing routines

Functions
=========
- `py:meth:galt_pointing_model_ha`
- `py:meth:galt_pointing_model_dec`
"""

from __future__ import annotations
from typing import List

import numpy as np
from skyfield.units import Angle


def galt_pointing_model_ha(
    ha_in: Angle,
    dec_in: Angle,
    a: List[float] = [-5.872, -0.5292, 5.458, -0.076, -0.707, 0.0, 0.0],
) -> Angle:
    """Calculate pointing correction in hour angle for the Galt Telescope
    See description of the pointing model by Lewis Knee: CHIME doclib #754

    Parameters
    ----------
    ha, dec : Skyfield Angle objects
        Target hour angle and declination

    a : list of floats
        List of coefficients (in arcmin) for the pointing model
        (NOTE: it is very unlikely that a user will want to change these
        from the defaults, which are taken from the pointing model as of
        2019-2-15)

    Returns
    -------
    correction : Angle
        Angular offset in hour angle
    """

    ha = ha_in.radians
    dec = dec_in.radians

    # hour angle pointing correction in arcmin
    delta_ha_cos_dec = (
        a[0]
        + a[1] * np.sin(dec)
        + a[2] * np.cos(dec)
        + a[3] * np.sin(ha) * np.sin(dec)
        + a[4] * np.cos(ha) * np.sin(dec)
        + a[5] * np.sin(ha) * np.cos(dec)
        + a[6] * np.cos(ha) * np.cos(dec)
    )

    return Angle(degrees=(delta_ha_cos_dec / np.cos(dec)) / 60.0)


def galt_pointing_model_dec(
    ha_in: Angle,
    dec_in: Angle,
    b: List[float] = [1.081, 0.707, -0.076, 0.0, 0.0, 0.0, 0.0],
) -> Angle:
    """Calculate pointing correction in declination for the Galt Telescope
    See description of the pointing model by Lewis Knee: CHIME doclib #754

    Parameters
    ----------
    ha, dec : Skyfield Angle objects
        Target hour angle and declination

    b : list of floats
        List of coefficients (in arcmin) for the pointing model
        (NOTE: it is very unlikely that a user will want to change these
        from the defaults, which are taken from the pointing model as of
        2019-2-15)

    Returns
    -------
    correction : Angle
        Angular offset in hour angle
    """

    ha = ha_in.radians
    dec = dec_in.radians

    # declination pointing correction in arcmin
    delta_dec = (
        b[0]
        + b[1] * np.sin(ha)
        + b[2] * np.cos(ha)
        + b[3] * np.sin(dec)
        + b[4] * np.cos(dec)
        + b[5] * np.sin(dec) * np.cos(ha)
        + b[6] * np.sin(dec) * np.sin(ha)
    )

    return Angle(degrees=delta_dec / 60.0)
