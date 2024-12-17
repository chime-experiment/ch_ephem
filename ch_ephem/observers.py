"""Observers for CHIME instruments

This module provides `caput.time.Observer` objects for CHIME instruments.
Position data for the instruments comes from the `instruments.yaml` file
provided with this module.

Any instrument defined in `instruments.yaml` will automatically have an
importable object created for it in `ch_ephem.observers`.  To retrieve a
particular Observer, import it from `observers`:

    >>> from ch_ephem.observers import chime

To get a list of all available instrument observers, use the :py:meth:`all`
function, which returns a dict of Observer objects:

    >>> import ch_ephem.observers
    >>> ch_ephem.observers.all()
    {'chime': <ch_ephem.observers.Observer object at 0x7f2225b264d0>,
    'galt': <ch_ephem.observers.Observer object at 0x7f2225e76110>,
    ...
    }

Observers objects found in this module subclass the standard
`caput.time.Observer` to additionally provide rotation, roll and tangent-space
offset for the instruments:

    >>> from ch_ephem.observers import pathfinder
    >>> print(pathfinder.rotation)
    1.986
    >>> print(pathfinder.roll)
    0.0
    >>> print(pathfinder.offset.x)
    373.754961
    >>> print(list(pathfinder.offset))
    [373.754961, -54.649866, 0.0]

Functions
=========

- :py:meth:`all`
- :py:meth:`reset`
"""

from __future__ import annotations

import datetime
import warnings
import pathlib
import yaml
from collections import namedtuple

from caput.time import Observer as CaputObserver

# LSD start time.  This is the same for all observers
_lsd_start = datetime.datetime(2013, 11, 15)

# Initialised on first use
_observers = None

# tangent-space offset representation.
Offset = namedtuple("Offset", ["x", "y", "z"])


class Observer(CaputObserver):
    """Representation of a local observer.

    This class extends `caput.time.Observer` to add CHIME-specific geometry
    parameters:

    Observer.rotation:
        rotation of the cylinder(s) in degrees anti-clockwise looking down from
        above.  (e.g. West from North)
    Observer.roll:
        roll of the cylinder focal-line in degree eastward from vertical
        (e.g. clockwise looking North along the focal line).
    Observer.offset:
        offset in metres of the antenna array zero point in the right-handed
        tangent-space anchored at (lat,lon,alt).  A positive X offset is
        eastward.  A positive Y offset is northward.  A positive Z offset
        is upward.
    """

    def __init__(
        self,
        lon=0.0,
        lat=0.0,
        alt=0.0,
        rot=0.0,
        roll=0.0,
        offset=[0.0, 0.0, 0.0],
        lsd_start=None,
        sf_wrapper=None,
    ):
        super().__init__(lon, lat, alt, lsd_start, sf_wrapper)
        self.rotation = rot
        self.roll = roll
        self.offset = Offset(*offset)


def all() -> dict[str, Observer]:
    """Returns a dict of all available Observers, keyed by name."""
    global _observers
    if _observers is None:
        with pathlib.Path(__file__).with_name("instruments.yaml").open() as f:
            data = yaml.safe_load(f)

        _observers = dict()

        for inst in data:
            # Vet YAML record
            missing = list()
            for key in "latitude", "longitude", "altitude", "rotation", "roll":
                if key not in data[inst]:
                    missing.append(key)

            if "offset" in data[inst]:
                for key in "x", "y", "z":
                    if key not in data[inst]["offset"]:
                        missing.append("offset." + key)
            else:
                missing.append("offset")

            if missing:
                warnings.warn(
                    f'Unable to create observer for "{inst}": missing {missing}'
                )
                continue

            _observers[inst] = Observer(
                lat=float(data[inst]["latitude"]),
                lon=float(data[inst]["longitude"]),
                alt=float(data[inst]["altitude"]),
                rot=float(data[inst]["rotation"]),
                roll=float(data[inst]["roll"]),
                offset=(
                    float(data[inst]["offset"]["x"]),
                    float(data[inst]["offset"]["y"]),
                    float(data[inst]["offset"]["z"]),
                ),
                lsd_start=_lsd_start,
            )

    return _observers


def reset() -> None:
    """Unload all the observers.

    Use this function if you've changed `instruments.yaml` on disk
    to clear the module's cache of Observers.
    """
    global _observers
    _observers = None


def __getattr__(name: str) -> Observer:
    """Retreive the `Observer` for the instrument named.

    Don't call this function directly; just import the
    instrument you want from the module:

        >>> from ch_ephem.observers import pathfinder

    Parameters
    ----------
    name:
        name of the instrument Observer to retrieve

    Returns
    -------
    observer:
        The Observer object for the named instrument

    Raises
    ------
    AttributeError:
        Data for the instrument named could not be found.
    """

    # Load observers from disk, if necessary
    observers = all()

    try:
        return observers[name]
    except KeyError as e:
        raise AttributeError(f"Unknown instrument: {name}") from e
