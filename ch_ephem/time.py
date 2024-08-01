"""CHIME ephemeris time routines.

Time Functions
==============

- :py:meth:`parse_date`
- :py:meth:`utc_lst_to_mjd`
- :py:meth:`chime_local_datetime`
"""

from __future__ import annotations
from typing import Optional, Any

import re
import datetime
from pytz import timezone

from caput.time import datetime_to_unix, unix_to_skyfield_time, Observer


def parse_date(datestring: str) -> datetime.datetime:
    """Convert date string to a datetime object.

    Parameters
    ----------
    datestring : str
        Date as YYYYMMDD-AAA, where AAA is one of [UTC, EST, EDT, PST, PDT]

    Returns
    -------
    date : datetime.datetime
        A python datetime object in UTC.
    """
    rm = re.match("([0-9]{8})-([A-Z]{3})", datestring)
    if rm is None:
        msg = (
            "Wrong format for datestring: {0}.".format(datestring)
            + "\nShould be YYYYMMDD-AAA, "
            + "where AAA is one of [UTC,EST,EDT,PST,PDT]"
        )
        raise ValueError(msg)

    datestring = rm.group(1)
    tzoffset = 0.0
    tz = rm.group(2)

    tzs = {"PDT": -7.0, "PST": -8.0, "EDT": -4.0, "EST": -5.0, "UTC": 0.0}

    if tz is not None:
        try:
            tzoffset = tzs[tz.upper()]
        except KeyError:
            print("Time zone {} not known. Known time zones:".format(tz))
            for key, value in tzs.items():
                print(key, value)
            print("Using UTC{:+.1f}.".format(tzoffset))

    return datetime.datetime.strptime(datestring, "%Y%m%d") - datetime.timedelta(
        hours=tzoffset
    )


def utc_lst_to_mjd(
    datestring: str, lst: float, obs: Optional[Observer] = None
) -> float:
    """Convert date string and LST to corresponding modified Julian Day

    Parameters
    ----------
    datestring : str
        Date.  The string is parsed using `parse_date` (q.v.).
    lst : float
        Local sidereal time at the observer in decimal hours
    obs : caput.time.Observer or None
        The observer.  If not given, or None, the CHIME
        instrument will be used.

    Returns
    -------
    mjd : float
        Modified Julian Date corresponding to the given time.
    """
    if obs is None:
        from .observers import chime as obs

    return (
        unix_to_skyfield_time(
            obs.lsa_to_unix(lst * 360 / 24, datetime_to_unix(parse_date(datestring)))
        ).tt
        - 2400000.5
    )


def chime_local_datetime(*args: Any) -> datetime.datetime:
    """Create a :class:`datetime.datetime` object in Canada/Pacific timezone.

    Parameters
    ----------
    *args
        Any valid arguments to the constructor of :class:`datetime.datetime`
        except *tzinfo*. Local date and time at CHIME.

    Returns
    -------
    dt : :class:`datetime.datetime`
        Timezone naive date and time but converted to UTC.

    """

    tz = timezone("Canada/Pacific")
    dt_naive = datetime.datetime(*args)
    if dt_naive.tzinfo:
        raise ValueError("Time zone should not be supplied.")
    dt_aware = tz.localize(dt_naive)
    return dt_aware.replace(tzinfo=None) - dt_aware.utcoffset()
