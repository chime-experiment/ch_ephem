"""Test ch_ephem.observers."""

import pickle

from ch_ephem import observers


def test_pickle_observers():
    """Observers should be picklable."""
    # These are all small powers of two to
    # avoid floating-point shennanigans
    test = observers.Observer(lat=1, lon=2, alt=4, rot=8, roll=16, offset=(32, 64, 128))

    pickled = pickle.dumps(test)
    assert pickled is not None

    result = pickle.loads(pickled)

    assert result.latitude == 1
    assert result.longitude == 2
    assert result.altitude == 4
    assert result.rotation == 8
    assert result.roll == 16
    assert result.offset.x == 32
    assert result.offset.y == 64
    assert result.offset.z == 128

    repickled = pickle.dumps(result)

    assert pickled == repickled
