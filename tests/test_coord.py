"""Test ch_ephem.coord."""

import pytest
import ch_ephem.coord

def test_hadec_to_bmxy():
   """Test coord.hadec_to_bmxy."""

  result = hadec_to_bmxy(0, 0)
  assert result[0] == pytest.approx(-0.053844273720663124)
  assert result[1] == pytest.approx(-49.32065803774193)


def test_bmxy_to_hadec():
   """Test coord.bmxy_to_hadec."""

  result = bmxy_to_hadec(0, 0)
  assert result[0] == pytest.approx(-6.669734114397701e-18)
  assert result[1] == pytest.approx(49.3207092194)
