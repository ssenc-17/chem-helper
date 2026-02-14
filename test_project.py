import pytest
from project import parse_equation, count_atoms, balance_equation

def test_parse_equation():
    r, p = parse_equation("H2 + O2 -> H2O")
    assert r == ["H2", "O2"]
    assert p == ["H2O"]

def test_parse_equation_invalid():
    with pytest.raises(ValueError):
        parse_equation("H2 + O2 H2O")

def test_count_atoms():
    assert count_atoms("H2O") == {"H": 2, "O": 1}
    assert count_atoms("CO2") == {"C": 1, "O": 2}

def test_balance_equation():
    assert balance_equation("H2 + O2 -> H2O") == "2H2 + O2 -> 2H2O"