import pytest
from project import format_compound, check_parentheses, parse_equation, count_atoms, calculate_molar_mass, balance_equation

def test_format_compound():
    assert format_compound("H2O") == "H₂O"
    assert format_compound("C6H12O6") == "C₆H₁₂O₆"
    assert format_compound("NaCl") == "NaCl"

def test_check_parentheses():
    assert check_parentheses("(H2O)") is True
    assert check_parentheses("((NH4)2SO4)") is True
    assert check_parentheses("(H2O") is False
    assert check_parentheses("H2O)") is False
    assert check_parentheses("H2O") is True

def test_parse_equation():
    reactants, products = parse_equation("H2 + O2", "H2O")
    assert reactants == ["H2", "O2"]
    assert products == ["H2O"]

    with pytest.raises(ValueError):
        parse_equation("", "H2O")

    with pytest.raises(ValueError):
        parse_equation("H2 + O2", "")

def test_count_atoms_simple():
    assert count_atoms("H2O") == {"H": 2, "O": 1}
    assert count_atoms("C6H12O6") == {"C": 6, "H": 12, "O": 6}
    assert count_atoms("NaCl") == {"Na": 1, "Cl": 1}

def test_count_atoms_with_parentheses():
    assert count_atoms("(NH4)2SO4") == {"N": 2, "H": 8, "S": 1, "O": 4}
    assert count_atoms("Al2(SO4)3") == {"Al": 2, "S": 3, "O": 12}

def test_count_atoms_invalid():
    with pytest.raises(ValueError):
        count_atoms("H2O)")
    with pytest.raises(ValueError):
        count_atoms("Xy2")

def test_calculate_molar_mass():
    assert abs(calculate_molar_mass("H2O") - 18.016) < 0.001
    assert abs(calculate_molar_mass("CO2") - 44.01) < 0.01

def test_balance_equation_simple():
    balanced, coeffs, compounds = balance_equation("H2 + O2", "H2O")
    assert balanced == "2H₂ + O₂ -> 2H₂O"
    assert coeffs == [2, 1, 2]
    assert compounds == ["H2", "O2", "H2O"]

def test_balance_equation_complex():
    balanced, coeffs, compounds = balance_equation("C3H8 + O2", "CO2 + H2O")
    assert balanced == "C₃H₈ + 5O₂ -> 3CO₂ + 4H₂O"
    assert coeffs == [1, 5, 3, 4]
    assert compounds == ["C3H8", "O2", "CO2", "H2O"]

def test_balance_equation_unbalanceable():
    with pytest.raises(ValueError):
        balance_equation("H2 + O2", "CO2")