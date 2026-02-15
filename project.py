import re
from sympy import Matrix, lcm

SUBSCRIPT_MAP = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")

VALID_ELEMENTS = {
    "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
    "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar",
    "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni",
    "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr",
    "Rb", "Sr", "Ag", "Cd", "Sn", "Sb", "I", "Ba",
    "Pt", "Au", "Hg", "Pb"
}

def format_compound(compound):
    return compound.translate(SUBSCRIPT_MAP)

def check_parentheses(formula):
    stack = []
    for char in formula:
        if char == "(":
            stack.append(char)
        elif char == ")":
            if not stack:
                return False
            stack.pop()
    return len(stack) == 0

def main():
    equation = input("enter your unbalanced chemical equation!! ")

    try:
        balanced = balance_equation(equation)
        print("\nbalanced equation!")
        print(balanced)

    except Exception as e:
        print("\ninvalid equation :(")
        print(e)


def parse_equation(equation):
    if "->" not in equation:
        raise ValueError("equation must contain '->'")

    left, right = equation.split("->")

    reactants = [c.strip() for c in left.split("+")]
    products = [c.strip() for c in right.split("+")]

    if not reactants or not products:
        raise ValueError("equation must have reactants and products")

    if any(c == "" for c in reactants + products):
        raise ValueError("empty compound detected")

    return reactants, products

def count_atoms(formula):
    if not check_parentheses(formula):
        raise ValueError(f"missing parentheses in {formula}")

    stack = [{}]
    i = 0

    while i < len(formula):

        if formula[i] == "(":
            stack.append({})
            i += 1

        elif formula[i] == ")":
            i += 1
            num = ""
            while i < len(formula) and formula[i].isdigit():
                num += formula[i]
                i += 1

            multiplier = int(num) if num else 1
            top = stack.pop()

            for element in top:
                top[element] *= multiplier

            for element, count in top.items():
                stack[-1][element] = stack[-1].get(element, 0) + count

        else:
            match = re.match(r"^([A-Z][a-z]?)(\d*)", formula[i:])
            if not match:
                raise ValueError(f"invalid formula structure in {formula}")

            element, num = match.groups()

            if element not in VALID_ELEMENTS:
                raise ValueError(f"invalid element symbol : {element}")

            count = int(num) if num else 1
            stack[-1][element] = stack[-1].get(element, 0) + count

            i += len(match.group(0))

    return stack[0]

def balance_equation(equation):
    reactants, products = parse_equation(equation)
    compounds = reactants + products

    compound_atom_counts = [count_atoms(c) for c in compounds]

    reactant_elements = set()
    product_elements = set()

    for i, compound in enumerate(compounds):
        if i < len(reactants):
            reactant_elements.update(compound_atom_counts[i].keys())
        else:
            product_elements.update(compound_atom_counts[i].keys())

    if reactant_elements != product_elements:
        raise ValueError(f"element mismatch!")

    elements = list(reactant_elements)

    if not elements:
        raise ValueError("no valid elements detected")

    matrix = []

    for element in elements:
        row = []

        for i in range(len(reactants)):
            row.append(compound_atom_counts[i].get(element, 0))

        for i in range(len(reactants), len(compounds)):
            row.append(-compound_atom_counts[i].get(element, 0))

        matrix.append(row)

    m = Matrix(matrix)
    nullspace = m.nullspace()

    if not nullspace:
        raise ValueError("can't be balanced - no solution :(")

    coeffs = nullspace[0]
    lcm_val = lcm([term.q for term in coeffs])
    coeffs = coeffs * lcm_val
    coeffs = [int(c) for c in coeffs]

    if all(c <= 0 for c in coeffs):
        coeffs = [-c for c in coeffs]

    if any(c == 0 for c in coeffs):
        raise ValueError("can't be balanced - zero coefficient found :(")

    result = []

    for i, compound in enumerate(compounds):
        coefficient = coeffs[i]
        formatted = format_compound(compound)

        if coefficient == 1:
            result.append(formatted)

        else:
            result.append(f"{coefficient}{formatted}")

    return (
        " + ".join(result[:len(reactants)])
        + " -> "
        + " + ".join(result[len(reactants):])
    )

if __name__ == "__main__":
    main()