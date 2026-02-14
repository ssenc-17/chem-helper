import re
from sympy import Matrix, lcm

def main():
    equation = input("enter your unbalanced chemical equation!! ")
    
    try:
        balanced = balance_equation(equation)
        print("balanced equation!")
        print(balanced)

    except Exception:
        print("not valid :(")

def parse_equation(equation):
    if "->" not in equation:
        raise ValueError("equation must contain a '->' btww")
    
    left, right = equation.split("->")
    reactants = [c.strip() for c in left.split("+")]
    products = [c.strip() for c in right.split("+")]

    return reactants, products

def count_atoms(compound):
    pattern = r"([A-Z][a-z]?)(\d*)"
    matches = re.findall(pattern, compound)

    atoms = {}
    for element, count in matches:
        count = int(count) if count else 1
        atoms[element] = atoms.get(element, 0) + count

    return atoms

def balance_equation(equation):
    reactants, products = parse_equation(equation)
    compounds = reactants + products

    elements = set()

    for compound in compounds:
        elements.update(count_atoms(compound).keys())
    elements = list(elements)

    matrix = []

    for element in elements:
        row = []
        for compound in reactants:
            row.append(count_atoms(compound).get(element, 0))
        for compound in products:
            row.append(-count_atoms(compound).get(element, 0))
        matrix.append(row)

    m = Matrix(matrix)
    nullspace = m.nullspace()

    if not nullspace:
        raise ValueError("whoops - no solution found, sorry!!")

    coeffs = nullspace[0]
    lcm_val = lcm([term.q for term in coeffs])
    coeffs = coeffs * lcm_val
    coeffs = [abs(int(c)) for c in coeffs]

    result = []

    for i, compound in enumerate(compounds):
        coefficient = coeffs[i]
        
        if coefficient == 1:
            result.append(compound)
        else:
            result.append(f"{coefficient}{compound}")

    return " + ".join(result[:len(reactants)]) + " -> " + " + ".join(result[len(reactants):])

if __name__ == "__main__":
    main()