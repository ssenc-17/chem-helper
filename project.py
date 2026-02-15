import re
from sympy import Matrix, lcm

SUBSCRIPT_MAP = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")

def format_compound(compound):
    return compound.translate(SUBSCRIPT_MAP)

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

def count_atoms(formula):
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
            match = re.match(r"([A-Z][a-z]?)(\d*)", formula[i:])
            if not match:
                raise ValueError("invalid format :(")

            element, num = match.groups()
            count = int(num) if num else 1

            stack[-1][element] = stack[-1].get(element, 0) + count

            i += len(match.group(0))

    return stack[0]

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
        formatted = format_compound(compound)

        if coefficient == 1:
            result.append(formatted)
        
        else:
            result.append(f"{coefficient}{formatted}")

    return " + ".join(result[:len(reactants)]) + " -> " + " + ".join(result[len(reactants):])

if __name__ == "__main__":
    main()