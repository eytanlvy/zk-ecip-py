from src.divisor import Divisor
from src.polynomial import Polynomial, test_colinearity
from src.field import BaseFieldElement as Felt, BaseField
from src.curve import G1Point, POINT_AT_INFINITY, Fp


def generate_principal_divisor_sequence(len):
    base_point = G1Point.gen_random_point()
    points = [base_point]
    for _ in range(1, len):
        next_point = G1Point.gen_random_point()
        points.append(next_point)
        base_point = base_point - next_point
    return points

def two_points(points):
    assert (len(points) == 2);
    q = -points[0] - points[1]
    div = Divisor({points[0]: 1, points[1]: 1, q: 1, POINT_AT_INFINITY: -3})
    return div

def recursive_sum(points, start_index=0):
    n = len(points)
    if n - start_index == 1:
        # Base case: single point
        return points[start_index]
    elif n - start_index == 2:
        # Base case: two points (Q_i)
        return -points[start_index] - points[start_index + 1]
    else:
        # Recursive case
        mid_index = start_index + (n - start_index) // 2
        sum_left = recursive_sum(points[:mid_index], start_index)
        sum_right = recursive_sum(points[mid_index:], mid_index)
        sum_q = -sum_left - sum_right

        # Calculate Q'_i
        q_prime = -(-points[start_index] - points[start_index + 1] - points[mid_index] - points[mid_index + 1])

        # Return the recursive sum P_{4i} + P_{4i+1} + P_{4i+2} + P_{4i+3} + Q'_i - 5O
        return Divisor({sum_left: 1, sum_right: 1, q_prime: 1, POINT_AT_INFINITY: -5})
    
def incremental_construction(points):
    assert (len(points) == 1, "One point can't form a principal divisor");
    if (len(points) == 2):
        assert(points[0] == -points[1], "Two points must be opposite");
        return two_points(points)
    elif (len(points) == 3):
        assert(test_colinearity(points), "Three points must be collinear");
        return two_points([points[0], points[1]])
    else:
        return recursive_sum(points)
    
if __name__ == "__main__":
    print(3//2)
    p = G1Point.gen_random_point()
    q = G1Point.gen_random_point()
    r = G1Point(Fp(1), Fp(2))
    D = two_points([p, q])
    print(f"D: {D.points} \ndegree: {D.degree}\n")
    domain = [p.x, q.x]
    values = [p.y, q.y]
    polynomial = Polynomial.interpolate_domain(domain, values)
    print(f"Équation de la droite entre p et q : {polynomial}")

    print("Génération de la séquence de diviseurs principaux...")
    sequence = generate_principal_divisor_sequence(6)
    sum_of_points = G1Point.zero()
    for point in sequence:
        sum_of_points = sum_of_points + point
    print(f"La somme des points est égale à zéro : {sum_of_points == G1Point.zero()}")

    print("Algo 3.1.1...")
    sequence = [G1Point.gen_random_point() for _ in range(6)]
    div = incremental_construction(sequence)
    print(f"div: {div.points} \ndegree: {div.degree}\n")

