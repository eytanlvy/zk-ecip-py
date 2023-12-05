from src.divisor import Divisor
from src.polynomial import Polynomial
from src.field import BaseFieldElement as Felt, BaseField
from src.curve import G1Point, POINT_AT_INFINITY, Fp

def two_points(points):
    assert (len(points) == 2);
    q = -points[0] - points[1]
    div = Divisor({points[0]: 1, points[1]: 1, q: 1, POINT_AT_INFINITY: -3})
    return div


p = G1Point.gen_random_point()
q = G1Point.gen_random_point()
D = two_points([p, q])
print ("HEEEEEERE")
print(f"D: {D.points} \ndegree: {D.degree}\n")
domain = [p.x, q.x]
values = [p.y, q.y]

polynomial = Polynomial.interpolate_domain(domain, values)
print(f"Équation de la droite entre p et q : {polynomial}")

