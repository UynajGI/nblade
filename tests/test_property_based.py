"""nblade Python 属性测试

使用 hypothesis 随机化验证几何代数的代数不变性。
"""

import math
from hypothesis import assume, given, settings, strategies as st

from nblade import Algebra

EPS = 1e-12  # relaxed tolerance for property tests (accumulated float error)


def assert_zero(mv, tol=EPS):
    """Assert multivector is effectively zero."""
    assert mv.norm() < tol, f"norm={mv.norm():.2e} > {tol}: {mv}"


def assert_equal(a, b, tol=EPS):
    """Assert two multivectors are effectively equal."""
    diff = a - b
    assert diff.norm() < tol, f"diff norm={diff.norm():.2e}: {diff}"


# ── Strategies ─────────────────────────────────────────────────────

dim_st = st.integers(min_value=2, max_value=5)
float_st = st.floats(-10.0, 10.0)
float_nz_st = st.floats(-10.0, 10.0).filter(lambda x: abs(x) > 0.1)
angle_st = st.floats(0.2, math.pi - 0.2)


@st.composite
def alg_and_vec(draw):
    dim = draw(dim_st)
    alg = Algebra.euclidean(dim)
    coeffs = draw(st.lists(float_st, min_size=dim, max_size=dim))
    v = alg.vector(coeffs)
    assume(not v.is_zero())
    return alg, v


@st.composite
def alg_and_two_vecs(draw):
    dim = draw(dim_st)
    alg = Algebra.euclidean(dim)
    a = draw(st.lists(float_st, min_size=dim, max_size=dim))
    b = draw(st.lists(float_st, min_size=dim, max_size=dim))
    va = alg.vector(a)
    vb = alg.vector(b)
    assume(not va.is_zero() and not vb.is_zero())
    return alg, va, vb


@st.composite
def alg_and_two_indices(draw):
    dim = draw(dim_st)
    alg = Algebra.euclidean(dim)
    i = draw(st.integers(0, dim - 1))
    j = draw(st.integers(0, dim - 1).filter(lambda x: x != i))
    return alg, i, j


# ── Geometric Product ──────────────────────────────────────────────

class TestGeometricProduct:

    @given(alg_and_two_vecs(), st.lists(float_st, min_size=1, max_size=5))
    @settings(max_examples=200)
    def test_associativity(self, data, c_coeffs):
        alg, a, b = data
        # pad c to match dimension
        coeffs = (c_coeffs + [0.0] * max(0, alg.dimension - len(c_coeffs)))[:alg.dimension]
        c = alg.vector(coeffs)
        assume(not c.is_zero())
        assert_zero((a * b) * c - a * (b * c))

    @given(alg_and_two_vecs(), st.lists(float_st, min_size=1, max_size=5))
    @settings(max_examples=200)
    def test_distributivity(self, data, c_coeffs):
        alg, a, b = data
        coeffs = (c_coeffs + [0.0] * max(0, alg.dimension - len(c_coeffs)))[:alg.dimension]
        c = alg.vector(coeffs)
        assume(not c.is_zero())
        assert_zero(a * (b + c) - (a * b + a * c))


# ── Outer Product ──────────────────────────────────────────────────

class TestOuterProduct:

    @given(alg_and_two_vecs())
    @settings(max_examples=200)
    def test_antisymmetry(self, data):
        _, a, b = data
        assert_zero((a ^ b) + (b ^ a))

    @given(alg_and_vec())
    @settings(max_examples=200)
    def test_self_wedge_zero(self, data):
        _, v = data
        assert (v ^ v).is_zero()  # exact — should be identically zero

    @given(alg_and_two_vecs(), st.lists(float_st, min_size=1, max_size=5))
    @settings(max_examples=200)
    def test_associativity(self, data, c_coeffs):
        alg, a, b = data
        coeffs = (c_coeffs + [0.0] * max(0, alg.dimension - len(c_coeffs)))[:alg.dimension]
        c = alg.vector(coeffs)
        assume(not c.is_zero())
        # (a^b)^c - a^(b^c) should be 0, but float error accumulates
        diff = ((a ^ b) ^ c) - (a ^ (b ^ c))
        assert diff.norm() < 1e-10, f"dim={alg.dimension} norm={diff.norm():.2e}"


# ── Inner Product ──────────────────────────────────────────────────

class TestInnerProduct:

    @given(alg_and_two_indices())
    @settings(max_examples=200)
    def test_orthogonality(self, data):
        alg, i, j = data
        ei = alg.basis_vector(i)
        ej = alg.basis_vector(j)
        assert (ei | ej).is_zero()  # orthogonal — should be exactly zero

    @given(alg_and_two_indices())
    @settings(max_examples=200)
    def test_same_basis_square(self, data):
        alg, i, _ = data
        ei = alg.basis_vector(i)
        inner = ei | ei
        assert abs(inner.scalar_part() - 1.0) < 1e-10


# ── Involutions ────────────────────────────────────────────────────

class TestInvolutions:

    @given(alg_and_vec())
    @settings(max_examples=200)
    def test_grade_involution_idempotent(self, data):
        _, v = data
        assert_equal(v.grade_involution().grade_involution(), v)

    @given(alg_and_vec())
    @settings(max_examples=200)
    def test_reversion_idempotent(self, data):
        _, v = data
        assert_equal(v.reversion().reversion(), v)

    @given(alg_and_vec())
    @settings(max_examples=200)
    def test_clifford_conjugate_idempotent(self, data):
        _, v = data
        assert_equal(v.clifford_conjugate().clifford_conjugate(), v)

    @given(alg_and_vec())
    @settings(max_examples=200)
    def test_norm_invariance(self, data):
        _, v = data
        n = v.norm()
        assert abs(v.grade_involution().norm() - n) < 1e-10
        assert abs(v.reversion().norm() - n) < 1e-10
        assert abs(v.clifford_conjugate().norm() - n) < 1e-10


# ── Dual ───────────────────────────────────────────────────────────

class TestDual:

    @given(alg_and_vec())
    @settings(max_examples=200)
    def test_roundtrip_vector(self, data):
        """dual(dual(v)) = v (odd dim) or -v (even dim) for Euclidean."""
        alg, v = data
        roundtrip = v.dual().dual()
        # odd dim -> +v, even dim -> -v
        expected = v if alg.dimension % 2 == 1 else v.scale(-1.0)
        assert_equal(roundtrip, expected)


# ── Inverse ────────────────────────────────────────────────────────

class TestInverse:

    @given(alg_and_vec())
    @settings(max_examples=200)
    def test_vector_inverse(self, data):
        _, v = data
        assume(v.is_invertible())
        inv = v.inverse()
        product = v * inv
        assert abs(product.scalar_part() - 1.0) < 1e-10

    @given(alg_and_two_indices())
    @settings(max_examples=200)
    def test_bivector_inverse(self, data):
        alg, i, j = data
        B = alg.basis_vector(i) ^ alg.basis_vector(j)
        assume(not B.is_zero())
        inv = B.inverse()
        product = B * inv
        assert abs(product.scalar_part() - 1.0) < 1e-10

    @given(alg_and_vec())
    @settings(max_examples=200)
    def test_inverse_of_inverse(self, data):
        _, v = data
        assume(v.is_invertible())
        inv = v.inverse()
        assert_equal(inv.inverse(), v)


# ── Rotation ───────────────────────────────────────────────────────

class TestRotation:

    @given(alg_and_two_indices(), angle_st, st.lists(float_st, min_size=1, max_size=5))
    @settings(max_examples=100)
    def test_preserves_norm(self, data, angle, coeffs):
        alg, i, j = data
        plane = alg.basis_vector(i) ^ alg.basis_vector(j)
        rotor = alg.rotor(plane, angle)
        dim = alg.dimension
        c = (coeffs + [0.0] * max(0, dim - len(coeffs)))[:dim]
        v = alg.vector(c)
        assume(not v.is_zero())
        rotated = v.rotate_by(rotor)
        assert abs(rotated.norm() - v.norm()) < 1e-10

    @given(alg_and_two_indices(), angle_st, st.lists(float_st, min_size=1, max_size=5))
    @settings(max_examples=100)
    def test_double_rotation(self, data, angle, coeffs):
        alg, i, j = data
        plane = alg.basis_vector(i) ^ alg.basis_vector(j)
        rotor = alg.rotor(plane, angle)
        rotor2 = alg.rotor(plane, angle * 2)
        dim = alg.dimension
        c = (coeffs + [0.0] * max(0, dim - len(coeffs)))[:dim]
        v = alg.vector(c)
        assume(not v.is_zero())
        once = v.rotate_by(rotor)
        twice = once.rotate_by(rotor)
        expected = v.rotate_by(rotor2)
        assert_equal(twice, expected)


# ── Reflection ─────────────────────────────────────────────────────

class TestReflection:

    @given(alg_and_two_indices(), st.lists(float_st, min_size=1, max_size=5))
    @settings(max_examples=200)
    def test_double_reflection_identity(self, data, coeffs):
        alg, i, _ = data
        n = alg.basis_vector(i)
        dim = alg.dimension
        c = (coeffs + [0.0] * max(0, dim - len(coeffs)))[:dim]
        v = alg.vector(c)
        assume(not v.is_zero())
        assert_equal(v.reflect_in(n).reflect_in(n), v)


# ── Commutator ─────────────────────────────────────────────────────

class TestCommutator:

    @given(alg_and_two_vecs())
    @settings(max_examples=200)
    def test_anticommutativity(self, data):
        _, a, b = data
        assert_zero(a.commutator(b) + b.commutator(a))

    @given(alg_and_vec())
    @settings(max_examples=200)
    def test_self_zero(self, data):
        _, v = data
        assert v.commutator(v).is_zero()  # exact — should be identically zero


# ── Norm ───────────────────────────────────────────────────────────

class TestNorm:

    @given(alg_and_vec())
    @settings(max_examples=200)
    def test_scalar_product_self_equals_norm_sq(self, data):
        _, v = data
        assert abs(v.scalar_product(v) - v.norm_squared()) < 1e-10

    @given(alg_and_two_vecs())
    @settings(max_examples=200)
    def test_scalar_product_symmetric(self, data):
        _, a, b = data
        assert abs(a.scalar_product(b) - b.scalar_product(a)) < 1e-10
