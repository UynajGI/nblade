"""
nblade Python 绑定测试

测试 Python API 的完整性和正确性。
"""

import pytest
import math
from nblade import Algebra, MultiVector, AlgebraConfig, create_rotor, __version__


class TestVersion:
    """版本测试"""

    def test_version_exists(self):
        """测试版本号存在"""
        assert __version__ is not None
        assert isinstance(__version__, str)
        assert len(__version__) > 0


class TestAlgebra:
    """Algebra 类测试"""

    def test_create_euclidean(self):
        """创建欧几里得代数"""
        algebra = Algebra.euclidean(3)
        assert algebra.dimension == 3
        assert algebra.signature == (3, 0, 0)
        assert algebra.basis_count == 3

    def test_create_spacetime(self):
        """创建时空代数"""
        algebra = Algebra.spacetime(4)
        assert algebra.dimension == 4
        assert algebra.signature == (1, 3, 0)

    def test_create_cga(self):
        """创建共形几何代数"""
        algebra = Algebra.cga()
        assert algebra.dimension == 5
        assert algebra.signature == (4, 1, 0)

    def test_create_custom(self):
        """创建自定义签名代数"""
        algebra = Algebra(4, p=2, q=1, r=1)
        assert algebra.dimension == 4
        assert algebra.signature == (2, 1, 1)

    def test_basis_vectors(self):
        """测试基向量"""
        algebra = Algebra.euclidean(3)
        e1, e2, e3 = algebra.basis_vectors()

        assert e1.config.dimension == 3
        assert e2.config.dimension == 3
        assert e3.config.dimension == 3

    def test_scalar(self):
        """测试标量创建"""
        algebra = Algebra.euclidean(3)
        scalar = algebra.scalar(5.0)
        assert scalar.scalar_part() == 5.0

    def test_one_and_zeros(self):
        """测试单位和零多向量"""
        algebra = Algebra.euclidean(3)

        one = algebra.one()
        assert one.scalar_part() == 1.0

        zero = algebra.zeros()
        assert zero.is_zero()

    def test_repr_and_str(self):
        """测试字符串表示"""
        algebra = Algebra.euclidean(3)

        repr_str = repr(algebra)
        assert "Algebra" in repr_str
        assert "3" in repr_str

        str_str = str(algebra)
        assert "Geometric Algebra" in str_str


class TestMultiVector:
    """MultiVector 类测试"""

    @pytest.fixture
    def algebra_3d(self):
        """3D 欧几里得代数"""
        return Algebra.euclidean(3)

    def test_basis_vector_square(self, algebra_3d):
        """测试基向量平方"""
        e1 = algebra_3d.basis_vector(0)
        e2 = algebra_3d.basis_vector(1)

        e1_sq = e1 * e1
        e2_sq = e2 * e2

        assert abs(e1_sq.scalar_part() - 1.0) < 1e-10
        assert abs(e2_sq.scalar_part() - 1.0) < 1e-10

    def test_geometric_product(self, algebra_3d):
        """测试几何积"""
        e1 = algebra_3d.basis_vector(0)
        e2 = algebra_3d.basis_vector(1)

        # e1 * e2 = e12
        product = e1 * e2

        # e2 * e1 = -e12
        product_rev = e2 * e1

        # 应该是反对称的
        sum_prod = product + product_rev
        assert sum_prod.is_zero()

    def test_outer_product(self, algebra_3d):
        """测试外积"""
        e1 = algebra_3d.basis_vector(0)
        e2 = algebra_3d.basis_vector(1)

        # e1 ^ e2 = e12
        wedge = e1 ^ e2

        # e1 ^ e1 = 0
        wedge_self = e1 ^ e1
        assert wedge_self.is_zero()

    def test_inner_product(self, algebra_3d):
        """测试内积"""
        e1 = algebra_3d.basis_vector(0)
        e2 = algebra_3d.basis_vector(1)

        # 正交向量内积为零
        inner = e1 | e2
        assert inner.is_zero()

        # 相同向量内积为平方
        inner_self = e1 | e1
        assert abs(inner_self.scalar_part() - 1.0) < 1e-10

    def test_grade_involution(self, algebra_3d):
        """测试阶次对合"""
        e1 = algebra_3d.basis_vector(0)
        e2 = algebra_3d.basis_vector(1)

        # 向量（奇阶次）变号
        star_e1 = ~e1
        expected = -e1

        # 检查是否相等（需要比较系数）
        diff = star_e1 - expected
        assert diff.is_zero()

    def test_reversion(self, algebra_3d):
        """测试反转"""
        e1 = algebra_3d.basis_vector(0)
        e2 = algebra_3d.basis_vector(1)

        plane = e1 * e2
        rev = plane.reversion()

        # (e1e2)† = e2e1 = -e1e2
        expected = -plane
        diff = rev - expected
        assert diff.is_zero()

    def test_clifford_conjugate(self, algebra_3d):
        """测试 Clifford 共轭"""
        e1 = algebra_3d.basis_vector(0)

        conj = e1.clifford_conjugate()
        expected = -e1

        diff = conj - expected
        assert diff.is_zero()

    def test_dual(self, algebra_3d):
        """测试对偶"""
        e1, e2, e3 = algebra_3d.basis_vectors()

        # e1 的对偶应该是 e23
        dual_e1 = e1.dual()

        # 检查阶次：向量的对偶是二重向量
        grade_2_part = dual_e1.grade(2)
        assert not grade_2_part.is_zero()

    def test_norm(self, algebra_3d):
        """测试范数"""
        e1 = algebra_3d.basis_vector(0)
        e2 = algebra_3d.basis_vector(1)

        v = 3 * e1 + 4 * e2

        norm_sq = v.norm_squared()
        norm = v.norm()

        assert abs(norm_sq - 25.0) < 1e-10
        assert abs(norm - 5.0) < 1e-10

    def test_inverse(self, algebra_3d):
        """测试逆"""
        e1 = algebra_3d.basis_vector(0)

        inv = e1.inverse()
        product = e1 * inv

        # e1 * e1⁻¹ = 1
        assert abs(product.scalar_part() - 1.0) < 1e-10

    def test_is_zero(self, algebra_3d):
        """测试零检查"""
        zero = algebra_3d.zeros()
        assert zero.is_zero()

        e1 = algebra_3d.basis_vector(0)
        assert not e1.is_zero()

    def test_is_invertible(self, algebra_3d):
        """测试可逆检查"""
        e1 = algebra_3d.basis_vector(0)
        assert e1.is_invertible()

        zero = algebra_3d.zeros()
        assert not zero.is_invertible()


class TestRotation:
    """旋转测试"""

    @pytest.fixture
    def algebra_2d(self):
        """2D 欧几里得代数"""
        return Algebra.euclidean(2)

    @pytest.fixture
    def algebra_3d(self):
        """3D 欧几里得代数"""
        return Algebra.euclidean(3)

    def test_create_rotor(self, algebra_2d):
        """测试转子创建"""
        e1, e2 = algebra_2d.basis_vectors()
        plane = e1 ^ e2

        rotor = create_rotor(plane, math.pi / 2)

        # 转子应该是单位转子
        rev = rotor.reversion()
        product = rotor * rev
        assert abs(product.scalar_part() - 1.0) < 1e-10

    def test_rotate_90_degrees(self, algebra_2d):
        """测试 90 度旋转"""
        e1, e2 = algebra_2d.basis_vectors()
        plane = e1 ^ e2

        rotor = algebra_2d.rotor(plane, math.pi / 2)
        rotated = e1.rotate_by(rotor)

        # e1 旋转 90 度应该接近 e2
        diff = rotated - e2
        assert diff.norm() < 1e-10

    def test_rotate_180_degrees(self, algebra_2d):
        """测试 180 度旋转"""
        e1 = algebra_2d.basis_vector(0)
        e2 = algebra_2d.basis_vector(1)
        plane = e1 ^ e2

        rotor = algebra_2d.rotor(plane, math.pi)
        rotated = e1.rotate_by(rotor)

        # e1 旋转 180 度应该接近 -e1
        diff = rotated - (-e1)
        assert diff.norm() < 1e-10

    def test_rotate_3d(self, algebra_3d):
        """测试 3D 旋转"""
        e1, e2, e3 = algebra_3d.basis_vectors()
        plane = e1 ^ e2

        rotor = algebra_3d.rotor(plane, math.pi / 2)

        # e1 旋转后应该接近 e2
        rotated_e1 = e1.rotate_by(rotor)
        assert (rotated_e1 - e2).norm() < 1e-10

        # e2 旋转后应该接近 -e1
        rotated_e2 = e2.rotate_by(rotor)
        assert (rotated_e2 - (-e1)).norm() < 1e-10

        # e3 垂直于旋转平面，应该不变
        rotated_e3 = e3.rotate_by(rotor)
        assert (rotated_e3 - e3).norm() < 1e-10


class TestProjection:
    """投影和拒绝测试"""

    @pytest.fixture
    def algebra_3d(self):
        """3D 欧几里得代数"""
        return Algebra.euclidean(3)

    def test_project_to_plane(self, algebra_3d):
        """测试投影到平面"""
        e1, e2, e3 = algebra_3d.basis_vectors()
        plane = e1 ^ e2

        # e1 在 e1-e2 平面内的投影是 e1 本身
        proj_e1 = e1.project_to(plane)
        assert (proj_e1 - e1).norm() < 1e-10

        # e3 在 e1-e2 平面内的投影是 0
        proj_e3 = e3.project_to(plane)
        assert proj_e3.is_zero()

    def test_reject_from_plane(self, algebra_3d):
        """测试拒绝"""
        e1, e2, e3 = algebra_3d.basis_vectors()
        plane = e1 ^ e2

        # e3 相对于 e1-e2 平面的拒绝是 e3 本身
        rej_e3 = e3.reject_from(plane)
        assert (rej_e3 - e3).norm() < 1e-10

        # e1 相对于 e1-e2 平面的拒绝是 0
        rej_e1 = e1.reject_from(plane)
        assert rej_e1.is_zero()


class TestReflection:
    """反射测试"""

    @pytest.fixture
    def algebra_3d(self):
        """3D 欧几里得代数"""
        return Algebra.euclidean(3)

    def test_reflect_in_plane(self, algebra_3d):
        """测试在平面中反射"""
        e1, e2, e3 = algebra_3d.basis_vectors()
        plane = e1 ^ e2

        # e3 在 e1-e2 平面中的反射是 -e3
        refl_e3 = e3.reflect_in(plane)
        assert (refl_e3 - (-e3)).norm() < 1e-10

        # e1 在 e1-e2 平面中的反射是 e1（不变）
        refl_e1 = e1.reflect_in(plane)
        assert (refl_e1 - e1).norm() < 1e-10


class TestScalarProduct:
    """标量积测试"""

    @pytest.fixture
    def algebra_3d(self):
        """3D 欧几里得代数"""
        return Algebra.euclidean(3)

    def test_scalar_product_orthogonal(self, algebra_3d):
        """测试正交向量的标量积"""
        e1, e2 = algebra_3d.basis_vectors()

        sp = e1.scalar_product(e2)
        assert abs(sp) < 1e-10

    def test_scalar_product_same(self, algebra_3d):
        """测试相同向量的标量积"""
        e1 = algebra_3d.basis_vector(0)

        sp = e1.scalar_product(e1)
        assert abs(sp - 1.0) < 1e-10

    def test_scalar_product_general(self, algebra_3d):
        """测试一般向量的标量积"""
        e1, e2 = algebra_3d.basis_vectors()

        a = 3 * e1 + 4 * e2
        b = 1 * e1 + 2 * e2

        sp = a.scalar_product(b)
        # 3*1 + 4*2 = 11
        assert abs(sp - 11.0) < 1e-10


class TestCommutator:
    """交换子测试"""

    @pytest.fixture
    def algebra_3d(self):
        """3D 欧几里得代数"""
        return Algebra.euclidean(3)

    def test_commutator_vectors(self, algebra_3d):
        """测试向量的交换子"""
        e1, e2 = algebra_3d.basis_vectors()

        # [e1, e2] = e1e2 - e2e1 = 2e1e2
        comm = e1.commutator(e2)
        expected = 2 * (e1 ^ e2)

        diff = comm - expected
        assert diff.norm() < 1e-10

    def test_commutator_same(self, algebra_3d):
        """测试相同向量的交换子"""
        e1 = algebra_3d.basis_vector(0)

        comm = e1.commutator(e1)
        assert comm.is_zero()


class TestFromCoefficients:
    """从系数创建多向量测试"""

    @pytest.fixture
    def algebra_2d(self):
        """2D 欧几里得代数"""
        return Algebra.euclidean(2)

    def test_from_coefficients(self, algebra_2d):
        """测试从系数创建"""
        # 在 2D 中：[标量，e1, e2, e12]
        coeffs = [1.0, 2.0, 3.0, 4.0]
        mv = algebra_2d.from_coefficients(coeffs)

        # 检查标量部分
        assert abs(mv.scalar_part() - 1.0) < 1e-10

        # 检查系数
        retrieved_coeffs = mv.coefficients()
        for i, (expected, actual) in enumerate(zip(coeffs, retrieved_coeffs)):
            assert abs(expected - actual) < 1e-10, f"Coefficient {i} mismatch"


class TestGradeProjection:
    """阶次投影测试"""

    @pytest.fixture
    def algebra_3d(self):
        """3D 欧几里得代数"""
        return Algebra.euclidean(3)

    def test_grade_projection(self, algebra_3d):
        """测试阶次投影"""
        # 创建一般多向量
        coeffs = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]
        mv = algebra_3d.from_coefficients(coeffs)

        # 提取各阶次部分
        grade_0 = mv.grade(0)
        grade_1 = mv.grade(1)
        grade_2 = mv.grade(2)
        grade_3 = mv.grade(3)

        # 检查标量部分
        assert abs(grade_0.scalar_part() - 1.0) < 1e-10

        # 检查向量部分
        assert not grade_1.is_zero()

        # 检查二重向量部分
        assert not grade_2.is_zero()

        # 检查三重向量部分
        assert not grade_3.is_zero()

    def test_even_odd_parts(self, algebra_3d):
        """测试偶部和奇部"""
        e1 = algebra_3d.basis_vector(0)
        plane = e1 ^ algebra_3d.basis_vector(1)

        # 偶部：标量 + 二重向量
        even = plane.even_part()
        assert not even.is_zero()

        # 奇部：向量 + 三重向量
        odd = e1.odd_part()
        assert not odd.is_zero()


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
