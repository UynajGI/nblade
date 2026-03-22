#!/usr/bin/env python3
"""
nblade Python 绑定测试套件 / nblade Python Bindings Test Suite

测试 _core Rust 扩展模块和 ga Python 包装器的所有公开 API。
Tests all public APIs of the _core Rust extension module and ga Python wrapper.

运行方式 / How to run:
    maturin develop && pytest python/test__core.py -v
"""

import math
import pytest
import numpy as np

# 导入 Rust 扩展模块 (通过 ga 包)

from nblade._core import AlgebraConfig, MultiVector, create_rotor, __version__

# 导入 Python 包装器
from nblade import Algebra


# =============================================================================
# Fixtures
# =============================================================================


@pytest.fixture
def config_2d():
    """2D 欧几里得代数配置 / 2D Euclidean algebra config"""
    return AlgebraConfig.euclidean(2)


@pytest.fixture
def config_3d():
    """3D 欧几里得代数配置 / 3D Euclidean algebra config"""
    return AlgebraConfig.euclidean(3)


@pytest.fixture
def config_4d():
    """4D 欧几里得代数配置 / 4D Euclidean algebra config"""
    return AlgebraConfig.euclidean(4)


@pytest.fixture
def algebra_3d():
    """3D 欧几里得代数 / 3D Euclidean algebra"""
    return Algebra.euclidean(3)


# =============================================================================
# AlgebraConfig Tests (Rust Extension)
# =============================================================================


class TestAlgebraConfig:
    """测试 AlgebraConfig 类 / Test AlgebraConfig class"""

    def test_euclidean_creation(self):
        """测试欧几里得代数创建 / Test Euclidean algebra creation"""
        config = AlgebraConfig.euclidean(3)
        assert config.dimension == 3
        assert config.basis_count == 8  # 2^3
        assert config.signature == (3, 0, 0)

    def test_custom_signature(self):
        """测试自定义签名创建 / Test custom signature creation"""
        config = AlgebraConfig(4, 2, 1, 1)
        assert config.dimension == 4
        assert config.signature == (2, 1, 1)

    def test_dimension_property(self):
        """测试维度属性 / Test dimension property"""
        for dim in [1, 2, 3, 4, 5, 10]:
            config = AlgebraConfig.euclidean(dim)
            assert config.dimension == dim

    def test_basis_count_property(self):
        """测试基数量属性 / Test basis count property"""
        for dim in [1, 2, 3, 4, 5]:
            config = AlgebraConfig.euclidean(dim)
            assert config.basis_count == 2**dim

    def test_repr(self):
        """测试字符串表示 / Test string representation"""
        config = AlgebraConfig.euclidean(3)
        repr_str = repr(config)
        assert "AlgebraConfig" in repr_str
        assert "dimension=3" in repr_str

    def test_str(self):
        """测试可读字符串 / Test readable string"""
        config = AlgebraConfig.euclidean(3)
        str_repr = str(config)
        # 格式是 "G(3,0,0) Algebra (3D, 8 basis vectors)"
        assert "G(3,0,0)" in str_repr


# =============================================================================
# MultiVector Creation Tests
# =============================================================================


class TestMultiVectorCreation:
    """测试 MultiVector 创建 / Test MultiVector creation"""

    def test_zeros(self, config_3d):
        """测试零多向量创建 / Test zero multivector creation"""
        mv = MultiVector.zeros(config_3d)
        assert mv.is_zero()
        assert mv.scalar_part() == 0.0

    def test_one(self, config_3d):
        """测试单位多向量创建 / Test one multivector creation"""
        mv = MultiVector.one(config_3d)
        assert not mv.is_zero()
        assert mv.scalar_part() == 1.0

    def test_from_scalar(self, config_3d):
        """测试标量多向量创建 / Test scalar multivector creation"""
        mv = MultiVector.from_scalar(config_3d, 5.0)
        assert mv.scalar_part() == 5.0

    def test_basis_vector(self, config_3d):
        """测试基向量创建 / Test basis vector creation"""
        e0 = MultiVector.basis_vector(config_3d, 0)
        e1 = MultiVector.basis_vector(config_3d, 1)
        e2 = MultiVector.basis_vector(config_3d, 2)

        # 基向量应该互相正交 / Basis vectors should be orthogonal
        assert (e0 * e1).scalar_part() == 0.0
        assert (e0 * e2).scalar_part() == 0.0
        assert (e1 * e2).scalar_part() == 0.0

    def test_from_coefficients(self, config_3d):
        """测试从系数创建 / Test creation from coefficients"""
        # 3D 代数有 8 个基 / 3D algebra has 8 basis elements
        coeffs = [1, 2, 3, 4, 5, 6, 7, 8]
        mv = MultiVector.from_coefficients(config_3d, coeffs)
        assert mv.scalar_part() == 1.0

    def test_config_property(self, config_3d):
        """测试配置属性 / Test config property"""
        mv = MultiVector.zeros(config_3d)
        mv_config = mv.config
        assert mv_config.dimension == 3


# =============================================================================
# Product Operations Tests
# =============================================================================


class TestProductOperations:
    """测试乘积运算 / Test product operations"""

    def test_geometric_product_orthogonal(self, config_3d):
        """测试正交向量的几何积 / Test geometric product of orthogonal vectors"""
        e1 = MultiVector.basis_vector(config_3d, 0)
        e2 = MultiVector.basis_vector(config_3d, 1)

        # e1 * e2 = e12 (外积)
        product = e1.geometric_product(e2)
        assert product.scalar_part() == 0.0

    def test_geometric_product_same_vector(self, config_3d):
        """测试同一向量的几何积 / Test geometric product of same vector"""
        e1 = MultiVector.basis_vector(config_3d, 0)

        # e1 * e1 = 1 (欧几里得代数)
        product = e1.geometric_product(e1)
        assert abs(product.scalar_part() - 1.0) < 1e-10

    def test_outer_product(self, config_3d):
        """测试外积 / Test outer product"""
        e1 = MultiVector.basis_vector(config_3d, 0)
        e2 = MultiVector.basis_vector(config_3d, 1)

        wedge = e1.outer_product(e2)
        assert wedge.scalar_part() == 0.0

        # 外积的反交换性: e1 ^ e2 = - (e2 ^ e1)
        wedge_rev = e2.outer_product(e1)
        # e2 ^ e1 = - (e1 ^ e2), 所以相加应该为零
        sum_wedge = wedge + wedge_rev
        assert sum_wedge.is_zero()

    def test_left_inner_orthogonal(self, config_3d):
        """测试正交向量的左内积 / Test left inner product of orthogonal vectors"""
        e1 = MultiVector.basis_vector(config_3d, 0)
        e2 = MultiVector.basis_vector(config_3d, 1)

        # 正交向量内积为零 / Inner product of orthogonal vectors is zero
        inner = e1.left_inner(e2)
        assert inner.is_zero()

    def test_left_inner_same_vector(self, config_3d):
        """测试同一向量的左内积 / Test left inner product of same vector"""
        e1 = MultiVector.basis_vector(config_3d, 0)

        # e1 | e1 = 1
        inner = e1.left_inner(e1)
        assert abs(inner.scalar_part() - 1.0) < 1e-10

    def test_right_inner(self, config_3d):
        """测试右内积 / Test right inner product"""
        e1 = MultiVector.basis_vector(config_3d, 0)

        inner = e1.right_inner(e1)
        assert abs(inner.scalar_part() - 1.0) < 1e-10

    def test_operator_overloads(self, config_3d):
        """测试运算符重载 / Test operator overloads"""
        e1 = MultiVector.basis_vector(config_3d, 0)
        e2 = MultiVector.basis_vector(config_3d, 1)

        # 几何积 *
        product = e1 * e2
        assert product.scalar_part() == 0.0

        # 外积 ^
        wedge = e1 ^ e2
        assert wedge.scalar_part() == 0.0

        # 左内积 |
        inner = e1 | e2
        assert inner.is_zero()


# =============================================================================
# Addition and Subtraction Tests
# =============================================================================


class TestAddSubtract:
    """测试加减运算 / Test addition and subtraction"""

    def test_addition(self, config_3d):
        """测试加法 / Test addition"""
        e1 = MultiVector.basis_vector(config_3d, 0)
        e2 = MultiVector.basis_vector(config_3d, 1)

        sum_mv = e1 + e2
        # 结果应该包含两个向量部分
        assert not sum_mv.is_zero()

    def test_subtraction(self, config_3d):
        """测试减法 / Test subtraction"""
        e1 = MultiVector.basis_vector(config_3d, 0)
        e2 = MultiVector.basis_vector(config_3d, 1)

        diff = e1 - e2
        assert not diff.is_zero()

    def test_negation(self, config_3d):
        """测试取负 / Test negation"""
        e1 = MultiVector.basis_vector(config_3d, 0)
        neg_e1 = -e1

        sum_result = e1 + neg_e1
        assert sum_result.is_zero()


# =============================================================================
# Involution Tests
# =============================================================================


class TestInvolutions:
    """测试对合运算 / Test involution operations"""

    def test_grade_involution_scalar(self, config_3d):
        """测试标量的阶次对合 / Test grade involution of scalar"""
        s = MultiVector.from_scalar(config_3d, 5.0)
        inv = s.grade_involution()
        assert inv.scalar_part() == 5.0  # 标量不变

    def test_grade_involution_vector(self, config_3d):
        """测试向量的阶次对合 / Test grade involution of vector"""
        e1 = MultiVector.basis_vector(config_3d, 0)
        inv = e1.grade_involution()
        # 向量（1阶）变号
        sum_result = e1 + inv
        assert sum_result.is_zero()

    def test_grade_involution_bivector(self, config_3d):
        """测试二重向量的阶次对合 / Test grade involution of bivector"""
        e1 = MultiVector.basis_vector(config_3d, 0)
        e2 = MultiVector.basis_vector(config_3d, 1)
        bivector = e1 ^ e2

        inv = bivector.grade_involution()
        # 二重向量（2阶）不变
        diff = bivector + inv
        # 由于阶次对合对偶数阶不变，应该得到两倍
        assert not diff.is_zero()

    def test_reversion_scalar(self, config_3d):
        """测试标量的反转 / Test reversion of scalar"""
        s = MultiVector.from_scalar(config_3d, 5.0)
        rev = s.reversion()
        assert rev.scalar_part() == 5.0

    def test_reversion_vector(self, config_3d):
        """测试向量的反转 / Test reversion of vector"""
        e1 = MultiVector.basis_vector(config_3d, 0)
        rev = e1.reversion()
        # 向量的反转是自身
        diff = e1 - rev
        assert diff.is_zero()

    def test_reversion_bivector(self, config_3d):
        """测试二重向量的反转 / Test reversion of bivector"""
        e1 = MultiVector.basis_vector(config_3d, 0)
        e2 = MultiVector.basis_vector(config_3d, 1)
        bivector = e1 ^ e2

        rev = bivector.reversion()
        # 二重向量的反转是负的
        sum_result = bivector + rev
        assert sum_result.is_zero()

    def test_clifford_conjugate(self, config_3d):
        """测试 Clifford 共轭 / Test Clifford conjugate"""
        e1 = MultiVector.basis_vector(config_3d, 0)
        conj = e1.clifford_conjugate()
        # 向量的 Clifford 共轭是负的
        sum_result = e1 + conj
        assert sum_result.is_zero()

    def test_invert_operator(self, config_3d):
        """测试取反运算符 / Test invert operator"""
        e1 = MultiVector.basis_vector(config_3d, 0)
        inv = ~e1
        # ~ 是阶次对合
        sum_result = e1 + inv
        assert sum_result.is_zero()


# =============================================================================
# Dual and Inverse Tests
# =============================================================================


class TestDualInverse:
    """测试对偶和逆运算 / Test dual and inverse operations"""

    def test_dual_vector_3d(self, config_3d):
        """测试 3D 向量的对偶 / Test dual of 3D vector"""
        e1 = MultiVector.basis_vector(config_3d, 0)
        dual = e1.dual()

        # 对偶应该是三重向量
        assert not dual.is_zero()

    def test_dual_bivector_3d(self, config_3d):
        """测试 3D 二重向量的对偶 / Test dual of 3D bivector"""
        e1 = MultiVector.basis_vector(config_3d, 0)
        e2 = MultiVector.basis_vector(config_3d, 1)
        bivector = e1 ^ e2

        dual = bivector.dual()
        # e12 的对偶应该是 e3
        assert not dual.is_zero()

    def test_inverse_vector(self, config_3d):
        """测试向量的逆 / Test inverse of vector"""
        e1 = MultiVector.basis_vector(config_3d, 0)

        inv = e1.inverse()
        # e1 * e1^{-1} = 1
        product = e1.geometric_product(inv)
        assert abs(product.scalar_part() - 1.0) < 1e-10

    def test_inverse_scalar(self, config_3d):
        """测试标量的逆 / Test inverse of scalar"""
        s = MultiVector.from_scalar(config_3d, 4.0)

        inv = s.inverse()
        product = s.geometric_product(inv)
        assert abs(product.scalar_part() - 1.0) < 1e-10

    def test_is_invertible(self, config_3d):
        """测试可逆性检查 / Test invertibility check"""
        e1 = MultiVector.basis_vector(config_3d, 0)
        zero = MultiVector.zeros(config_3d)

        assert e1.is_invertible()
        assert not zero.is_invertible()


# =============================================================================
# Norm and Scalar Product Tests
# =============================================================================


class TestNormScalarProduct:
    """测试范数和标量积 / Test norm and scalar product"""

    def test_norm_squared_vector(self, config_3d):
        """测试向量的范数平方 / Test norm squared of vector"""
        e1 = MultiVector.basis_vector(config_3d, 0)
        assert abs(e1.norm_squared() - 1.0) < 1e-10

    def test_norm_vector(self, config_3d):
        """测试向量的范数 / Test norm of vector"""
        e1 = MultiVector.basis_vector(config_3d, 0)
        assert abs(e1.norm() - 1.0) < 1e-10

    def test_norm_scalar(self, config_3d):
        """测试标量的范数 / Test norm of scalar"""
        s = MultiVector.from_scalar(config_3d, 4.0)
        assert abs(s.norm() - 4.0) < 1e-10

    def test_scalar_product_same(self, config_3d):
        """测试同一向量的标量积 / Test scalar product of same vector"""
        e1 = MultiVector.basis_vector(config_3d, 0)
        sp = e1.scalar_product(e1)
        assert abs(sp - 1.0) < 1e-10

    def test_scalar_product_orthogonal(self, config_3d):
        """测试正交向量的标量积 / Test scalar product of orthogonal vectors"""
        e1 = MultiVector.basis_vector(config_3d, 0)
        e2 = MultiVector.basis_vector(config_3d, 1)
        sp = e1.scalar_product(e2)
        assert abs(sp) < 1e-10


# =============================================================================
# Grade Projection Tests
# =============================================================================


class TestGradeProjection:
    """测试阶次投影 / Test grade projection"""

    def test_grade_scalar(self, config_3d):
        """测试标量的阶次投影 / Test grade projection of scalar"""
        s = MultiVector.from_scalar(config_3d, 5.0)
        g0 = s.grade(0)
        assert abs(g0.scalar_part() - 5.0) < 1e-10

    def test_grade_vector(self, config_3d):
        """测试向量的阶次投影 / Test grade projection of vector"""
        e1 = MultiVector.basis_vector(config_3d, 0)
        g1 = e1.grade(1)
        # 向量是 1 阶
        assert not g1.is_zero()

    def test_even_part(self, config_3d):
        """测试偶部 / Test even part"""
        s = MultiVector.from_scalar(config_3d, 5.0)
        even = s.even_part()
        assert abs(even.scalar_part() - 5.0) < 1e-10

    def test_odd_part(self, config_3d):
        """测试奇部 / Test odd part"""
        e1 = MultiVector.basis_vector(config_3d, 0)
        odd = e1.odd_part()
        assert not odd.is_zero()


# =============================================================================
# Geometry Operations Tests
# =============================================================================


class TestGeometryOperations:
    """测试几何运算 / Test geometry operations"""

    def test_project_to_vector(self, config_3d):
        """测试向量投影 / Test vector projection"""
        e1 = MultiVector.basis_vector(config_3d, 0)
        e2 = MultiVector.basis_vector(config_3d, 1)

        # 将 e2 投影到 e1 上 / Project e2 onto e1
        proj = e2.project_to(e1)
        assert proj.is_zero()  # 正交向量投影为零

    def test_project_to_self(self, config_3d):
        """测试自身投影 / Test projection onto self"""
        e1 = MultiVector.basis_vector(config_3d, 0)

        # e1 在 e1 上的投影应该等于 e1
        proj = e1.project_to(e1)
        # 验证投影结果与原向量相同
        diff = proj - e1
        assert diff.is_zero()

    def test_reject_from_orthogonal(self, config_3d):
        """测试正交向量的拒绝 / Test rejection from orthogonal"""
        e1 = MultiVector.basis_vector(config_3d, 0)
        e2 = MultiVector.basis_vector(config_3d, 1)

        # 从 e1 拒绝 e2 / Reject e2 from e1
        reject = e2.reject_from(e1)
        # 因为正交，拒绝等于原向量
        assert not reject.is_zero()

    def test_reflect_in_plane(self, config_3d):
        """测试平面反射 / Test reflection in plane"""
        e1 = MultiVector.basis_vector(config_3d, 0)
        e2 = MultiVector.basis_vector(config_3d, 1)
        plane = e1  # 用向量作为反射平面

        # 在 e1 中反射 e2 / Reflect e2 in e1
        reflected = e2.reflect_in(plane)
        assert not reflected.is_zero()

    def test_rotate_by_90_degrees(self, config_3d):
        """测试 90 度旋转 / Test 90 degree rotation"""
        e1 = MultiVector.basis_vector(config_3d, 0)
        e2 = MultiVector.basis_vector(config_3d, 1)
        plane = e1 ^ e2

        # 创建 90 度转子 / Create 90 degree rotor
        rotor = create_rotor(plane, math.pi / 2)

        # 旋转 e1 / Rotate e1
        rotated = e1.rotate_by(rotor)
        assert not rotated.is_zero()


# =============================================================================
# Commutator Tests
# =============================================================================


class TestCommutator:
    """测试交换子 / Test commutator"""

    def test_commutator_basis_vectors(self, config_3d):
        """测试基向量的交换子 / Test commutator of basis vectors"""
        e1 = MultiVector.basis_vector(config_3d, 0)
        e2 = MultiVector.basis_vector(config_3d, 1)

        comm = e1.commutator(e2)
        # [e1, e2] = e1*e2 - e2*e1
        assert not comm.is_zero()


# =============================================================================
# Python Wrapper (Algebra class) Tests
# =============================================================================


class TestAlgebraWrapper:
    """测试 Python 包装器 Algebra 类 / Test Python wrapper Algebra class"""

    def test_euclidean_creation(self):
        """测试欧几里得代数创建 / Test Euclidean algebra creation"""
        alg = Algebra.euclidean(3)
        assert alg.dimension == 3
        assert alg.signature == (3, 0, 0)

    def test_spacetime_creation(self):
        """测试时空代数创建 / Test spacetime algebra creation"""
        alg = Algebra.spacetime(4)
        assert alg.dimension == 4
        assert alg.signature == (1, 3, 0)

    def test_cga_creation(self):
        """测试共形几何代数创建 / Test conformal GA creation"""
        alg = Algebra.cga()
        assert alg.dimension == 5
        assert alg.signature == (4, 1, 0)

    def test_custom_creation(self):
        """测试自定义代数创建 / Test custom algebra creation"""
        alg = Algebra(4, p=2, q=1, r=1)
        assert alg.dimension == 4
        assert alg.signature == (2, 1, 1)

    def test_basis_vector_method(self, algebra_3d):
        """测试基向量方法 / Test basis vector method"""
        e1 = algebra_3d.basis_vector(0)
        assert not e1.is_zero()

    def test_basis_vectors_method(self, algebra_3d):
        """测试基向量列表方法 / Test basis vectors method"""
        vectors = algebra_3d.basis_vectors()
        assert len(vectors) == 3
        for v in vectors:
            assert not v.is_zero()

    def test_scalar_method(self, algebra_3d):
        """测试标量方法 / Test scalar method"""
        s = algebra_3d.scalar(5.0)
        assert abs(s.scalar_part() - 5.0) < 1e-10

    def test_one_method(self, algebra_3d):
        """测试单位方法 / Test one method"""
        one = algebra_3d.one()
        assert abs(one.scalar_part() - 1.0) < 1e-10

    def test_zeros_method(self, algebra_3d):
        """测试零方法 / Test zeros method"""
        zero = algebra_3d.zeros()
        assert zero.is_zero()

    def test_from_coefficients_method(self):
        """测试系数创建方法 / Test from_coefficients method"""
        alg = Algebra.euclidean(2)
        mv = alg.from_coefficients([1, 2, 3, 4])
        assert abs(mv.scalar_part() - 1.0) < 1e-10

    def test_str_representation(self, algebra_3d):
        """测试字符串表示 / Test string representation"""
        str_repr = str(algebra_3d)
        assert "G(3, 0, 0)" in str_repr or "3D" in str_repr

    def test_repr_representation(self, algebra_3d):
        """测试 repr 表示 / Test repr representation"""
        repr_str = repr(algebra_3d)
        assert "Algebra" in repr_str


# =============================================================================
# Module-Level Tests
# =============================================================================


class TestModule:
    """测试模块级别功能 / Test module-level features"""

    def test_version_exists(self):
        """测试版本号存在 / Test version exists"""
        assert __version__ is not None
        assert isinstance(__version__, str)

    def test_create_rotor_function(self, config_3d):
        """测试 create_rotor 函数 / Test create_rotor function"""
        e1 = MultiVector.basis_vector(config_3d, 0)
        e2 = MultiVector.basis_vector(config_3d, 1)
        plane = e1 ^ e2

        rotor = create_rotor(plane, math.pi / 4)
        assert not rotor.is_zero()


# =============================================================================
# Edge Cases and Error Handling Tests
# =============================================================================


class TestEdgeCases:
    """测试边缘情况 / Test edge cases"""

    def test_1d_algebra(self):
        """测试 1D 代数 / Test 1D algebra"""
        config = AlgebraConfig.euclidean(1)
        assert config.dimension == 1
        assert config.basis_count == 2

    def test_high_dimension_algebra(self):
        """测试高维代数 / Test high dimension algebra"""
        config = AlgebraConfig.euclidean(10)
        assert config.dimension == 10
        assert config.basis_count == 1024

    def test_zero_scalar(self, config_3d):
        """测试零标量 / Test zero scalar"""
        s = MultiVector.from_scalar(config_3d, 0.0)
        assert s.is_zero()

    def test_negative_scalar(self, config_3d):
        """测试负标量 / Test negative scalar"""
        s = MultiVector.from_scalar(config_3d, -5.0)
        assert abs(s.scalar_part() - (-5.0)) < 1e-10

    def test_large_scalar(self, config_3d):
        """测试大标量 / Test large scalar"""
        s = MultiVector.from_scalar(config_3d, 1e10)
        assert abs(s.scalar_part() - 1e10) < 1e-5

    def test_small_scalar(self, config_3d):
        """测试小标量 / Test small scalar"""
        s = MultiVector.from_scalar(config_3d, 1e-10)
        assert abs(s.scalar_part() - 1e-10) < 1e-15


# =============================================================================
# Integration Tests
# =============================================================================


class TestIntegration:
    """集成测试 / Integration tests"""

    def test_full_computation_pipeline(self, algebra_3d):
        """测试完整计算流程 / Test full computation pipeline"""
        # 创建向量
        e1, e2, e3 = algebra_3d.basis_vectors()

        # 创建平面
        plane = e1 ^ e2

        # 创建转子
        rotor = algebra_3d.rotor(plane, math.pi / 2)

        # 旋转向量
        rotated = e1.rotate_by(rotor)

        # 验证旋转结果
        assert not rotated.is_zero()

    def test_geometric_product_identity(self, config_3d):
        """测试几何积恒等式 / Test geometric product identity"""
        e1 = MultiVector.basis_vector(config_3d, 0)
        e2 = MultiVector.basis_vector(config_3d, 1)

        # ab = a·b + a∧b 对于正交向量
        geometric = e1.geometric_product(e2)
        inner = e1.left_inner(e2)
        outer = e1.outer_product(e2)

        # 对于正交向量，内积为零，几何积等于外积
        assert inner.is_zero()


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
