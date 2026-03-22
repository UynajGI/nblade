"""
nblade: N-dimensional Blade - 高性能几何代数 Python 库
=====================================================

nblade 是一个基于 Rust 实现的高性能几何代数库，支持任意维度（最多 64 维）、任意签名 (p, q, r)。

核心特性
--------
- **任意维度**: 支持最多 64 维向量空间
- **任意签名**: 支持 G(p, q, r) 任意度量签名
- **高性能**: 使用 Rust 实现，支持并行计算和 SIMD 优化
- **双表示**: 自动选择密集或稀疏表示以优化性能
- **完整运算**: 实现所有几何代数标准运算
- **NumPy 支持**: 直接从 numpy 数组或 Python 列表创建向量

快速开始
--------
>>> import nblade
>>>
>>> # 创建 3D 欧几里得几何代数 G(3,0,0)
>>> algebra = nblade.Algebra.euclidean(3)
>>>
>>> # 从 Python 列表创建向量
>>> v = algebra.vector([1, 2, 3])  # 等价于 1*e1 + 2*e2 + 3*e3
>>> print(v)
1.0*e1 + 2.0*e2 + 3.0*e3
>>>
>>> # 创建基向量
>>> e1 = algebra.basis_vector(0)
>>> e2 = algebra.basis_vector(1)
>>>
>>> # 几何积
>>> product = e1 * e2
>>> print(product)
1.0*e12
>>>
>>> # 外积（楔积）
>>> wedge = e1 ^ e2
>>> print(wedge)
1.0*e12
>>>
>>> # 左内积
>>> inner = e1 | e2
>>> print(inner)
0.0

数学背景
--------
几何代数（Geometric Algebra，又称 Clifford Algebra）是传统向量代数的推广，
由 William Clifford 在 19 世纪发明。它定义了：

1. **多重向量 (Multivector)**: 包含标量、向量、二重向量（平面）、三重向量（体积）等
2. **几何积 (Geometric Product)**: 统一的乘法运算，包含内积和外积
3. **几何解释**: 每个代数对象都有清晰的几何意义

基本关系式::

    ab = a·b + a∧b

其中：
- ab 是几何积
- a·b 是内积（对称部分）
- a∧b 是外积（反对称部分）

对于正交向量，内积为零，几何积等于外积。

许可证
------
MIT License
"""

from typing import List, Tuple, Union

from ._core import (
    AlgebraConfig,
    MultiVector,
    create_rotor,
    reciprocal_frame,
    verify_reciprocal_frame,
    metric_tensor,
    basis_expansion,
    basis_reconstruction,
)

__version__ = "0.1.0"

__all__ = [
    "Algebra",
    "AlgebraConfig",
    "MultiVector",
    "create_rotor",
    "reciprocal_frame",
    "verify_reciprocal_frame",
    "metric_tensor",
    "basis_expansion",
    "basis_reconstruction",
    "__version__",
]


class Algebra:
    """
    几何代数配置类

    用于创建和配置任意维度、任意签名的几何代数空间。
    支持欧几里得空间、时空代数、共形几何代数等。

    参数
    ----
    dimension : int
        向量空间维度
    p : int, optional
        正平方基向量数量 (e_i² = +1)，默认 0
    q : int, optional
        负平方基向量数量 (e_i² = -1)，默认 0
    r : int, optional
        零平方基向量数量 (e_i² = 0)，默认 0

    示例
    ----
    >>> import nblade
    >>> algebra = nblade.Algebra.euclidean(3)
    >>> print(algebra)
    Geometric Algebra G(3, 0, 0) (3D)
    """

    def __init__(self, dimension: int, p: int = 0, q: int = 0, r: int = 0) -> None:
        self._config = AlgebraConfig(dimension, p, q, r)

    @classmethod
    def euclidean(cls, dimension: int) -> "Algebra":
        """创建欧几里得几何代数 G(n, 0, 0)"""
        return cls(dimension, dimension, 0, 0)

    @classmethod
    def spacetime(cls, dimension: int) -> "Algebra":
        """创建时空代数 G(1, n-1, 0)"""
        return cls(dimension, 1, dimension - 1, 0)

    @classmethod
    def cga(cls) -> "Algebra":
        """创建共形几何代数 G(4, 1, 0)"""
        return cls(5, 4, 1, 0)

    @property
    def config(self) -> AlgebraConfig:
        return self._config

    @property
    def dimension(self) -> int:
        return self._config.dimension

    @property
    def basis_count(self) -> int:
        return self._config.basis_count

    @property
    def signature(self) -> Tuple[int, int, int]:
        return self._config.signature

    def basis_vector(self, i: int) -> MultiVector:
        """创建第 i 个基向量 e_i"""
        return MultiVector.basis_vector(self._config, i)

    def basis_vectors(self) -> List[MultiVector]:
        """获取所有基向量"""
        return [self.basis_vector(i) for i in range(self.dimension)]

    def scalar(self, value: float) -> MultiVector:
        """创建标量多向量"""
        return MultiVector.from_scalar(self._config, value)

    def from_scalar(self, value: float) -> MultiVector:
        """创建标量多向量（scalar 方法的别名）/ Create scalar multivector (alias for scalar)"""
        return self.scalar(value)

    def one(self) -> MultiVector:
        """创建单位多向量（标量 1）"""
        return MultiVector.one(self._config)

    def zeros(self) -> MultiVector:
        """创建零多向量"""
        return MultiVector.zeros(self._config)

    def vector(self, data: Union[List[float], Tuple[float, ...]]) -> MultiVector:
        """从列表或元组创建向量（1-向量）"""
        if len(data) != self.dimension:
            raise ValueError(
                f"Vector data length ({len(data)}) must equal algebra dimension ({self.dimension})"
            )

        result = self.zeros()
        for i, coeff in enumerate(data):
            if abs(coeff) > 1e-15:
                basis_vec = self.basis_vector(i)
                scaled_vec = basis_vec.scale(coeff)
                result = result.add(scaled_vec)
        return result

    def from_coefficients(self, coefficients: List[float]) -> MultiVector:
        """从系数数组创建多向量"""
        return MultiVector.from_coefficients(self._config, coefficients)

    def rotor(self, plane: MultiVector, angle: float) -> MultiVector:
        """创建转子（Rotator）"""
        return create_rotor(plane, angle)

    def __repr__(self) -> str:
        return f"Algebra(dimension={self.dimension}, signature={self.signature})"

    def __str__(self) -> str:
        p, q, r = self.signature
        return f"Geometric Algebra G({p}, {q}, {r}) ({self.dimension}D)"
