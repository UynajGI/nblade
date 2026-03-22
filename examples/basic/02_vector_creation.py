"""
从 Python 列表/数组创建向量示例

演示如何从 Python 列表、元组、NumPy 数组创建几何代数向量
"""

from nblade import Algebra
import math


def demo_vector_creation():
    """演示从 Python 数据结构创建向量"""
    print("=== 从 Python 列表/数组创建向量 ===\n")

    # 创建 3D 欧几里得代数
    alg = Algebra.euclidean(3)
    print(f"代数: {alg}\n")

    # 从 Python 列表创建向量
    print("1. 从 Python 列表创建向量:")
    v1 = alg.vector([1, 2, 3])
    print(f"   v1 = alg.vector([1, 2, 3]) -> {v1}")

    # 从 Python 元组创建向量
    print("\n2. 从 Python 元组创建向量:")
    v2 = alg.vector((4, 5, 6))
    print(f"   v2 = alg.vector((4, 5, 6)) -> {v2}")

    # 从标量创建标量多向量
    print("\n3. 从标量创建标量多向量:")
    s = alg.scalar(2.5)
    print(f"   s = alg.scalar(2.5) -> {s}")

    # 创建单位多向量
    print("\n4. 创建单位多向量:")
    one = alg.one()
    print(f"   one = alg.one() -> {one}")

    # 创建零多向量
    print("\n5. 创建零多向量:")
    zero = alg.zeros()
    print(f"   zero = alg.zeros() -> {zero}")

    # 获取基向量
    print("\n6. 获取基向量:")
    e1, e2, e3 = alg.basis_vectors()
    print(f"   e1 = {e1}")
    print(f"   e2 = {e2}")
    print(f"   e3 = {e3}")

    # 单独获取基向量
    print("\n7. 单独获取基向量:")
    e1_alt = alg.basis_vector(0)
    e2_alt = alg.basis_vector(1)
    e3_alt = alg.basis_vector(2)
    print(f"   e1_alt = {e1_alt}")
    print(f"   e2_alt = {e2_alt}")
    print(f"   e3_alt = {e3_alt}")

    # 验证相等性
    print("\n8. 验证相等性:")
    print(f"   e1 == e1_alt: {e1 == e1_alt}")
    print(f"   e2 == e2_alt: {e2 == e2_alt}")
    print(f"   e3 == e3_alt: {e3 == e3_alt}")

    # 创建更高维度的向量
    print("\n9. 高维向量创建 (5D):")
    alg_5d = Algebra.euclidean(5)
    v5d = alg_5d.vector([1, 2, 3, 4, 5])
    print(f"   5D 向量: {v5d}")

    # 创建混合多向量 (从系数)
    print("\n10. 从系数创建混合多向量:")
    # 在 3D 中，系数顺序为 [标量, e1, e2, e3, e12, e13, e23, e123]
    coeffs = [1.0, 2.0, 3.0, 4.0, 0.5, 0.0, 0.0, 0.1]
    mixed_mv = alg.from_coefficients(coeffs)
    print(f"   混合多向量: {mixed_mv}")

    print("\n向量创建示例完成！\n")


def demo_operations_with_created_vectors():
    """演示使用创建的向量进行运算"""
    print("=== 使用创建的向量进行运算 ===\n")

    alg = Algebra.euclidean(3)

    # 创建测试向量
    v1 = alg.vector([1, 0, 0])  # e1
    v2 = alg.vector([0, 1, 0])  # e2
    v3 = alg.vector([1, 1, 1])  # 一般向量

    print(f"v1 = {v1}")
    print(f"v2 = {v2}")
    print(f"v3 = {v3}\n")

    # 运算示例
    print("运算示例:")

    # 几何积
    prod = v1 * v2
    print(f"  v1 * v2 = {prod}")

    # 外积
    wedge = v1 ^ v2
    print(f"  v1 ^ v2 = {wedge}")

    # 内积
    inner = v1 | v2
    print(f"  v1 | v2 = {inner}")

    # 向量加法
    sum_v = v1 + v2
    print(f"  v1 + v2 = {sum_v}")

    # 标量乘法
    scaled = v1 * 3.0
    print(f"  v1 * 3.0 = {scaled}")

    # 范数
    norm_v3 = v3.norm()
    print(f"  |v3| = {norm_v3:.6f}")

    # 范数平方
    norm_sq_v3 = v3.norm_squared()
    print(f"  |v3|² = {norm_sq_v3:.6f}")

    # 逆 (如果存在)
    try:
        inv_v1 = v1.inverse()
        print(f"  v1⁻¹ = {inv_v1}")

        # 验证 v1 * v1⁻¹ = 1
        verify_inv = v1 * inv_v1
        print(f"  v1 * v1⁻¹ = {verify_inv}")
    except Exception as e:
        print(f"  v1 逆不存在: {e}")

    print("\n运算示例完成！\n")


if __name__ == "__main__":
    demo_vector_creation()
    demo_operations_with_created_vectors()
