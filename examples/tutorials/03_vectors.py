"""
nblade 向量操作教程
===================

详细介绍向量的创建、操作和运算：
- 多种创建向量的方法
- 向量的基本属性
- 向量的代数运算
- 向量的几何操作

作者: nblade 团队
"""

import nblade
import math


def section_1_creating_vectors():
    """第一节：创建向量"""
    print("\n" + "=" * 60)
    print("第一节：创建向量")
    print("=" * 60)

    alg = nblade.Algebra.euclidean(3)

    # 方法1: 从列表创建
    print("\n【方法1: 从 Python 列表创建】")
    v1 = alg.vector([1.0, 2.0, 3.0])
    print(f"v1 = alg.vector([1, 2, 3]) = {v1}")

    # 方法2: 从元组创建
    print("\n【方法2: 从 Python 元组创建】")
    v2 = alg.vector((4.0, 5.0, 6.0))
    print(f"v2 = alg.vector((4, 5, 6)) = {v2}")

    # 方法3: 获取基向量
    print("\n【方法3: 获取基向量】")
    e1 = alg.basis_vector(0)  # 第 0 个基向量
    e2 = alg.basis_vector(1)  # 第 1 个基向量
    e3 = alg.basis_vector(2)  # 第 2 个基向量
    print(f"e1 = alg.basis_vector(0) = {e1}")
    print(f"e2 = alg.basis_vector(1) = {e2}")
    print(f"e3 = alg.basis_vector(2) = {e3}")

    # 方法4: 获取所有基向量
    print("\n【方法4: 获取所有基向量】")
    basis = alg.basis_vectors()
    print(f"basis = alg.basis_vectors() = {len(basis)} 个基向量")
    for i, b in enumerate(basis):
        print(f"  e{i + 1} = {b}")

    # 方法5: 创建零向量
    print("\n【方法5: 创建零向量/标量】")
    zero = alg.zeros()
    print(f"zero = alg.zeros() = {zero}")

    one = alg.one()
    print(f"one = alg.one() = {one}")

    scalar = alg.scalar(3.14)
    print(f"scalar = alg.scalar(3.14) = {scalar}")


def section_2_vector_properties():
    """第二节：向量属性"""
    print("\n" + "=" * 60)
    print("第二节：向量属性")
    print("=" * 60)

    alg = nblade.Algebra.euclidean(3)
    v = alg.vector([3.0, 4.0, 0.0])

    print(f"\n向量 v = {v}")

    # 范数
    print("\n【范数】")
    norm = v.norm()
    norm_sq = v.norm_squared()
    print(f"范数 |v| = {norm:.6f}")
    print(f"范数平方 |v|² = {norm_sq:.6f}")

    # 标量部分
    print("\n【标量部分】")
    scalar = v.scalar_part()
    print(f"标量部分 = {scalar} (向量没有标量部分)")

    # 系数
    print("\n【系数数组】")
    coeffs = v.coefficients()
    print(f"系数 = {coeffs}")
    print("  (包含所有基的系数，大部分为 0)")

    # 阶次
    print("\n【阶次】")
    grade_1 = v.grade(1)
    print(f"1-阶次部分 = {grade_1} (向量本身)")
    grade_0 = v.grade(0)
    print(f"0-阶次部分 = {grade_0} (标量部分为空)")


def section_3_algebraic_operations():
    """第三节：代数运算"""
    print("\n" + "=" * 60)
    print("第三节：代数运算")
    print("=" * 60)

    alg = nblade.Algebra.euclidean(3)
    e1, e2, e3 = alg.basis_vectors()

    a = alg.vector([1.0, 2.0, 3.0])
    b = alg.vector([4.0, 5.0, 6.0])

    print(f"a = {a}")
    print(f"b = {b}")

    # 加法
    print("\n【加法】")
    sum_ab = a + b
    print(f"a + b = {sum_ab}")

    # 减法
    print("\n【减法】")
    diff_ab = a - b
    print(f"a - b = {diff_ab}")

    # 标量乘法
    print("\n【标量乘法】")
    scaled = a.scale(2.0)
    print(f"a.scale(2.0) = {scaled}")

    # 负号
    print("\n【取负】")
    neg = -a
    print(f"-a = {neg}")

    # 逆
    print("\n【逆向量】")
    inv = a.inverse()
    print(f"a⁻¹ = {inv}")
    print(f"验证 a * a⁻¹ = {(a * inv).scalar_part():.6f}")

    # 三种乘积
    print("\n【三种乘积】")
    geom = a * b
    outer = a ^ b
    inner = a | b

    print(f"几何积 a * b = {geom}")
    print(f"外积 a ^ b   = {outer}")
    print(f"内积 a | b   = {inner}")


def section_4_geometric_operations():
    """第四节：几何操作"""
    print("\n" + "=" * 60)
    print("第四节：几何操作")
    print("=" * 60)

    alg = nblade.Algebra.euclidean(3)
    e1, e2, e3 = alg.basis_vectors()

    # 投影
    print("\n【投影】")
    a = alg.vector([1.0, 2.0, 3.0])
    b = e1  # 投影到 x 轴

    proj = a.project_to(b)
    print(f"a = {a}")
    print(f"b = {b} (e1, x 轴方向)")
    print(f"a 在 b 上的投影 = {proj}")

    # 拒绝 (正交分量)
    print("\n【拒绝 (正交分量)】")
    reject = a.reject_from(b)
    print(f"a 在 b 上的拒绝 = {reject}")
    print(f"验证: proj + reject = {proj + reject}")

    # 反射
    print("\n【反射】")
    v = e1 + e2
    normal = e1  # 反射面法向量

    reflected = v.reflect_in(normal)
    print(f"v = {v}")
    print(f"反射面法向量 = {normal}")
    print(f"反射结果 = {reflected}")

    # 对偶
    print("\n【对偶】")
    v = alg.vector([1.0, 2.0, 3.0])
    dual_v = v.dual()
    print(f"v = {v}")
    print(f"dual(v) = {dual_v}")
    print("  向量的对偶是二重向量")

    # 逆对偶
    inv_dual = dual_v.inverse_dual()
    print(f"inverse_dual(dual(v)) = {inv_dual}")


def section_5_involutions():
    """第五节：对合运算"""
    print("\n" + "=" * 60)
    print("第五节：对合运算")
    print("=" * 60)

    alg = nblade.Algebra.euclidean(3)
    e1, e2, e3 = alg.basis_vectors()

    print("\n对合运算是将多向量映射到自身的运算")
    print("对于向量（奇阶次），三种对合的效果相同或相似")

    v = alg.vector([1.0, 2.0, 3.0])
    print(f"\n向量 v = {v}")

    # 阶次对合
    print("\n【阶次对合 (Grade Involution)】")
    # 对于 r 阶次元素，阶次对合乘以 (-1)^r
    # 向量是 1 阶次，所以乘以 (-1)^1 = -1
    involution = v.grade_involution()
    print(f"grade_involution(v) = {involution}")
    print("  公式: (-1)^r * v_r，向量 r=1，结果为 -v")

    # 反转
    print("\n【反转 (Reversion)】")
    # 对于 r 阶次元素，反转乘以 (-1)^(r(r-1)/2)
    # 向量 r=1，r(r-1)/2 = 0，所以乘以 1
    reversion = v.reversion()
    print(f"reversion(v) = {reversion}")
    print("  公式: (-1)^(r(r-1)/2) * v_r，向量 r=1，结果不变")

    # Clifford 共轭
    print("\n【Clifford 共轭】")
    # Clifford 共轭 = 阶次对合 + 反转
    conjugate = v.clifford_conjugate()
    print(f"clifford_conjugate(v) = {conjugate}")

    # 偶部和奇部
    print("\n【偶部和奇部】")
    print(f"even_part(v) = {v.even_part()}")
    print(f"odd_part(v) = {v.odd_part()}")

    # 对于混合多向量的对合效果
    print("\n【混合多向量的对合】")
    mixed = alg.scalar(1.0) + e1 + (e1 ^ e2) + (e1 ^ e2 ^ e3)
    print(f"混合多向量: {mixed}")
    print(f"阶次对合:   {mixed.grade_involution()}")
    print(f"反转:       {mixed.reversion()}")
    print(f"Clifford共轭: {mixed.clifford_conjugate()}")


def section_6_higher_dimensions():
    """第六节：高维向量"""
    print("\n" + "=" * 60)
    print("第六节：高维向量")
    print("=" * 60)

    # 5D 向量
    print("\n【5D 欧几里得向量】")
    alg_5d = nblade.Algebra.euclidean(5)
    v5 = alg_5d.vector([1, 2, 3, 4, 5])
    print(f"代数: {alg_5d}")
    print(f"向量: {v5}")

    basis5 = alg_5d.basis_vectors()
    print(f"基向量数量: {len(basis5)}")

    # 高维外积
    print("\n【高维外积】")
    v1 = alg_5d.vector([1, 0, 0, 0, 0])
    v2 = alg_5d.vector([0, 1, 0, 0, 0])
    v3 = alg_5d.vector([0, 0, 1, 0, 0])

    bivector = v1 ^ v2
    trivector = v1 ^ v2 ^ v3
    print(f"二重向量 v1∧v2 = {bivector}")
    print(f"三重向量 v1∧v2∧v3 = {trivector}")

    # 不同签名
    print("\n【时空代数 (Spacetime Algebra)】")
    alg_st = nblade.Algebra.spacetime(4)  # G(1,3,0)
    print(f"代数: {alg_st}")
    print("  1 个时间基向量 (平方 = +1)")
    print("  3 个空间基向量 (平方 = -1)")

    e0 = alg_st.basis_vector(0)  # 时间方向
    e1 = alg_st.basis_vector(1)  # 空间方向

    print(f"\ne0² = {(e0 * e0).scalar_part()}")  # 应为 +1
    print(f"e1² = {(e1 * e1).scalar_part()}")  # 应为 -1


def main():
    """主函数"""
    print("=" * 60)
    print("nblade 向量操作教程")
    print("=" * 60)

    section_1_creating_vectors()
    section_2_vector_properties()
    section_3_algebraic_operations()
    section_4_geometric_operations()
    section_5_involutions()
    section_6_higher_dimensions()

    print("\n" + "=" * 60)
    print("教程完成！")
    print("=" * 60)
    print("\n关键点:")
    print("  - 使用 alg.vector([...]) 创建向量")
    print("  - 使用 alg.basis_vectors() 获取基向量")
    print("  - 支持所有基本代数运算和几何操作")
    print("  - 支持任意维度 (最多 64D)")


if __name__ == "__main__":
    main()
