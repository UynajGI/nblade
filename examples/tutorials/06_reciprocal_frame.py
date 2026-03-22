"""
nblade 互逆标架教程
==================

本教程详细介绍互逆标架的概念、计算方法和应用：
- 什么是互逆标架
- 如何计算互逆标架
- 验证互逆条件
- 度量张量的计算
- 互逆标架的应用

作者: nblade 团队
"""

import nblade
from nblade import reciprocal_frame, verify_reciprocal_frame, metric_tensor


def main():
    """主函数 - 互逆标架教程"""
    print("=" * 70)
    print("nblade 互逆标架教程")
    print("=" * 70)

    section_1_what_is_reciprocal_frame()
    section_2_orthogonal_frame()
    section_3_non_orthogonal_frame()
    section_4_metric_tensor()
    section_5_applications()

    print("\n" + "=" * 70)
    print("教程完成！")
    print("=" * 70)
    print("\n关键点:")
    print("  - 互逆标架满足 aⁱ⌋aⱼ = δⁱⱼ")
    print("  - 正交基的互逆标架就是其本身")
    print("  - 非正交基的互逆标架需要计算")
    print("  - 互逆标架用于在任意标架下分解向量")


def section_1_what_is_reciprocal_frame():
    """第一节：什么是互逆标架"""
    print("\n" + "=" * 70)
    print("第一节：什么是互逆标架")
    print("=" * 70)

    print("""
互逆标架（Reciprocal Frame）是几何代数中的重要概念。

定义：
--------
给定一个标架 {a₁, a₂, ..., aₙ}，其互逆标架 {a¹, a², ..., aⁿ}
定义为满足以下条件的向量集合：

    aⁱ⌋aⱼ = δⁱⱼ

其中：
- ⌋ 表示左收缩（left contraction），对于向量就是内积
- δⁱⱼ 是 Kronecker delta（i=j 时为 1，否则为 0）
- i,j 是指标（上标表示互逆向量，下标表示原向量）

计算公式：
--------
aⁱ = (-1)^(i-1) × (a₁∧...∧âᵢ∧...∧aₙ) × a_N⁻¹

其中：
- âᵢ 表示省略 aᵢ
- a_N = a₁∧a₂∧...∧aₙ 是体积元素
- a_N⁻¹ 是体积元素的逆

直观理解：
--------
互逆向量 aⁱ 与原向量 aᵢ "配对"，与其他原向量正交。
就像标准正交基 {e₁, e₂, e₃} 满足 eᵢ·eⱼ = δᵢⱼ，
互逆标架将这种性质推广到任意（包括非正交）标架。
""")


def section_2_orthogonal_frame():
    """第二节：正交标架的互逆标架"""
    print("\n" + "=" * 70)
    print("第二节：正交标架的互逆标架")
    print("=" * 70)

    print("""
在正交归一标架中，互逆标架就是原标架本身。
这是最简单的情况，用于理解互逆标架的概念。
""")

    # 创建 3D 欧几里得几何代数
    alg = nblade.Algebra.euclidean(3)
    e1, e2, e3 = alg.basis_vectors()

    print("【正交归一基】")
    print(f"标准正交基: e₁ = {e1}")
    print(f"            e₂ = {e2}")
    print(f"            e₃ = {e3}")

    # 计算互逆标架
    frame = [e1, e2, e3]
    reciprocal = reciprocal_frame(frame)

    print("\n【互逆标架】")
    print(f"a¹ = {reciprocal[0]}")
    print(f"a² = {reciprocal[1]}")
    print(f"a³ = {reciprocal[2]}")

    print("\n【验证 aⁱ⌋aⱼ = δⁱⱼ】")
    # 验证互逆条件
    for i in range(3):
        for j in range(3):
            contraction = reciprocal[i] | frame[j]
            value = contraction.scalar_part()
            expected = 1.0 if i == j else 0.0
            status = "✓" if abs(value - expected) < 1e-10 else "✗"
            print(f"  a^{i + 1}⌋a_{j + 1} = {value:6.3f} (期望 {expected}) {status}")

    # 使用 verify_reciprocal_frame 函数
    print("\n【使用 verify_reciprocal_frame() 函数】")
    try:
        verify_reciprocal_frame(frame, reciprocal, 1e-10)
        print("验证通过: 满足 aⁱ⌋aⱼ = δⁱⱼ")
    except Exception as e:
        print(f"验证失败: {e}")

    print("\n结论：")
    print("  对于正交归一基，互逆标架 = 原标架")


def section_3_non_orthogonal_frame():
    """第三节：非正交标架的互逆标架"""
    print("\n" + "=" * 70)
    print("第三节：非正交标架的互逆标架")
    print("=" * 70)

    print("""
在非正交标架中，互逆标架需要计算得到。
本节演示如何创建非正交标架，并计算其互逆标架。
""")

    # 创建 3D 欧几里得几何代数
    alg = nblade.Algebra.euclidean(3)
    e1, e2, e3 = alg.basis_vectors()

    print("【创建非正交标架】")
    print("定义三个非正交的向量:")
    print("  a₁ = e₁           (沿 x 轴)")
    print("  a₂ = 0.5e₁ + e₂   (倾斜于 xy 平面)")
    print("  a₃ = e₃           (沿 z 轴)")

    # 创建非正交标架
    a1 = e1
    a2 = e1.scale(0.5) + e2
    a3 = e3

    print(f"\n原标架:")
    print(f"  a₁ = {a1}")
    print(f"  a₂ = {a2}")
    print(f"  a₃ = {a3}")

    # 验证非正交性
    print("\n【验证非正交性】")
    print("计算 a₁·a₂:", (a1 | a2).scalar_part())
    print("期望: 0.5 (非零，说明不正交)")

    # 计算互逆标架
    frame = [a1, a2, a3]
    reciprocal = reciprocal_frame(frame)

    print("\n【计算互逆标架】")
    print(f"a¹ = {reciprocal[0]}")
    print(f"a² = {reciprocal[1]}")
    print(f"a³ = {reciprocal[2]}")

    # 验证互逆条件
    print("\n【验证 aⁱ⌋aⱼ = δⁱⱼ】")
    try:
        verify_reciprocal_frame(frame, reciprocal, 1e-10)
        print("  验证通过: 满足互逆条件")
    except Exception as e:
        print(f"  验证失败: {e}")

    # 详细验证
    print("\n【详细验证】")
    for i in range(3):
        for j in range(3):
            contraction = reciprocal[i] | frame[j]
            value = contraction.scalar_part()
            expected = 1.0 if i == j else 0.0
            print(f"  a^{i + 1}⌋a_{j + 1} = {value:6.3f}")

    # 另一个例子：缩放向量
    print("\n【例子：缩放向量的互逆标架】")
    print("如果 a₁ = 2e₁, a₂ = 3e₂, a₃ = e₃")
    print("那么互逆标架应该是:")
    print("  a¹ = (1/2)e₁, a² = (1/3)e₂, a³ = e₃")

    b1 = e1.scale(2.0)
    b2 = e2.scale(3.0)
    b3 = e3

    frame2 = [b1, b2, b3]
    reciprocal2 = reciprocal_frame(frame2)

    print(f"\n实际计算结果:")
    print(f"  a¹ = {reciprocal2[0]}")
    print(f"  a² = {reciprocal2[1]}")
    print(f"  a³ = {reciprocal2[2]}")


def section_4_metric_tensor():
    """第四节：度量张量"""
    print("\n" + "=" * 70)
    print("第四节：度量张量")
    print("=" * 70)

    print("""
度量张量 gᵢⱼ = aᵢ·aⱼ 描述了标架的几何性质。
对于非正交标架，度量张量不是单位矩阵。
""")

    # 创建 3D 几何代数
    alg = nblade.Algebra.euclidean(3)
    e1, e2, e3 = alg.basis_vectors()

    print("【例子1：正交归一标架的度量张量】")
    frame_ortho = [e1, e2, e3]
    g_ortho = metric_tensor(frame_ortho)

    print("度量张量 gᵢⱼ:")
    print(f"  | {g_ortho[0][0]:6.3f}  {g_ortho[0][1]:6.3f}  {g_ortho[0][2]:6.3f} |")
    print(f"  | {g_ortho[1][0]:6.3f}  {g_ortho[1][1]:6.3f}  {g_ortho[1][2]:6.3f} |")
    print(f"  | {g_ortho[2][0]:6.3f}  {g_ortho[2][1]:6.3f}  {g_ortho[2][2]:6.3f} |")
    print("  这是单位矩阵，因为基是正交归一的")

    print("\n【例子2：非正交标架的度量张量】")
    print("标架: a₁ = e₁, a₂ = 0.5e₁ + e₂, a₃ = e₃")

    a1 = e1
    a2 = e1.scale(0.5) + e2
    a3 = e3

    frame_non_ortho = [a1, a2, a3]
    g_non_ortho = metric_tensor(frame_non_ortho)

    print("度量张量 gᵢⱼ = aᵢ·aⱼ:")
    print(
        f"  | {g_non_ortho[0][0]:6.3f}  {g_non_ortho[0][1]:6.3f}  {g_non_ortho[0][2]:6.3f} |"
    )
    print(
        f"  | {g_non_ortho[1][0]:6.3f}  {g_non_ortho[1][1]:6.3f}  {g_non_ortho[1][2]:6.3f} |"
    )
    print(
        f"  | {g_non_ortho[2][0]:6.3f}  {g_non_ortho[2][1]:6.3f}  {g_non_ortho[2][2]:6.3f} |"
    )

    print("\n解释:")
    print(f"  g₁₁ = a₁·a₁ = {g_non_ortho[0][0]} (a₁ 的长度平方)")
    print(f"  g₁₂ = a₁·a₂ = {g_non_ortho[0][1]} (a₁ 和 a₂ 的内积)")
    print(f"  g₂₂ = a₂·a₂ = {g_non_ortho[1][1]} (a₂ 的长度平方 = 1 + 0.5² = 1.25)")

    print("\n【度量张量与互逆标架的关系】")
    print("互逆标架满足: aⁱ = gⁱʲ aⱼ")
    print("其中 gⁱʲ 是度量张量的逆矩阵")
    print("这说明了互逆标架与度量张量的密切联系")


def section_5_applications():
    """第五节：互逆标架的应用"""
    print("\n" + "=" * 70)
    print("第五节：互逆标架的应用")
    print("=" * 70)

    print("""
互逆标架在几何代数中有许多重要应用：

1. 在任意标架下分解向量
2. 坐标系转换
3. 物理中的共变和逆变向量
4. 晶体学中的倒易晶格
""")

    # 创建几何代数
    alg = nblade.Algebra.euclidean(3)
    e1, e2, e3 = alg.basis_vectors()

    print("【应用：在任意标架下分解向量】")
    print("\n问题：给定一个非正交标架 {a₁, a₂, a₃}，")
    print("      如何将任意向量 v 表示为 v = v¹a₁ + v²a₂ + v³a₃？")

    # 创建非正交标架
    a1 = e1
    a2 = e1.scale(0.5) + e2
    a3 = e3

    print(f"\n定义标架:")
    print(f"  a₁ = {a1}")
    print(f"  a₂ = {a2}")
    print(f"  a₃ = {a3}")

    # 计算互逆标架
    frame = [a1, a2, a3]
    reciprocal = reciprocal_frame(frame)

    print(f"\n互逆标架:")
    print(f"  a¹ = {reciprocal[0]}")
    print(f"  a² = {reciprocal[1]}")
    print(f"  a³ = {reciprocal[2]}")

    # 定义一个待分解的向量
    v = alg.vector([3.0, 2.0, 1.0])
    print(f"\n待分解的向量 v = {v}")
    print(f"在标准基下: v = 3e₁ + 2e₂ + e₃")

    # 使用互逆标架分解
    print("\n【使用互逆标架计算系数】")
    print("公式: vⁱ = aⁱ⌋v")

    coeffs = []
    for i in range(3):
        coeff = (reciprocal[i] | v).scalar_part()
        coeffs.append(coeff)
        print(f"  v^{i + 1} = a^{i + 1}⌋v = {coeff:.4f}")

    print(f"\n所以: v = {coeffs[0]:.4f}a₁ + {coeffs[1]:.4f}a₂ + {coeffs[2]:.4f}a₃")

    # 验证分解
    print("\n【验证分解】")
    reconstruction = a1.scale(coeffs[0]) + a2.scale(coeffs[1]) + a3.scale(coeffs[2])
    print(
        f"重构: {coeffs[0]:.4f}a₁ + {coeffs[1]:.4f}a₂ + {coeffs[2]:.4f}a₃ = {reconstruction}"
    )
    print(f"原始: v = {v}")

    diff = reconstruction - v
    error = diff.norm()
    print(f"误差: {error:.2e}")

    if error < 1e-10:
        print("✓ 验证成功!")
    else:
        print("✗ 验证失败")

    print("\n【其他应用】")
    print("1. 晶体学：倒易晶格就是互逆标架")
    print("2. 广义相对论：共变和逆变向量的转换")
    print("3. 数值分析：曲线坐标系中的基向量")
    print("4. 计算机图形学：非正交坐标系的变换")


if __name__ == "__main__":
    main()
