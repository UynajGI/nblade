"""
nblade 对偶运算教程
==================

详细介绍几何代数中的对偶运算（Hodge 对偶）：
- 对偶运算的基本概念和计算方法
- 对偶与叉积的关系
- 逆对偶运算
- 不同维度下的对偶特性

对偶运算将 r 阶元素映射到 (n-r) 阶元素，其中 n 是空间维度。
在几何代数中，对偶具有重要的几何意义：
- 向量的对偶表示其"正交补空间"
- 在 3D 中，叉积可以通过对偶与外积的组合表示

作者: nblade 团队
"""

import nblade
import math


def section_1_basic_dual():
    """第一节：基本对偶运算"""
    print("\n" + "=" * 60)
    print("第一节：基本对偶运算")
    print("=" * 60)

    alg = nblade.Algebra.euclidean(3)
    e1, e2, e3 = alg.basis_vectors()

    print("\n对偶运算（Hodge Dual）将 r 阶元素映射到 (n-r) 阶元素")
    print("在 n 维空间中，对偶运算将元素的「阶次」与其「补」互换")

    # 向量的对偶
    print("\n【向量的对偶 → 二重向量】")
    v = alg.vector([1.0, 2.0, 3.0])
    v_dual = v.dual()
    print(f"向量 v = {v}")
    print(f"dual(v) = {v_dual}")
    print("  在 3D 中，向量的对偶是二重向量（表示垂直于原向量的平面）")

    # 基向量的对偶
    print("\n【基向量的对偶】")
    print(f"dual(e1) = {e1.dual()}")
    print(f"dual(e2) = {e2.dual()}")
    print(f"dual(e3) = {e3.dual()}")
    print("  注意: dual(e1) = e2∧e3, dual(e2) = e3∧e1, dual(e3) = e1∧e2")

    # 二重向量的对偶
    print("\n【二重向量的对偶 → 向量】")
    B = e1 ^ e2  # xy 平面
    B_dual = B.dual()
    print(f"二重向量 B = e1∧e2 = {B}")
    print(f"dual(B) = {B_dual}")
    print("  在 3D 中，二重向量的对偶是向量（表示平面的法向量）")

    # 三重向量的对偶
    print("\n【三重向量的对偶 → 标量】")
    I = e1 ^ e2 ^ e3  # 伪标量
    I_dual = I.dual()
    print(f"伪标量 I = e1∧e2∧e3 = {I}")
    print(f"dual(I) = {I_dual}")
    print("  在 3D 中，伪标量的对偶是标量")

    # 标量的对偶
    print("\n【标量的对偶 → 伪标量】")
    s = alg.scalar(2.0)
    s_dual = s.dual()
    print(f"标量 s = {s}")
    print(f"dual(s) = {s_dual}")
    print("  在 3D 中，标量的对偶是伪标量（整个空间的体积元素）")


def section_2_cross_product_via_dual():
    """第二节：通过对偶计算叉积"""
    print("\n" + "=" * 60)
    print("第二节：对偶与叉积的关系")
    print("=" * 60)

    alg = nblade.Algebra.euclidean(3)
    e1, e2, e3 = alg.basis_vectors()

    print("\n关键公式：在 3D 中，a × b = (a ∧ b)*")
    print("其中 * 表示对偶运算")
    print("叉积可以通过外积后再取对偶来实现")

    # 示例 1: e1 × e2
    print("\n【示例 1: e1 × e2】")
    a = e1
    b = e2

    wedge = a ^ b
    cross_via_dual = wedge.dual()

    print(f"a = {a}")
    print(f"b = {b}")
    print(f"a ∧ b = {wedge}")
    print(f"(a ∧ b)* = {cross_via_dual}")
    print(f"  结果: e1 × e2 = e3")

    # 示例 2: e2 × e3
    print("\n【示例 2: e2 × e3】")
    a = e2
    b = e3

    wedge = a ^ b
    cross_via_dual = wedge.dual()

    print(f"a = {a}")
    print(f"b = {b}")
    print(f"a ∧ b = {wedge}")
    print(f"(a ∧ b)* = {cross_via_dual}")
    print(f"  结果: e2 × e3 = e1")

    # 示例 3: 一般向量
    print("\n【示例 3: 一般向量的叉积】")
    a = alg.vector([1.0, 2.0, 0.0])
    b = alg.vector([0.0, 1.0, 1.0])

    # 几何代数方法
    wedge = a ^ b
    cross_ga = wedge.dual()

    # 标准叉积（手动计算）
    ax, ay, az = 1.0, 2.0, 0.0
    bx, by, bz = 0.0, 1.0, 1.0
    cx = ay * bz - az * by
    cy = az * bx - ax * bz
    cz = ax * by - ay * bx
    cross_std = alg.vector([cx, cy, cz])

    print(f"向量 a = {a}")
    print(f"向量 b = {b}")
    print(f"a ∧ b = {wedge}")
    print(f"(a ∧ b)* (几何代数) = {cross_ga}")
    print(f"a × b (标准叉积)   = {cross_std}")

    # 验证
    diff = cross_ga - cross_std
    print(f"\n验证: 两者相等，误差范数 = {diff.norm():.2e}")

    # 叉积的几何意义
    print("\n【叉积的几何意义】")
    print("  a × b 是一个向量，满足:")
    print("  1. 方向：垂直于 a 和 b 张成的平面")
    print("  2. 大小：|a × b| = |a||b|sin(θ) = |a ∧ b|")
    print("  3. 遵循右手定则")

    # 面积关系
    area_bivector = a ^ b
    area_scalar = area_bivector.norm()
    cross_norm = cross_ga.norm()

    print(f"\n  |a ∧ b| (二重向量范数) = {area_scalar:.6f}")
    print(f"  |a × b| (叉积范数)     = {cross_norm:.6f}")
    print(f"  两者相等，验证通过")


def section_3_inverse_dual():
    """第三节：逆对偶运算"""
    print("\n" + "=" * 60)
    print("第三节：逆对偶运算")
    print("=" * 60)

    alg = nblade.Algebra.euclidean(3)
    e1, e2, e3 = alg.basis_vectors()

    print("\n逆对偶运算是 dual() 的逆运算")
    print("关系：a.dual().inverse_dual() == a")
    print("在 3D 欧几里得空间中，dual 和 inverse_dual 相差一个符号")

    # 向量的逆对偶
    print("\n【向量的逆对偶】")
    v = alg.vector([1.0, 2.0, 3.0])
    v_dual = v.dual()
    v_back = v_dual.inverse_dual()

    print(f"原向量 v = {v}")
    print(f"dual(v) = {v_dual}")
    print(f"inverse_dual(dual(v)) = {v_back}")

    # 验证
    diff = v - v_back
    print(f"\n验证: v - inverse_dual(dual(v)) = {diff}")
    print(f"误差范数 = {diff.norm():.2e}")

    # 二重向量的逆对偶
    print("\n【二重向量的逆对偶】")
    B = e1 ^ e2 + 2.0 * (e2 ^ e3)
    B_dual = B.dual()
    B_back = B_dual.inverse_dual()

    print(f"原二重向量 B = {B}")
    print(f"dual(B) = {B_dual}")
    print(f"inverse_dual(dual(B)) = {B_back}")

    diff = B - B_back
    print(f"\n验证: B - inverse_dual(dual(B)) = {diff}")
    print(f"误差范数 = {diff.norm():.2e}")

    # 伪标量的逆对偶
    print("\n【伪标量的逆对偶】")
    I = e1 ^ e2 ^ e3
    I_dual = I.dual()
    I_back = I_dual.inverse_dual()

    print(f"原伪标量 I = {I}")
    print(f"dual(I) = {I_dual}")
    print(f"inverse_dual(dual(I)) = {I_back}")

    diff = I - I_back
    print(f"\n验证: I - inverse_dual(dual(I)) = {diff}")
    print(f"误差范数 = {diff.norm():.2e}")

    # 双重对偶的性质
    print("\n【双重对偶的性质】")
    print("在 n 维空间中，dual(dual(A)) = ±A")
    print("符号取决于 n 和 A 的阶次")

    v_double = v.dual().dual()
    print(f"v.dual().dual() = {v_double}")
    print(f"在 3D 中，向量的双重对偶等于原向量的负值: -v = {(-v)}")

    # 验证双重对偶
    diff_double = v_double + v  # 应为 0，因为 dual(dual(v)) = -v
    print(f"v.dual().dual() + v = {diff_double}")
    print(f"误差范数 = {diff_double.norm():.2e}")


def section_4_dual_in_different_dimensions():
    """第四节：不同维度下的对偶运算"""
    print("\n" + "=" * 60)
    print("第四节：不同维度下的对偶运算")
    print("=" * 60)

    print("\n对偶运算的效果取决于空间的维度")
    print("对偶将 r 阶元素映射到 (n-r) 阶元素")

    # 2D 对偶
    print("\n【2D 欧几里得空间 G(2,0,0)】")
    alg_2d = nblade.Algebra.euclidean(2)
    e1_2d, e2_2d = alg_2d.basis_vectors()

    print(f"基向量: e1, e2")
    print(f"dual(e1) = {e1_2d.dual()}")
    print(f"dual(e2) = {e2_2d.dual()}")
    print("  在 2D 中，向量的对偶是向量（旋转 90 度）")

    v_2d = alg_2d.vector([1.0, 0.0])
    print(f"向量 v = {v_2d}")
    print(f"dual(v) = {v_2d.dual()}")
    print("  dual([1, 0]) = [0, 1]，相当于逆时针旋转 90 度")

    # 验证双重对偶在 2D
    v_2d_double = v_2d.dual().dual()
    print(f"v.dual().dual() = {v_2d_double}")
    print(f"  在 2D 中，双重对偶等于负原向量: -v = {-v_2d}")

    # 3D 对偶
    print("\n【3D 欧几里得空间 G(3,0,0)】")
    alg_3d = nblade.Algebra.euclidean(3)
    e1_3d, e2_3d, e3_3d = alg_3d.basis_vectors()

    v_3d = alg_3d.vector([1.0, 2.0, 3.0])
    print(f"向量 v = {v_3d}")
    print(f"dual(v) = {v_3d.dual()}")
    print("  在 3D 中，向量的对偶是二重向量")

    B_3d = e1_3d ^ e2_3d
    print(f"二重向量 B = {B_3d}")
    print(f"dual(B) = {B_3d.dual()}")
    print("  在 3D 中，二重向量的对偶是向量")

    # 4D 对偶
    print("\n【4D 欧几里得空间 G(4,0,0)】")
    alg_4d = nblade.Algebra.euclidean(4)
    e1_4d, e2_4d, e3_4d, e4_4d = alg_4d.basis_vectors()

    v_4d = alg_4d.vector([1.0, 2.0, 3.0, 4.0])
    print(f"向量 v = {v_4d}")
    print(f"dual(v) = {v_4d.dual()}")
    print("  在 4D 中，向量的对偶是三重向量")

    B_4d = e1_4d ^ e2_4d
    print(f"二重向量 B = e1∧e2 = {B_4d}")
    print(f"dual(B) = {B_4d.dual()}")
    print("  在 4D 中，二重向量的对偶仍是二重向量")

    # 维度与对偶阶次的关系表
    print("\n【维度与对偶阶次关系】")
    print("  维度 | 向量(r=1)对偶 | 二重向量(r=2)对偶")
    print("  -----|---------------|------------------")
    print("  2D   | 向量(r=1)     | 标量(r=0)")
    print("  3D   | 二重向量(r=2) | 向量(r=1)")
    print("  4D   | 三重向量(r=3) | 二重向量(r=2)")


def section_5_volume_element():
    """第五节：对偶与体积元素的关系"""
    print("\n" + "=" * 60)
    print("第五节：对偶与体积元素的关系")
    print("=" * 60)

    alg = nblade.Algebra.euclidean(3)
    e1, e2, e3 = alg.basis_vectors()

    print("\n对偶运算与伪标量（体积元素）密切相关")
    print("对偶 A* = A · I，其中 I 是伪标量（单位体积元素）")

    # 获取体积元素
    config = alg.config
    I = config.volume_element()
    print(f"\n体积元素 I = {I}")
    print(f"I = e1 ∧ e2 ∧ e3 = {e1 ^ e2 ^ e3}")

    # 体积元素的平方
    I_squared = config.volume_element_squared()
    print(f"\nI² = {I_squared}")
    print("  在 3D 欧几里得空间中，I² = -1")

    # 验证对偶公式
    print("\n【验证对偶公式 A* = A · I】")
    v = alg.vector([1.0, 2.0, 3.0])

    # 方法 1: 使用 dual()
    dual_method1 = v.dual()
    # 方法 2: 使用左内积
    dual_method2 = v.left_inner(I)

    print(f"向量 v = {v}")
    print(f"v.dual() = {dual_method1}")
    print(f"v | I = {dual_method2}")

    diff = dual_method1 - dual_method2
    print(f"\n验证: dual(v) = v | I，误差范数 = {diff.norm():.2e}")

    # 逆对偶与体积元素逆的关系
    print("\n【逆对偶与体积元素逆的关系】")
    print("inverse_dual(A) = A · I⁻¹")

    I_inv = config.volume_element_inverse()
    print(f"I⁻¹ = {I_inv}")

    # 方法 1: 使用 inverse_dual()
    inv_dual_method1 = v.dual().inverse_dual()
    # 方法 2: 使用体积元素逆
    inv_dual_method2 = v.dual().left_inner(I_inv)

    print(f"\nv.dual().inverse_dual() = {inv_dual_method1}")
    print(f"dual(v) | I⁻¹ = {inv_dual_method2}")

    # 验证逆对偶返回原向量
    diff_back = v - inv_dual_method1
    print(f"\n验证: v = inverse_dual(dual(v))，误差范数 = {diff_back.norm():.2e}")


def section_6_applications():
    """第六节：对偶运算的应用"""
    print("\n" + "=" * 60)
    print("第六节：对偶运算的应用")
    print("=" * 60)

    alg = nblade.Algebra.euclidean(3)
    e1, e2, e3 = alg.basis_vectors()

    print("\n对偶运算在几何和物理中有广泛应用")

    # 应用 1: 计算法向量
    print("\n【应用 1: 计算平面的法向量】")
    plane = e1 ^ e2  # xy 平面
    normal = plane.dual()
    print(f"平面 B = e1∧e2 = {plane}")
    print(f"法向量 n = dual(B) = {normal}")
    print("  法向量垂直于原平面")

    # 应用 2: 平行四边形面积与法向量
    print("\n【应用 2: 平行四边形的面积和法向量】")
    a = alg.vector([2.0, 0.0, 0.0])
    b = alg.vector([1.0, 1.0, 0.0])

    # 外积得到有向面积（二重向量）
    area_bivector = a ^ b
    # 对偶得到法向量
    normal_vector = area_bivector.dual()
    # 范数得到面积大小
    area = area_bivector.norm()

    print(f"边向量 a = {a}")
    print(f"边向量 b = {b}")
    print(f"有向面积 (二重向量) = {area_bivector}")
    print(f"法向量 = {normal_vector}")
    print(f"面积大小 = {area:.6f}")

    # 应用 3: 判断点是否在平面上
    print("\n【应用 3: 利用对偶判断共面性】")
    v1 = alg.vector([1.0, 0.0, 0.0])
    v2 = alg.vector([0.0, 1.0, 0.0])
    v3 = alg.vector([1.0, 1.0, 0.0])  # 与 v1, v2 共面
    v4 = alg.vector([1.0, 1.0, 1.0])  # 不与 v1, v2 共面

    plane = v1 ^ v2
    print(f"参考平面 = v1∧v2 = {plane}")

    # 向量与平面的内积
    check_v3 = v3 | plane
    check_v4 = v4 | plane

    print(f"v3 | plane = {check_v3}")
    print(f"  v3 在平面上 (结果为 0)")
    print(f"v4 | plane = {check_v4}")
    print(f"  v4 不在平面上 (结果非 0)")


def main():
    """主函数"""
    print("=" * 60)
    print("nblade 对偶运算教程")
    print("Hodge 对偶及其应用")
    print("=" * 60)

    section_1_basic_dual()
    section_2_cross_product_via_dual()
    section_3_inverse_dual()
    section_4_dual_in_different_dimensions()
    section_5_volume_element()
    section_6_applications()

    print("\n" + "=" * 60)
    print("教程完成！")
    print("=" * 60)
    print("\n关键点:")
    print("  - 对偶将 r 阶元素映射到 (n-r) 阶元素")
    print("  - 在 3D 中，a × b = (a ∧ b)*")
    print("  - dual(inverse_dual(A)) = A")
    print("  - 对偶与体积元素相关: A* = A · I")
    print("  - 不同维度下，对偶的效果不同")


if __name__ == "__main__":
    main()
