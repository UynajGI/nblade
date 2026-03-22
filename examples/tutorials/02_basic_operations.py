"""
nblade 基本运算教程
===================

详细介绍几何代数的三种基本乘积：
1. 几何积 (Geometric Product)
2. 外积 / 楔积 (Outer/Wedge Product)
3. 内积 (Inner Product)

作者: nblade 团队
"""

import nblade
import math


def section_1_geometric_product():
    """第一节：几何积"""
    print("\n" + "=" * 60)
    print("第一节：几何积 (Geometric Product)")
    print("=" * 60)

    alg = nblade.Algebra.euclidean(3)
    e1, e2, e3 = alg.basis_vectors()

    print("\n几何积是几何代数的基础运算，定义为:")
    print("  ab = a·b + a∧b")
    print("其中 a·b 是内积（对称部分），a∧b 是外积（反对称部分）")

    # 正交向量的几何积
    print("\n【正交向量的几何积】")
    prod_orth = e1 * e2
    print(f"e1 * e2 = {prod_orth}")
    print("  正交向量的内积为 0，所以几何积等于外积")

    # 平行向量的几何积
    print("\n【平行向量的几何积】")
    prod_para = e1 * e1
    print(f"e1 * e1 = {prod_para}")
    print("  平行向量的外积为 0，所以几何积等于内积（= 1，因为 e1² = 1）")

    # 一般向量的几何积
    print("\n【一般向量的几何积】")
    a = alg.vector([1.0, 2.0, 0.0])
    b = alg.vector([2.0, 1.0, 0.0])

    geom = a * b
    print(f"a = {a}")
    print(f"b = {b}")
    print(f"a * b = {geom}")

    # 分解验证
    inner = a | b
    outer = a ^ b
    print(f"\n分解验证:")
    print(f"  a·b (内积) = {inner}")
    print(f"  a∧b (外积) = {outer}")
    print(f"  a·b + a∧b = {inner + outer}")
    print(f"  与 a*b 比较: {(geom - (inner + outer)).norm():.2e}")

    # 几何积的性质
    print("\n【几何积的性质】")
    print(f"结合律: (e1*e2)*e3 = e1*(e2*e3)")
    left = (e1 * e2) * e3
    right = e1 * (e2 * e3)
    print(f"  左边 = {left}")
    print(f"  右边 = {right}")
    print(f"  相等: {(left - right).norm() < 1e-10}")

    print(f"\n非交换性: e1*e2 ≠ e2*e1")
    print(f"  e1*e2 = {e1 * e2}")
    print(f"  e2*e1 = {e2 * e1}")
    print(f"  (注意符号相反)")


def section_2_outer_product():
    """第二节：外积（楔积）"""
    print("\n" + "=" * 60)
    print("第二节：外积 / 楔积 (Outer/Wedge Product)")
    print("=" * 60)

    alg = nblade.Algebra.euclidean(3)
    e1, e2, e3 = alg.basis_vectors()

    print("\n外积表示张成的子空间，具有以下性质:")
    print("  - 反对称性: a∧b = -b∧a")
    print("  - 结合性: (a∧b)∧c = a∧(b∧c)")
    print("  - 线性性: (αa + βb)∧c = α(a∧c) + β(b∧c)")

    # 向量的外积 → 二重向量
    print("\n【向量的外积 → 二重向量 (平面)】")
    bivector = e1 ^ e2
    print(f"e1 ^ e2 = {bivector}")
    print("  表示 e1 和 e2 张成的平面 (xy 平面)")

    # 反对称性
    print("\n【反对称性】")
    print(f"e1 ^ e2 = {e1 ^ e2}")
    print(f"e2 ^ e1 = {e2 ^ e1}")
    print(f"两者相差负号")

    # 向量自身的楔积为零
    print("\n【向量与自身的外积为零】")
    print(f"e1 ^ e1 = {e1 ^ e1}")
    print("向量与自身线性相关，张成的面积为 0")

    # 三个向量的外积 → 三重向量
    print("\n【三个向量的外积 → 三重向量 (体积)】")
    trivector = e1 ^ e2 ^ e3
    print(f"e1 ^ e2 ^ e3 = {trivector}")
    print("  表示三个向量张成的体积 (伪标量)")

    # 计算平行四边形面积
    print("\n【应用：计算平行四边形面积】")
    a = alg.vector([2.0, 0.0, 0.0])
    b = alg.vector([1.0, 1.0, 0.0])

    area_bivector = a ^ b
    area = area_bivector.norm()
    print(f"a = {a}")
    print(f"b = {b}")
    print(f"a ∧ b = {area_bivector}")
    print(f"面积 |a∧b| = {area:.6f} (应为 2)")

    # 高维外积
    print("\n【高维中的外积】")
    alg_4d = nblade.Algebra.euclidean(4)
    e1_4d, e2_4d, e3_4d, e4_4d = alg_4d.basis_vectors()

    blade_3 = e1_4d ^ e2_4d ^ e3_4d
    print(f"4D 中的三重向量: e1∧e2∧e3 = {blade_3}")

    blade_4 = e1_4d ^ e2_4d ^ e3_4d ^ e4_4d
    print(f"4D 中的四重向量 (伪标量): e1∧e2∧e3∧e4 = {blade_4}")


def section_3_inner_product():
    """第三节：内积"""
    print("\n" + "=" * 60)
    print("第三节：内积 (Inner Product)")
    print("=" * 60)

    alg = nblade.Algebra.euclidean(3)
    e1, e2, e3 = alg.basis_vectors()

    print("\n内积是几何代数中的收缩运算，nblade 支持:")
    print("  - 左内积 (左收缩): a | b 或 a.left_inner(b)")
    print("  - 右内积 (右收缩): a.right_inner(b)")

    # 向量的内积
    print("\n【向量的内积 → 标量】")
    a = alg.vector([1.0, 2.0, 3.0])
    b = alg.vector([4.0, 5.0, 6.0])

    inner = a | b
    print(f"a = {a}")
    print(f"b = {b}")
    print(f"a | b = {inner}")
    print("  结果是标量，与传统的点积相同")

    # 正交性检验
    print("\n【正交性检验】")
    print(f"e1 | e2 = {e1 | e2} (正交，结果为 0)")
    print(f"e1 | e1 = {e1 | e1} (非正交，结果为 1)")

    # 计算夹角
    print("\n【计算向量夹角】")
    v1 = alg.vector([1.0, 0.0, 0.0])
    v2 = alg.vector([1.0, 1.0, 0.0])

    # cos(θ) = (a·b) / (|a||b|)
    dot = (v1 | v2).scalar_part()
    norm1 = v1.norm()
    norm2 = v2.norm()
    cos_theta = dot / (norm1 * norm2)
    theta = math.acos(cos_theta)

    print(f"v1 = {v1}")
    print(f"v2 = {v2}")
    print(f"内积 v1·v2 = {dot}")
    print(f"cos(θ) = {cos_theta:.6f}")
    print(f"夹角 θ = {math.degrees(theta):.2f}°")

    # 向量与二重向量的内积
    print("\n【向量与二重向量的内积 → 向量】")
    B = e1 ^ e2  # xy 平面
    v = e1

    result = v | B
    print(f"二重向量 B = {B}")
    print(f"向量 v = {v}")
    print(f"v | B = {result}")
    print("  向量与二重向量的内积降阶，结果是向量")


def section_4_relationships():
    """第四节：三种乘积的关系"""
    print("\n" + "=" * 60)
    print("第四节：三种乘积的关系")
    print("=" * 60)

    alg = nblade.Algebra.euclidean(3)
    e1, e2, e3 = alg.basis_vectors()

    print("\n核心公式:")
    print("  几何积 = 内积 + 外积")
    print("  ab = a·b + a∧b")

    a = alg.vector([1.0, 2.0, 0.0])
    b = alg.vector([3.0, 1.0, 0.0])

    print(f"\na = {a}")
    print(f"b = {b}")

    geom = a * b
    inner = a | b
    outer = a ^ b

    print(f"\n几何积 a*b = {geom}")
    print(f"内积 a|b   = {inner}")
    print(f"外积 a^b   = {outer}")

    # 验证关系
    sum_parts = inner + outer
    print(f"\n内积 + 外积 = {sum_parts}")

    diff = geom - sum_parts
    print(f"与几何积的差 = {diff}")
    print(f"差值范数 = {diff.norm():.2e} (应接近 0)")

    # 对称性分析
    print("\n【对称性分析】")
    print(f"内积 a|b = b|a: {(a | b) - (b | a)}")
    print("  内积是对称的")
    print(f"\n外积 a^b = -(b^a): {(a ^ b) + (b ^ a)}")
    print("  外积是反对称的")
    print(f"\n几何积 ab ≠ ba:")
    print(f"  a*b = {a * b}")
    print(f"  b*a = {b * a}")
    print("  几何积既不对称也不反对称")


def section_5_practical_examples():
    """第五节：实际应用示例"""
    print("\n" + "=" * 60)
    print("第五节：实际应用示例")
    print("=" * 60)

    alg = nblade.Algebra.euclidean(3)
    e1, e2, e3 = alg.basis_vectors()

    # 应用1：投影
    print("\n【应用1：向量投影】")
    a = alg.vector([1.0, 2.0, 0.0])
    b = alg.vector([1.0, 0.0, 0.0])

    # 投影公式: proj_b(a) = (a·b / b·b) * b
    dot_ab = (a | b).scalar_part()
    dot_bb = (b | b).scalar_part()
    proj_a_on_b = b.scale(dot_ab / dot_bb)

    print(f"向量 a = {a}")
    print(f"向量 b = {b}")
    print(f"a 在 b 上的投影 = {proj_a_on_b}")

    # 应用2：计算距离
    print("\n【应用2：点到原点的距离】")
    p = alg.vector([3.0, 4.0, 0.0])
    distance = p.norm()
    print(f"点 p = {p}")
    print(f"到原点的距离 |p| = {distance:.6f} (应为 5)")

    # 应用3：计算面积
    print("\n【应用3：三角形面积】")
    v1 = alg.vector([1.0, 0.0, 0.0])
    v2 = alg.vector([0.0, 1.0, 0.0])

    area_bivector = v1 ^ v2
    area = 0.5 * area_bivector.norm()
    print(f"三角形的边: {v1}, {v2}")
    print(f"三角形面积 = 0.5 * |v1∧v2| = {area:.6f}")

    # 应用4：判断共面
    print("\n【应用4：判断三个向量是否共面】")
    v1 = alg.vector([1.0, 0.0, 0.0])
    v2 = alg.vector([0.0, 1.0, 0.0])
    v3 = alg.vector([1.0, 1.0, 0.0])  # 与 v1, v2 共面

    triple = v1 ^ v2 ^ v3
    print(f"v1 = {v1}")
    print(f"v2 = {v2}")
    print(f"v3 = {v3}")
    print(f"v1∧v2∧v3 = {triple}")
    print(f"范数 = {triple.norm():.2e} (接近 0 表示共面)")


def main():
    """主函数"""
    print("=" * 60)
    print("nblade 基本运算教程")
    print("三种乘积：几何积、外积、内积")
    print("=" * 60)

    section_1_geometric_product()
    section_2_outer_product()
    section_3_inner_product()
    section_4_relationships()
    section_5_practical_examples()

    print("\n" + "=" * 60)
    print("教程完成！")
    print("=" * 60)
    print("\n总结:")
    print("  几何积 (*): 最基本的运算，包含完整信息")
    print("  外积 (^): 表示张成的子空间，反对称")
    print("  内积 (|): 表示投影/收缩，对称")
    print("\n核心公式: ab = a·b + a∧b")


if __name__ == "__main__":
    main()
