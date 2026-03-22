"""
nblade 旋转教程
===============

详细介绍如何使用几何代数进行旋转：
- 2D 旋转与复数
- 3D 旋转与四元数/转子
- 旋转的组合与插值
- 实际应用示例

作者: nblade 团队
"""

import nblade
import math


def section_1_2d_rotations():
    """第一节：2D 旋转与复数"""
    print("\n" + "=" * 60)
    print("第一节：2D 旋转与复数")
    print("=" * 60)

    alg = nblade.Algebra.euclidean(2)
    e1, e2 = alg.basis_vectors()

    print("\n在 2D 几何代数中，伪标量 I = e1∧e2 具有性质 I² = -1")
    print("这与复数中的虚数单位 i 相同")

    # 伪标量
    I = e1 ^ e2
    print(f"\n伪标量 I = e1∧e2 = {I}")
    print(f"I² = {(I * I).scalar_part()} (等于 -1)")

    # 复数表示
    print("\n【复数表示】")
    # z = a + bI 类似于复数 a + bi
    z1 = alg.scalar(3.0) + I.scale(2.0)  # 3 + 2I
    print(f"复数 z = 3 + 2I = {z1}")

    # 旋转算子
    print("\n【旋转算子】")
    angle = math.pi / 4  # 45 度
    cos_a = math.cos(angle)
    sin_a = math.sin(angle)

    # e^(Iθ) = cos(θ) + I·sin(θ)
    rotor_2d = alg.scalar(cos_a) + I.scale(sin_a)
    print(f"旋转算子 R = cos(45°) + I·sin(45°)")
    print(f"           = {rotor_2d}")

    # 旋转向量
    print("\n【旋转向量】")
    v = e1
    print(f"原始向量: v = {v}")

    # 2D 旋转: v' = R·v (单侧乘法)
    rotated = rotor_2d * v * rotor_2d.reversion()  # 或简化为 rotor * v (对于2D偶元素)
    print(f"旋转 45° 后: v' = {rotated}")

    # 理论验证
    expected = e1.scale(cos_a) + e2.scale(sin_a)
    print(f"理论预期: {expected}")
    print(f"误差: {(rotated - expected).norm():.2e}")


def section_2_3d_rotations_rotors():
    """第二节：3D 旋转与转子"""
    print("\n" + "=" * 60)
    print("第二节：3D 旋转与转子")
    print("=" * 60)

    alg = nblade.Algebra.euclidean(3)
    e1, e2, e3 = alg.basis_vectors()

    print("\n3D 旋转使用转子 (Rotor)，定义在旋转平面上")
    print("转子: R = exp(-B·θ/2)，其中 B 是单位二重向量（旋转平面）")

    # 创建转子
    print("\n【创建转子】")
    plane = e1 ^ e2  # xy 平面
    angle = math.pi / 3  # 60 度

    rotor = alg.rotor(plane, angle)
    print(f"旋转平面: {plane} (xy 平面)")
    print(f"旋转角度: {math.degrees(angle)}°")
    print(f"转子 R = {rotor}")

    # 转子性质
    print("\n【转子性质】")
    rotor_rev = rotor.reversion()
    print(f"R† (反转) = {rotor_rev}")

    # R·R† 应该是 1 (归一化)
    rr = rotor * rotor_rev
    print(f"R·R† = {rr.scalar_part():.6f} (应接近 1)")

    # 旋转向量
    print("\n【旋转向量】")
    v = e1
    print(f"原始向量: v = {v}")

    # 旋转公式: v' = R·v·R†
    rotated = v.rotate_by(rotor)
    print(f"旋转后: v' = {rotated}")

    # 手动计算验证
    rotated_manual = rotor * v * rotor.reversion()
    print(f"手动计算: R·v·R† = {rotated_manual}")

    # 理论值
    cos_a = math.cos(angle)
    sin_a = math.sin(angle)
    expected = e1.scale(cos_a) + e2.scale(sin_a)
    print(f"理论预期: cos(60°)e1 + sin(60°)e2 = {expected}")

    # 旋转任意向量
    print("\n【旋转任意向量】")
    v_arb = alg.vector([1.0, 2.0, 3.0])
    rotated_arb = v_arb.rotate_by(rotor)
    print(f"v = {v_arb}")
    print(f"旋转 60° 后 = {rotated_arb}")


def section_3_multiple_rotations():
    """第三节：多个旋转的组合"""
    print("\n" + "=" * 60)
    print("第三节：多个旋转的组合")
    print("=" * 60)

    alg = nblade.Algebra.euclidean(3)
    e1, e2, e3 = alg.basis_vectors()

    print("\n多个旋转的组合通过转子的几何积实现")
    print("R_combined = R2 · R1 (先应用 R1，再应用 R2)")

    # 第一个旋转: 绕 z 轴 (xy 平面) 旋转 45°
    r1_plane = e1 ^ e2
    r1_angle = math.pi / 4
    R1 = alg.rotor(r1_plane, r1_angle)

    # 第二个旋转: 绕 x 轴 (yz 平面) 旋转 30°
    r2_plane = e2 ^ e3
    r2_angle = math.pi / 6
    R2 = alg.rotor(r2_plane, r2_angle)

    print(f"\nR1: xy 平面旋转 45°")
    print(f"R2: yz 平面旋转 30°")

    # 组合旋转
    R_combined = R2 * R1
    print(f"\n组合转子 R_combined = R2·R1 = {R_combined}")

    # 应用组合旋转
    v = e1
    print(f"\n原始向量: v = {v}")

    # 方法1: 依次旋转
    v_rot1 = v.rotate_by(R1)
    v_rot2 = v_rot1.rotate_by(R2)
    print(f"依次旋转: R2(R1(v)) = {v_rot2}")

    # 方法2: 组合旋转
    v_combined = v.rotate_by(R_combined)
    print(f"组合旋转: R_combined(v) = {v_combined}")

    # 验证
    error = (v_rot2 - v_combined).norm()
    print(f"\n误差: {error:.2e}")


def section_4_rotation_interpolation():
    """第四节：旋转插值 (球面线性插值)"""
    print("\n" + "=" * 60)
    print("第四节：旋转插值 (SLERP)")
    print("=" * 60)

    alg = nblade.Algebra.euclidean(3)
    e1, e2, e3 = alg.basis_vectors()

    print("\n转子可以用于平滑插值两个旋转")

    # 两个旋转
    R1 = alg.rotor(e1 ^ e2, 0)  # 初始姿态
    R2 = alg.rotor(e1 ^ e2, math.pi / 2)  # 最终姿态 (旋转 90°)

    print(f"R1 (初始): {R1}")
    print(f"R2 (最终): {R2}")

    # 简单的线性插值 (对于归一化转子)
    print("\n【插值演示】")
    num_steps = 5
    for i in range(num_steps + 1):
        t = i / num_steps

        # 简化插值: R(t) = (1-t)R1 + tR2 (需要重新归一化)
        # 更精确的方法是使用 SLERP
        R_interp = R1.scale(1 - t) + R2.scale(t)

        # 归一化 (如果需要)
        # norm = (R_interp * R_interp.reversion()).scalar_part() ** 0.5
        # R_interp = R_interp.scale(1 / norm)

        v_rotated = e1.rotate_by(R_interp)
        print(f"t = {t:.1f}: {v_rotated}")


def section_5_quaternions():
    """第五节：转子与四元数的关系"""
    print("\n" + "=" * 60)
    print("第五节：转子与四元数")
    print("=" * 60)

    alg = nblade.Algebra.euclidean(3)
    e1, e2, e3 = alg.basis_vectors()

    print("\n四元数可以自然地嵌入到 3D 几何代数中")
    print("四元数单位 i, j, k 对应二重向量:")

    # 四元数单位
    i = e2 ^ e3
    j = e3 ^ e1
    k = e1 ^ e2

    print(f"  i = e2∧e3 = {i}")
    print(f"  j = e3∧e1 = {j}")
    print(f"  k = e1∧e2 = {k}")

    # 验证四元数乘法规则
    print("\n【四元数乘法规则验证】")
    print(f"i² = {(i * i).scalar_part()} (应为 -1)")
    print(f"j² = {(j * j).scalar_part()} (应为 -1)")
    print(f"k² = {(k * k).scalar_part()} (应为 -1)")

    print(f"\nij = {i * j}")
    print(f"jk = {j * k}")
    print(f"ki = {k * i}")

    ijk = i * j * k
    print(f"\nijk = {ijk} (应为 -1)")

    # 四元数旋转
    print("\n【四元数旋转示例】")
    # 四元数 q = cos(θ/2) + sin(θ/2)·(ax·i + ay·j + az·k)
    # 对应转子 R = cos(θ/2) + sin(θ/2)·(ax·e2∧e3 + ay·e3∧e1 + az·e1∧e2)

    # 绕 z 轴旋转 60°
    angle = math.pi / 3
    axis = alg.vector([0, 0, 1])  # z 轴

    # 使用 nblade 的转子创建
    plane = e1 ^ e2  # xy 平面 (垂直于 z 轴)
    rotor = alg.rotor(plane, angle)

    print(f"绕 z 轴旋转 60°")
    print(f"转子: {rotor}")

    v = alg.vector([1, 0, 0])
    rotated = v.rotate_by(rotor)
    print(f"\n{v} 旋转后 = {rotated}")


def section_6_practical_examples():
    """第六节：实际应用示例"""
    print("\n" + "=" * 60)
    print("第六节：实际应用示例")
    print("=" * 60)

    alg = nblade.Algebra.euclidean(3)
    e1, e2, e3 = alg.basis_vectors()

    # 示例1: 旋转一个三角形
    print("\n【示例1: 旋转三角形】")
    p1 = e1
    p2 = e2
    p3 = e3 * 0.5 + e1

    print("原始三角形顶点:")
    print(f"  p1 = {p1}")
    print(f"  p2 = {p2}")
    print(f"  p3 = {p3}")

    # 绕 z 轴旋转 90°
    rotor = alg.rotor(e1 ^ e2, math.pi / 2)

    p1_rot = p1.rotate_by(rotor)
    p2_rot = p2.rotate_by(rotor)
    p3_rot = p3.rotate_by(rotor)

    print("\n旋转 90° 后:")
    print(f"  p1' = {p1_rot}")
    print(f"  p2' = {p2_rot}")
    print(f"  p3' = {p3_rot}")

    # 示例2: 计算旋转轴
    print("\n【示例2: 从两个向量确定旋转】")
    v1 = alg.vector([1, 0, 0])
    v2 = alg.vector([0, 1, 0])

    print(f"向量 v1 = {v1}")
    print(f"向量 v2 = {v2}")

    # 计算旋转平面
    # 从 v1 到 v2 的旋转平面是 v1∧v2
    rotation_plane = v1 ^ v2
    print(f"\n旋转平面 v1∧v2 = {rotation_plane}")

    # 旋转角度
    cos_angle = (v1 | v2).scalar_part() / (v1.norm() * v2.norm())
    angle = math.acos(cos_angle)
    print(f"旋转角度: {math.degrees(angle):.2f}°")


def main():
    """主函数"""
    print("=" * 60)
    print("nblade 旋转教程")
    print("从 2D 到 3D，理解几何代数中的旋转")
    print("=" * 60)

    section_1_2d_rotations()
    section_2_3d_rotations_rotors()
    section_3_multiple_rotations()
    section_4_rotation_interpolation()
    section_5_quaternions()
    section_6_practical_examples()

    print("\n" + "=" * 60)
    print("教程完成！")
    print("=" * 60)
    print("\n关键点:")
    print("  - 2D 旋转: 使用伪标量 I = e1∧e2")
    print("  - 3D 旋转: 使用转子 R = exp(-B·θ/2)")
    print("  - 旋转公式: v' = R·v·R†")
    print("  - 多个旋转: R_combined = R2·R1")
    print("  - 四元数: i,j,k 对应 e2∧e3, e3∧e1, e1∧e2")


if __name__ == "__main__":
    main()
