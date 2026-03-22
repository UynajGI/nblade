"""
nblade 转子插值示例
==================

演示几何代数中的转子插值技术：
- 线性插值 (LERP) 及其局限性
- 球面线性插值 (SLERP) - 正确方法
- 转子 SLERP 的实现
- 平滑动画应用
- 与四元数 SLERP 的对比

转子插值是动画、机器人学和计算机图形学的核心技术。

作者: nblade 团队
"""

import nblade
import math


def section_1_rotor_review():
    """第一节：转子快速回顾"""
    print("\n" + "=" * 60)
    print("第一节：转子快速回顾")
    print("=" * 60)

    alg = nblade.Algebra.euclidean(3)
    e1, e2, e3 = alg.basis_vectors()

    print("\n转子 R 表示旋转，具有以下性质：")
    print("  1. R 是偶元素（标量 + 二重向量）")
    print("  2. R·R† = 1（归一化条件）")
    print("  3. 旋转公式: v' = R·v·R†")

    # 创建转子
    plane = e1 ^ e2
    angle = math.pi / 4
    R = alg.rotor(plane, angle)

    print(f"\n示例转子:")
    print(f"  旋转平面: {plane} (xy 平面)")
    print(f"  旋转角度: {math.degrees(angle)}°")
    print(f"  转子 R = {R}")

    # 验证归一化
    R_rev = R.reversion()
    RR = R * R_rev
    print(f"\n归一化验证:")
    print(f"  R·R† = {RR.scalar_part():.10f} (应等于 1)")


def section_2_lerp_limitations():
    """第二节：线性插值 (LERP) 及其局限性"""
    print("\n" + "=" * 60)
    print("第二节：线性插值 (LERP) 及其局限性")
    print("=" * 60)

    alg = nblade.Algebra.euclidean(3)
    e1, e2, e3 = alg.basis_vectors()

    print("\n线性插值 (LERP) 是最简单的插值方法：")
    print("  R(t) = (1-t)·R1 + t·R2")
    print("\n但 LERP 对转子有两个问题：")
    print("  1. 插值结果不保持归一化")
    print("  2. 角速度不恒定（非等速旋转）")

    # 两个转子
    R1 = alg.rotor(e1 ^ e2, 0)  # 初始姿态
    R2 = alg.rotor(e1 ^ e2, math.pi / 2)  # 最终姿态 (90°)

    print(f"\n【LERP 演示】")
    print(f"R1 (初始): {R1}")
    print(f"R2 (最终): {R2}")

    print("\n插值过程:")
    print("  t    |R(t)|²    问题")
    print("  " + "-" * 40)

    for i in range(6):
        t = i / 5

        # LERP
        R_lerp = R1.scale(1 - t) + R2.scale(t)

        # 检查归一化
        norm_sq = (R_lerp * R_lerp.reversion()).scalar_part()

        # 角度（从标量部分计算）
        scalar_part = R_lerp.scalar_part()
        angle = 2 * math.acos(
            min(1, max(-1, abs(scalar_part) / max(norm_sq**0.5, 1e-10)))
        )

        print(
            f"  {t:.1f}  {norm_sq:.6f}   {'归一化' if abs(norm_sq - 1) < 0.01 else '未归一化!'}"
        )

    print("\n结论: LERP 破坏了转子的归一化，需要重新归一化")
    print("      即使重新归一化，角速度仍然不恒定")


def section_3_slerp_theory():
    """第三节：球面线性插值 (SLERP) 理论"""
    print("\n" + "=" * 60)
    print("第三节：球面线性插值 (SLERP) 理论")
    print("=" * 60)

    print("\nSLERP 在单位四维球面上进行插值，保持等速旋转")
    print("\n公式:")
    print("  R(t) = sin((1-t)θ)/sin(θ) · R1 + sin(tθ)/sin(θ) · R2")
    print("\n其中 θ 是两个转子之间的角度")

    print("\n【计算转子间角度】")
    print("对于归一化转子 R1 和 R2:")
    print("  cos(θ/2) = <R1, R2> = (R1†·R2)₀ 的标量部分")

    alg = nblade.Algebra.euclidean(3)
    e1, e2, e3 = alg.basis_vectors()

    R1 = alg.rotor(e1 ^ e2, 0)
    R2 = alg.rotor(e1 ^ e2, math.pi / 2)

    # 计算角度
    # R1†·R2 的标量部分
    R1_rev = R1.reversion()
    R1R2 = R1_rev * R2
    cos_half_theta = R1R2.scalar_part()

    print(f"\nR1 = {R1}")
    print(f"R2 = {R2}")
    print(f"R1†·R2 = {R1R2}")
    print(f"cos(θ/2) = {cos_half_theta}")

    theta = 2 * math.acos(cos_half_theta)
    print(f"\n两个转子间的角度 θ = {math.degrees(theta):.2f}°")


def slerp_rotor(R1, R2, t, alg):
    """
    转子的球面线性插值 (SLERP)

    参数:
        R1, R2: 起始和结束转子（归一化）
        t: 插值参数 [0, 1]
        alg: 代数对象

    返回:
        插值后的归一化转子
    """
    # 计算点积 <R1, R2>
    R1_rev = R1.reversion()
    dot = (R1_rev * R2).scalar_part()

    # 处理数值稳定性
    dot = min(1.0, max(-1.0, dot))

    # 计算角度
    theta = math.acos(abs(dot))

    # 处理特殊情况
    if abs(theta) < 1e-10:
        # 角度太小，使用 LERP
        return R1.scale(1 - t) + R2.scale(t)

    # SLERP 公式
    sin_theta = math.sin(theta)
    coeff1 = math.sin((1 - t) * theta) / sin_theta
    coeff2 = math.sin(t * theta) / sin_theta

    # 如果点积为负，需要反转一个转子
    if dot < 0:
        coeff1 = -coeff1

    R_interp = R1.scale(coeff1) + R2.scale(coeff2)

    return R_interp


def section_4_slerp_implementation():
    """第四节：SLERP 实现"""
    print("\n" + "=" * 60)
    print("第四节：SLERP 实现与验证")
    print("=" * 60)

    alg = nblade.Algebra.euclidean(3)
    e1, e2, e3 = alg.basis_vectors()

    print("\n【SLERP 插值演示】")

    # 两个转子：从 0° 到 90°
    R1 = alg.rotor(e1 ^ e2, 0)
    R2 = alg.rotor(e1 ^ e2, math.pi / 2)

    print(f"R1 (初始): 旋转 0°")
    print(f"R2 (最终): 旋转 90°")

    print("\n插值过程:")
    print("  t    角度      |R(t)|²    验证")
    print("  " + "-" * 50)

    for i in range(11):
        t = i / 10

        # SLERP
        R_slerp = slerp_rotor(R1, R2, t, alg)

        # 检查归一化
        norm_sq = (R_slerp * R_slerp.reversion()).scalar_part()

        # 计算实际旋转角度
        scalar_part = R_slerp.scalar_part()
        angle = 2 * math.acos(min(1, max(-1, abs(scalar_part))))

        # 理论角度
        expected_angle = t * math.pi / 2

        print(
            f"  {t:.1f}  {math.degrees(angle):6.2f}°   {norm_sq:.10f}   {'✓' if abs(norm_sq - 1) < 1e-6 else '✗'}"
        )

    print("\n结论: SLERP 保持了归一化，角度线性增长")


def section_5_smooth_animation():
    """第五节：平滑动画应用"""
    print("\n" + "=" * 60)
    print("第五节：平滑动画应用")
    print("=" * 60)

    alg = nblade.Algebra.euclidean(3)
    e1, e2, e3 = alg.basis_vectors()

    print("\n模拟向量从 v1 平滑旋转到 v2")

    # 起始和结束向量
    v1 = e1
    v2 = e2

    print(f"起始向量: v1 = {v1}")
    print(f"结束向量: v2 = {v2}")

    # 计算从 v1 到 v2 的转子
    # 方法：通过旋转平面和角度
    rotation_plane = v1 ^ v2

    # 角度
    cos_angle = (v1 | v2).scalar_part() / (v1.norm() * v2.norm())
    angle = math.acos(cos_angle)

    # 转子
    R1 = alg.rotor(rotation_plane, 0)
    R2 = alg.rotor(rotation_plane, angle)

    print(f"\n旋转平面: {rotation_plane}")
    print(f"旋转角度: {math.degrees(angle):.2f}°")

    # 动画帧
    print("\n【动画帧 (10 帧)】")
    positions = []
    for i in range(11):
        t = i / 10
        R = slerp_rotor(R1, R2, t, alg)
        v = v1.rotate_by(R)
        positions.append(v)
        print(f"  帧 {i:2d}: v = {v}")

    # 验证端点
    print(f"\n验证:")
    print(f"  起点: {positions[0]} (应接近 {v1})")
    print(f"  终点: {positions[-1]} (应接近 {v2})")


def section_6_quaternion_comparison():
    """第六节：与四元数 SLERP 对比"""
    print("\n" + "=" * 60)
    print("第六节：转子 SLERP 与四元数 SLERP 对比")
    print("=" * 60)

    print("\n四元数和转子本质上是相同的：")
    print("  四元数 q = w + xi + yj + zk")
    print("  转子   R = a + b·e2∧e3 + c·e3∧e1 + d·e1∧e2")
    print("\n对应关系:")
    print("  i ↔ e2∧e3,  j ↔ e3∧e1,  k ↔ e1∧e2")
    print("  w ↔ 标量部分")

    alg = nblade.Algebra.euclidean(3)
    e1, e2, e3 = alg.basis_vectors()

    print("\n【四元数 SLERP 公式】")
    print("  q(t) = sin((1-t)θ)/sin(θ) · q1 + sin(tθ)/sin(θ) · q2")
    print("\n这与转子 SLERP 公式完全相同！")

    print("\n【示例：绕 z 轴旋转】")

    # 创建绕 z 轴旋转 120° 的转子
    R = alg.rotor(e1 ^ e2, 2 * math.pi / 3)

    # 提取四元数分量
    # R = cos(θ/2) + sin(θ/2)·(e1∧e2)
    w = R.scalar_part()

    # 二重向量部分
    bivector_part = R - alg.scalar(w)
    # 对于绕 z 轴旋转，只有 e1∧e2 分量
    k_coeff = (bivector_part | (e1 ^ e2)).scalar_part() / (e1 ^ e2).norm_squared()

    print(f"\n转子 R = {R}")
    print(f"四元数形式: q = {w:.4f} + {k_coeff:.4f}k")
    print(f"  其中 k 对应 e1∧e2")

    # 归一化验证
    norm_sq = (R * R.reversion()).scalar_part()
    print(f"\n归一化验证: |R|² = {norm_sq:.10f}")

    print("\n【优势对比】")
    print("┌─────────────────┬──────────────────┐")
    print("│ 转子 SLERP      │ 四元数 SLERP     │")
    print("├─────────────────┼──────────────────┤")
    print("│ 几何意义明确    │ 需要抽象理解     │")
    print("│ 任意维度通用    │ 仅限 3D          │")
    print("│ 与其他 GA 运算  │ 需要单独实现     │")
    print("│ 统一框架        │                  │")
    print("└─────────────────┴──────────────────┘")


def section_7_practical_tips():
    """第七节：实用技巧"""
    print("\n" + "=" * 60)
    print("第七节：实用技巧与注意事项")
    print("=" * 60)

    alg = nblade.Algebra.euclidean(3)
    e1, e2, e3 = alg.basis_vectors()

    print("\n【技巧 1：处理接近的转子】")
    print("当两个转子非常接近时，sin(θ) ≈ 0，需要特殊处理")
    print("解决方案：当 θ < ε 时使用 LERP")

    R1 = alg.rotor(e1 ^ e2, 0.001)  # 非常小的旋转
    R2 = alg.rotor(e1 ^ e2, 0.002)

    R_interp = slerp_rotor(R1, R2, 0.5, alg)
    print(f"\n小角度插值成功: {R_interp}")

    print("\n【技巧 2：选择最短路径】")
    print("转子的双重覆盖性质：R 和 -R 表示相同旋转")
    print("选择点积为正的方向以确保最短路径")

    R1 = alg.rotor(e1 ^ e2, 0)
    R2 = alg.rotor(e1 ^ e2, math.pi * 1.5)  # 270°

    # 检查点积
    dot = (R1.reversion() * R2).scalar_part()
    print(f"\n点积 cos(θ/2) = {dot:.4f}")

    if dot < 0:
        print("点积为负，应使用 -R2 以获得最短路径")
        R2_alt = R2.scale(-1)
        dot_alt = (R1.reversion() * R2_alt).scalar_part()
        print(f"使用 -R2 后的点积 = {dot_alt:.4f}")

    print("\n【技巧 3：恒定角速度】")
    print("SLERP 保证恒定角速度，适合动画关键帧")
    print("对于 n 帧，每帧角度变化 = θ/n")


def main():
    """主函数"""
    print("=" * 60)
    print("nblade 转子插值示例")
    print("SLERP: 球面线性插值的正确方法")
    print("=" * 60)

    section_1_rotor_review()
    section_2_lerp_limitations()
    section_3_slerp_theory()
    section_4_slerp_implementation()
    section_5_smooth_animation()
    section_6_quaternion_comparison()
    section_7_practical_tips()

    print("\n" + "=" * 60)
    print("示例完成！")
    print("=" * 60)
    print("\n关键点:")
    print("  - LERP 破坏归一化，不适用于转子")
    print("  - SLERP 在单位球面上插值，保持等速")
    print("  - SLERP 公式: R(t) = sin((1-t)θ)/sin(θ)·R1 + sin(tθ)/sin(θ)·R2")
    print("  - 转子 SLERP 与四元数 SLERP 本质相同")
    print("  - 注意处理小角度和最短路径")


if __name__ == "__main__":
    main()
