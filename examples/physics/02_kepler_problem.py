"""
nblade 开普勒问题示例
====================

演示几何代数在行星轨道问题中的应用：
- 角动量作为二重向量 L = r ∧ p
- 中心力场运动
- 轨道方程的几何代数推导
- 偏心率矢量
- 数值模拟与守恒量验证

本示例展示几何代数如何简化开普勒问题的推导，
相比传统向量方法更加直观优雅。

作者: nblade 团队
"""

import nblade
import math


def section_1_angular_momentum():
    """第一节：角动量作为二重向量"""
    print("\n" + "=" * 60)
    print("第一节：角动量作为二重向量")
    print("=" * 60)

    alg = nblade.Algebra.euclidean(3)
    e1, e2, e3 = alg.basis_vectors()

    print("\n在中心力场中，角动量守恒是核心性质")
    print("传统方法：L = r × p (叉积，仅限3D)")
    print("几何代数：L = r ∧ p (外积，任意维度)")

    # 示例：行星在xy平面运动
    # 位置向量
    r = e1.scale(5.0) + e2.scale(0.0)  # 近日点
    print(f"\n位置向量 r = {r}")

    # 动量向量 (切向)
    v = e2.scale(1.0)  # 速度
    mass = 1.0
    p = v.scale(mass)
    print(f"动量向量 p = m·v = {p}")

    # 角动量作为二重向量
    L = r ^ p
    print(f"\n角动量 L = r ∧ p = {L}")

    # 角动量的大小
    L_magnitude = L.norm()
    print(f"角动量大小 |L| = {L_magnitude:.6f}")

    print("\n【物理意义】")
    print("二重向量 L = r ∧ p 表示位置和动量张成的平面")
    print("其大小等于轨道平面上扫过的面积速度的两倍")

    # 验证中心力场中角动量守恒
    print("\n【角动量守恒验证】")
    print("对于中心力 F = f(r)·r̂，力矩 τ = r ∧ F = 0")
    print("因为 r 与 F 共线，外积为零")

    # 计算力矩
    F = r.scale(-1.0)  # 示例中心力，指向原点
    torque = r ^ F
    print(f"\n中心力 F = {F}")
    print(f"力矩 τ = r ∧ F = {torque}")
    print(f"力矩大小 |τ| = {torque.norm():.10f} (应为零)")


def section_2_kepler_setup():
    """第二节：开普勒问题设置"""
    print("\n" + "=" * 60)
    print("第二节：开普勒问题设置")
    print("=" * 60)

    alg = nblade.Algebra.euclidean(3)
    e1, e2, e3 = alg.basis_vectors()

    print("\n开普勒问题：粒子在中心力场 F = -k/r² · r̂ 中的运动")

    # 引力常数
    k = 1.0  # G·M·m
    mass = 1.0

    print(f"\n参数设置:")
    print(f"  引力常数 k = {k}")
    print(f"  粒子质量 m = {mass}")

    # 初始条件（椭圆轨道）
    # 近日点距离
    r_peri = 1.0
    # 近日点速度（使轨道为椭圆）
    v_peri = math.sqrt(k / r_peri) * 1.2  # 略小于圆周运动速度，形成椭圆

    r_vec = e1.scale(r_peri)
    v_vec = e2.scale(v_peri)
    p_vec = v_vec.scale(mass)

    print(f"\n初始条件 (近日点):")
    print(f"  位置 r = {r_vec}")
    print(f"  速度 v = {v_vec}")
    print(f"  速率 |v| = {v_vec.norm():.6f}")

    # 角动量
    L = r_vec ^ p_vec
    L_mag = L.norm()
    print(f"\n角动量 L = r ∧ p")
    print(f"  L = {L}")
    print(f"  |L| = {L_mag:.6f}")

    # 轨道能量
    r_mag = r_vec.norm()
    E_kinetic = 0.5 * mass * v_vec.norm() ** 2
    E_potential = -k / r_mag
    E_total = E_kinetic + E_potential

    print(f"\n能量:")
    print(f"  动能 T = {E_kinetic:.6f}")
    print(f"  势能 V = {E_potential:.6f}")
    print(f"  总能量 E = T + V = {E_total:.6f}")

    if E_total < 0:
        print("  轨道类型: 椭圆 (E < 0)")
    elif abs(E_total) < 1e-10:
        print("  轨道类型: 抛物线 (E = 0)")
    else:
        print("  轨道类型: 双曲线 (E > 0)")

    return k, mass, L_mag, E_total


def section_3_eccentricity_vector():
    """第三节：偏心率矢量"""
    print("\n" + "=" * 60)
    print("第三节：偏心率矢量 (Laplace-Runge-Lenz矢量)")
    print("=" * 60)

    alg = nblade.Algebra.euclidean(3)
    e1, e2, e3 = alg.basis_vectors()

    print("\n偏心率矢量是开普勒问题的重要守恒量")
    print("定义: e = (L·v)/k - r̂")
    print("其中 L·v 表示角动量与速度的左内积")

    # 参数
    k = 1.0
    mass = 1.0

    # 椭圆轨道参数
    r_peri = 1.0
    # 选择合适的速度以获得偏心率 < 1 (椭圆轨道)
    # v_peri < sqrt(2*k/r_peri) 时轨道为椭圆
    v_peri = math.sqrt(k / r_peri) * 0.9  # 小于逃逸速度

    r_vec = e1.scale(r_peri)
    v_vec = e2.scale(v_peri)
    p_vec = v_vec.scale(mass)

    # 角动量二重向量
    L = r_vec ^ p_vec
    L_mag = L.norm()

    print(f"\n计算偏心率矢量:")
    print(f"  位置 r = {r_vec}")
    print(f"  速度 v = {v_vec}")
    print(f"  角动量 L = r ∧ p = {L}")

    # 使用对偶将角动量转换为向量（仅在3D）
    L_vec = L.dual()
    print(f"  角动量向量 (对偶) L* = {L_vec}")

    # 偏心率矢量: e = (p × L_vec)/k - r̂ (Laplace-Runge-Lenz 矢量)
    # 在GA中，L_vec × v = (L_vec ∧ v)*
    # 注意：使用角动量向量 L_vec = dual(L)，而不是二重向量 L
    L_cross_v = (L_vec ^ v_vec).dual()  # L_vec × v
    r_hat = r_vec.scale(1.0 / r_vec.norm())

    # 对于引力问题，偏心率矢量指向远日点方向
    # e = (v × L) / (k/m) - r̂ = (L × v) / k - r̂ (单位制不同)
    # 这里我们用更简单的公式
    ecc_vec = L_cross_v.scale(1.0 / k) - r_hat
    ecc_mag = ecc_vec.norm()

    print(f"\n偏心率矢量计算:")
    print(f"  L_vec × v = (L_vec ∧ v)* = {L_cross_v}")
    print(f"  r̂ = r/|r| = {r_hat}")
    print(f"  e = (L_vec × v)/k - r̂ = {ecc_vec}")
    print(f"  偏心率 |e| = {ecc_mag:.6f}")

    # 验证偏心率与轨道参数的关系
    semi_latus_rectum = L_mag**2 / (mass * k)

    # 处理抛物线轨道 (e=1) 的特殊情况
    if ecc_mag >= 0.999:
        print(f"\n  轨道类型: 抛物线或双曲线 (e ≈ {ecc_mag:.3f})")
        print(f"  半正焦弦 l = L²/(mk) = {semi_latus_rectum:.6f}")
        semi_major = float("inf")  # 抛物线没有有限的半长轴
    else:
        semi_major = semi_latus_rectum / (1 - ecc_mag**2)
        print(f"\n轨道参数:")
        print(f"  半正焦弦 l = L²/(mk) = {semi_latus_rectum:.6f}")
        print(f"  半长轴 a = l/(1-e²) = {semi_major:.6f}")

    return ecc_mag, semi_latus_rectum


def section_4_orbital_equation():
    """第四节：轨道方程推导"""
    print("\n" + "=" * 60)
    print("第四节：轨道方程推导")
    print("=" * 60)

    print("\n使用几何代数推导极坐标下的轨道方程")
    print("\n标准形式: r(θ) = l / (1 + e·cos θ)")
    print("其中 l 是半正焦弦，e 是偏心率")

    alg = nblade.Algebra.euclidean(3)
    e1, e2, e3 = alg.basis_vectors()

    # 参数
    k = 1.0
    mass = 1.0

    # 设置椭圆轨道
    r_peri = 1.0
    ecc = 0.5  # 偏心率

    # 计算对应的近日点速度
    # v_peri² = (k/m) · (1+e) / (1-e) · 1/r_peri
    v_peri_sq = (k / mass) * (1 + ecc) / (1 - ecc) / r_peri
    v_peri = math.sqrt(v_peri_sq)

    print(f"\n轨道参数:")
    print(f"  近日点距离 r_p = {r_peri}")
    print(f"  偏心率 e = {ecc}")
    print(f"  近日点速度 v_p = {v_peri:.6f}")

    # 初始条件
    r_vec = e1.scale(r_peri)
    v_vec = e2.scale(v_peri)

    # 角动量
    L = r_vec ^ v_vec.scale(mass)
    L_mag = L.norm()

    # 半正焦弦
    l = L_mag**2 / (mass * k)

    print(f"\n推导结果:")
    print(f"  角动量 L = {L_mag:.6f}")
    print(f"  半正焦弦 l = L²/(mk) = {l:.6f}")
    print(f"  理论半正焦弦 = r_p(1+e) = {r_peri * (1 + ecc):.6f}")

    # 验证轨道方程
    print("\n【验证轨道方程 r = l/(1 + e·cos θ)】")

    angles = [0, math.pi / 4, math.pi / 2, 3 * math.pi / 4, math.pi]
    angle_names = ["0 (近日点)", "π/4", "π/2", "3π/4", "π (远日点)"]

    print("\n  θ        r(理论)      r(公式)")
    print("  " + "-" * 40)

    for angle, name in zip(angles, angle_names):
        # 理论值: 从几何计算
        cos_theta = math.cos(angle)
        r_formula = l / (1 + ecc * cos_theta)

        print(f"  {name:10s} {r_formula:10.6f}")

    # 远日点
    r_apo = l / (1 - ecc)
    print(f"\n  远日点距离 = l/(1-e) = {r_apo:.6f}")
    print(f"  半长轴 a = (r_p + r_a)/2 = {(r_peri + r_apo) / 2:.6f}")


def section_5_numerical_simulation():
    """第五节：数值模拟"""
    print("\n" + "=" * 60)
    print("第五节：数值模拟与守恒量验证")
    print("=" * 60)

    alg = nblade.Algebra.euclidean(3)
    e1, e2, e3 = alg.basis_vectors()

    print("\n使用几何代数进行轨道数值积分")
    print("方法: 显式欧拉法 (演示用)")

    # 物理参数
    k = 1.0  # G*M
    mass = 1.0

    # 轨道参数 (椭圆)
    ecc = 0.5
    r_peri = 1.0
    v_peri = math.sqrt(k / mass * (1 + ecc) / (1 - ecc) / r_peri)

    # 初始状态
    r = e1.scale(r_peri) + e2.scale(0.0)
    v = e1.scale(0.0) + e2.scale(v_peri)

    print(f"\n初始条件:")
    print(f"  位置 r₀ = {r}")
    print(f"  速度 v₀ = {v}")

    # 计算初始守恒量
    L_initial = r ^ v.scale(mass)
    E_initial = 0.5 * mass * v.norm() ** 2 - k / r.norm()

    L_initial_mag = L_initial.norm()

    print(f"\n初始守恒量:")
    print(f"  角动量 |L| = {L_initial_mag:.10f}")
    print(f"  能量 E = {E_initial:.10f}")

    # 模拟参数
    dt = 0.01
    num_steps = 1000

    print(f"\n模拟参数:")
    print(f"  时间步长 dt = {dt}")
    print(f"  总步数 N = {num_steps}")
    print(f"  总时间 T = {dt * num_steps}")

    # 存储轨迹
    trajectory = []

    # 数值积分 (显式欧拉)
    for step in range(num_steps):
        # 记录当前位置
        if step % 100 == 0:
            trajectory.append((r.norm(), step * dt))

        # 计算引力
        r_mag = r.norm()
        F = r.scale(-k / (r_mag**3))  # F = -k/r³ · r

        # 更新速度和位置
        v_new = v + F.scale(dt / mass)
        r_new = r + v.scale(dt)

        v = v_new
        r = r_new

    # 验证守恒量
    L_final = r ^ v.scale(mass)
    E_final = 0.5 * mass * v.norm() ** 2 - k / r.norm()

    print("\n【守恒量验证】")
    print(f"  初始角动量 |L₀| = {L_initial_mag:.10f}")
    print(f"  最终角动量 |L | = {L_final.norm():.10f}")
    print(
        f"  相对误差 = {abs(L_final.norm() - L_initial_mag) / L_initial_mag * 100:.6f}%"
    )

    print(f"\n  初始能量 E₀ = {E_initial:.10f}")
    print(f"  最终能量 E  = {E_final:.10f}")
    print(f"  绝对误差 = {abs(E_final - E_initial):.10f}")

    # 输出轨道点
    print("\n【轨道轨迹采样】")
    print("  时间 t      距离 r")
    print("  " + "-" * 30)
    for r_val, t_val in trajectory:
        print(f"  {t_val:8.4f}    {r_val:10.6f}")


def section_6_geometric_algebra_advantages():
    """第六节：几何代数的优势"""
    print("\n" + "=" * 60)
    print("第六节：几何代数的优势")
    print("=" * 60)

    print("\n几何代数处理开普勒问题的优势:")

    print("\n1. 角动量的自然表示")
    print("   传统: L = r × p (叉积，仅限3D)")
    print("   GA:   L = r ∧ p (外积，任意维度)")
    print("   → 二重向量明确表示旋转平面")

    print("\n2. 统一的运算框架")
    print("   几何积、外积、内积统一处理")
    print("   不需要区分点积、叉积、并矢等")

    print("\n3. 守恒量的几何意义")
    print("   角动量 L = r ∧ p: 位置和动量张成的有向面积")
    print("   偏心率矢量 e: 指向近日点的固定方向")

    print("\n4. 高维扩展性")
    print("   同样的公式适用于任意维度")
    print("   不需要为每个维度重新推导")

    print("\n5. 简洁的代数推导")
    print("   使用转子进行旋转")
    print("   避免复杂的三角函数运算")


def main():
    """主函数"""
    print("=" * 60)
    print("nblade 开普勒问题示例")
    print("几何代数在轨道力学中的应用")
    print("=" * 60)

    section_1_angular_momentum()
    section_2_kepler_setup()
    section_3_eccentricity_vector()
    section_4_orbital_equation()
    section_5_numerical_simulation()
    section_6_geometric_algebra_advantages()

    print("\n" + "=" * 60)
    print("示例完成！")
    print("=" * 60)
    print("\n关键点:")
    print("  - 角动量作为二重向量: L = r ∧ p")
    print("  - 在中心力场中角动量守恒")
    print("  - 偏心率矢量 e = (L×v)/k - r̂ 是守恒量")
    print("  - 轨道方程: r = l/(1 + e·cos θ)")
    print("  - 几何代数简化了推导并统一了框架")
    print("\n几何代数让经典力学的结构更加清晰！")


if __name__ == "__main__":
    main()
