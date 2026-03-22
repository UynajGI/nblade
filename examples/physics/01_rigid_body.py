"""
nblade 刚体物理示例
===================

演示几何代数在刚体动力学中的应用：
- 角动量作为二重向量
- 惯性张量
- 刚体旋转
- 欧拉方程的几何代数形式

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

    print("\n在传统向量分析中，角动量 L = r × p")
    print("在几何代数中，角动量自然表示为二重向量 L = r ∧ p")

    # 位置向量
    r = alg.vector([3.0, 4.0, 0.0])
    print(f"\n位置向量 r = {r}")

    # 动量向量
    p = alg.vector([2.0, -1.0, 5.0])
    print(f"动量向量 p = {p}")

    # 角动量作为二重向量
    L_bivector = r ^ p
    print(f"\n角动量二重向量 L = r∧p = {L_bivector}")

    # 与传统叉积的关系
    print("\n【与传统叉积的关系】")
    # 在 3D 中，L_dual = *(r∧p) 对应于传统的 L = r×p
    L_dual = L_bivector.dual()
    print(f"角动量对偶 (对应叉积) = {L_dual}")

    # 分解角动量的分量
    print("\n【角动量分量】")
    # L = L₁₂(e₁∧e₂) + L₁₃(e₁∧e₃) + L₂₃(e₂∧e₃)
    L_xy = (L_bivector | (e1 ^ e2)).scalar_part()
    L_xz = (L_bivector | (e1 ^ e3)).scalar_part()
    L_yz = (L_bivector | (e2 ^ e3)).scalar_part()

    print(f"  L_xy (xy 平面分量): {L_xy:.6f}")
    print(f"  L_xz (xz 平面分量): {L_xz:.6f}")
    print(f"  L_yz (yz 平面分量): {L_yz:.6f}")


def section_2_inertia_tensor():
    """第二节：惯性张量"""
    print("\n" + "=" * 60)
    print("第二节：惯性张量")
    print("=" * 60)

    alg = nblade.Algebra.euclidean(3)
    e1, e2, e3 = alg.basis_vectors()

    print("\n惯性张量将角速度映射到角动量")
    print("在几何代数中，我们可以使用双重向量运算")

    # 简单立方体的惯性张量 (对角)
    # 对于边长为 a 的立方体，关于中心的惯性矩
    # I_xx = I_yy = I_zz = (1/6) m a²
    mass = 1.0
    a = 1.0
    I_diag = mass * a * a / 6

    print(f"\n假设单位质量的立方体 (m=1, a=1)")
    print(f"惯性矩 I = {I_diag:.6f}")

    # 角速度向量
    omega = alg.vector([0.0, 0.0, 1.0])  # 绕 z 轴旋转
    print(f"\n角速度 ω = {omega}")

    # 角动量 L = I·ω (对于简单情况)
    L = omega.scale(I_diag)
    print(f"角动量 L = I·ω = {L}")

    # 一般惯性张量
    print("\n【一般惯性张量】")
    print("对于任意形状的刚体，惯性张量是 3×3 矩阵")
    # 使用度量张量函数演示
    vectors = [e1, e2, e3]
    g = nblade.metric_tensor(vectors)
    print("度量张量 g_ij:")
    for row in g:
        print(f"  {row}")


def section_3_rigid_body_rotation():
    """第三节：刚体旋转动力学"""
    print("\n" + "=" * 60)
    print("第三节：刚体旋转动力学")
    print("=" * 60)

    alg = nblade.Algebra.euclidean(3)
    e1, e2, e3 = alg.basis_vectors()

    print("\n刚体的旋转可以用转子描述")
    print("角速度对应于转子的时间导数")

    # 初始姿态
    initial_rotor = alg.rotor(e1 ^ e2, 0)
    print(f"\n初始姿态: R₀ = {initial_rotor}")

    # 角速度 (绕 z 轴)
    omega_magnitude = 1.0  # rad/s
    omega = alg.vector([0, 0, omega_magnitude])
    print(f"角速度 ω = {omega}")

    # 计算旋转轴
    # 在几何代数中，角速度平面是 ω 的对偶
    omega_plane = omega.dual()  # e1 ∧ e2
    print(f"角速度平面 = {omega_plane}")

    # 时间步进模拟
    print("\n【时间步进模拟】")
    dt = 0.1  # 时间步长
    num_steps = 10

    rotor = initial_rotor
    for step in range(num_steps):
        # 更新转子
        # 小角度近似: ΔR ≈ (1 + ω·dt/2)
        # 精确: R(t+dt) = exp(-ω·dt/2) · R(t)
        delta_angle = omega_magnitude * dt
        delta_rotor = alg.rotor(e1 ^ e2, delta_angle)
        rotor = delta_rotor * rotor

        # 计算旋转后的向量
        rotated = e1.rotate_by(rotor)
        print(f"t = {step * dt:.1f}s: e1 → {rotated}")


def section_4_euler_equations():
    """第四节：欧拉方程的几何代数形式"""
    print("\n" + "=" * 60)
    print("第四节：欧拉方程的几何代数形式")
    print("=" * 60)

    alg = nblade.Algebra.euclidean(3)
    e1, e2, e3 = alg.basis_vectors()

    print("\n欧拉方程描述无外力矩刚体的旋转:")
    print("  dL/dt + ω × L = 0")
    print("\n在几何代数中，可以写成:")
    print("  dL/dt + ω ∧ L = 0 (使用外积)")

    # 惯性主轴
    I1, I2, I3 = 1.0, 2.0, 3.0  # 不同的主惯性矩

    print(f"\n假设不对称刚体:")
    print(f"  主惯性矩 I₁={I1}, I₂={I2}, I₃={I3}")

    # 角速度分量
    omega1, omega2, omega3 = 1.0, 0.5, 0.2
    omega = e1.scale(omega1) + e2.scale(omega2) + e3.scale(omega3)
    print(f"\n角速度 ω = {omega1}e₁ + {omega2}e₂ + {omega3}e₃")

    # 角动量 L = I₁ω₁e₁ + I₂ω₂e₂ + I₃ω₃e₃
    L = e1.scale(I1 * omega1) + e2.scale(I2 * omega2) + e3.scale(I3 * omega3)
    print(f"角动量 L = {L}")

    # ω × L 的几何代数表示
    omega_cross_L = (omega ^ L).dual()
    print(f"\nω × L = {omega_cross_L}")

    print("\n注意: 这表示角动量变化率，对于自由旋转的刚体")
    print("角动量矢量在体坐标系中会随时间变化")


def section_5_energy():
    """第五节：转动动能"""
    print("\n" + "=" * 60)
    print("第五节：转动动能")
    print("=" * 60)

    alg = nblade.Algebra.euclidean(3)
    e1, e2, e3 = alg.basis_vectors()

    print("\n转动动能 T = (1/2) ω · L = (1/2) Σ Iᵢωᵢ²")

    # 主惯性矩
    I1, I2, I3 = 1.0, 2.0, 3.0

    # 角速度
    omega1, omega2, omega3 = 1.0, 0.5, 0.2
    omega = e1.scale(omega1) + e2.scale(omega2) + e3.scale(omega3)

    # 角动量
    L = e1.scale(I1 * omega1) + e2.scale(I2 * omega2) + e3.scale(I3 * omega3)

    # 动能: T = (1/2) ω · L
    T = 0.5 * (omega | L).scalar_part()

    print(f"角速度 ω = {omega}")
    print(f"角动量 L = {L}")
    print(f"转动动能 T = {T:.6f}")

    # 验证: T = (1/2) Σ Iᵢωᵢ²
    T_verify = 0.5 * (I1 * omega1**2 + I2 * omega2**2 + I3 * omega3**2)
    print(f"验证: (1/2)ΣIᵢωᵢ² = {T_verify:.6f}")


def section_6_simulation():
    """第六节：简单刚体模拟"""
    print("\n" + "=" * 60)
    print("第六节：简单刚体模拟")
    print("=" * 60)

    alg = nblade.Algebra.euclidean(3)
    e1, e2, e3 = alg.basis_vectors()

    print("\n模拟一个立方体的自由旋转")

    # 参数
    mass = 1.0
    side = 1.0
    I = mass * side**2 / 6  # 惯性矩

    # 初始角速度
    omega = e1.scale(1.0) + e2.scale(0.5)
    print(f"初始角速度 ω = {omega}")
    print(f"惯性矩 I = {I:.6f}")

    # 初始角动量
    L = omega.scale(I)
    print(f"初始角动量 L = {L}")

    # 初始姿态 (转子)
    rotor = alg.one()

    # 模拟参数
    dt = 0.05
    num_steps = 20

    print(f"\n模拟 {num_steps} 步 (dt = {dt}s):")
    print("-" * 40)

    # 固定角动量，更新角速度和姿态
    for step in range(num_steps):
        # 当前角速度 (L/I)
        omega_current = L.scale(1.0 / I)

        # 计算旋转平面和角度
        omega_mag = omega_current.norm()
        if omega_mag > 1e-10:
            # 旋转轴方向
            axis = omega_current.scale(1.0 / omega_mag)

            # 小角度旋转
            delta_angle = omega_mag * dt

            # 构造增量转子
            # 使用角速度平面的对偶
            omega_plane = omega_current.dual()

            try:
                delta_rotor = alg.rotor(omega_plane, delta_angle)
                rotor = delta_rotor * rotor
            except:
                # 如果 rotor 创建失败，跳过
                pass

        # 打印状态
        if step % 5 == 0:
            rotated_e1 = e1.rotate_by(rotor)
            print(f"t = {step * dt:.2f}s: e1 体坐标 = {rotated_e1}")

    print("\n模拟完成！")
    print("注意: 角动量守恒，但体坐标系中的角速度方向可能变化")


def main():
    """主函数"""
    print("=" * 60)
    print("nblade 刚体物理示例")
    print("几何代数在经典力学中的应用")
    print("=" * 60)

    section_1_angular_momentum()
    section_2_inertia_tensor()
    section_3_rigid_body_rotation()
    section_4_euler_equations()
    section_5_energy()
    section_6_simulation()

    print("\n" + "=" * 60)
    print("示例完成！")
    print("=" * 60)
    print("\n关键点:")
    print("  - 角动量作为二重向量: L = r ∧ p")
    print("  - 角速度与旋转平面的对偶关系")
    print("  - 转子描述刚体姿态")
    print("  - 几何代数简化了旋转动力学的表达")


if __name__ == "__main__":
    main()
