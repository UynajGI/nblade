"""
nblade 电磁学示例
================

演示几何代数在电磁学中的应用：
- 电磁场作为二重向量 F = E + iB
- 场不变量及其物理意义
- 洛伦兹力与带电粒子运动
- 麦克斯韦方程组的几何代数形式
- 平面波解

在时空代数(STA)中，电磁场优雅地表示为一个二重向量，
麦克斯韦方程组简化为单一方程：∇F = J

本示例使用3D欧几里得代数展示基本概念，
并在适当位置讨论时空代数 formulation。

作者: nblade 团队
"""

import nblade
import math


def section_1_electromagnetic_field():
    """第一节：电磁场作为二重向量"""
    print("\n" + "=" * 60)
    print("第一节：电磁场作为二重向量")
    print("=" * 60)

    alg = nblade.Algebra.euclidean(3)
    e1, e2, e3 = alg.basis_vectors()
    I = alg.config.volume_element()  # 伪标量

    print("\n在传统向量分析中，电磁场用两个分离的向量表示：")
    print("  电场 E (极向量)")
    print("  磁场 B (轴向量/伪向量)")
    print("\n在几何代数中，电磁场统一为一个二重向量：")
    print("  F = E + IB")
    print("其中 I 是伪标量，将磁场转换为二重向量")

    # 定义电场分量
    E_x, E_y, E_z = 1.0, 2.0, 3.0
    E = e1.scale(E_x) + e2.scale(E_y) + e3.scale(E_z)

    print(f"\n【电场】")
    print(f"  E = {E_x}e₁ + {E_y}e₂ + {E_z}e₃")
    print(f"  E = {E}")

    # 定义磁场分量
    B_x, B_y, B_z = 0.5, 1.0, 1.5
    B = e1.scale(B_x) + e2.scale(B_y) + e3.scale(B_z)

    print(f"\n【磁场】")
    print(f"  B = {B_x}e₁ + {B_y}e₂ + {B_z}e₃")
    print(f"  B = {B}")

    # 伪标量 I = e₁ ∧ e₂ ∧ e₃
    print(f"\n【伪标量】")
    print(f"  I = e₁ ∧ e₂ ∧ e₃ (体积元素)")
    print(f"  I² = -1 (在欧几里得空间)")

    # 构造电磁场二重向量
    # 在3D中，IB 将向量 B 转换为二重向量
    IB = I * B  # 伪标量与向量的几何积

    print(f"\n【磁场转换为二重向量】")
    print(f"  IB = I · B = {IB}")

    # 完整的电磁场二重向量 F = E + IB
    # 注意：这里 E 是向量，需要先转换为二重向量
    # 在3D欧几里得代数中，电场作为向量，磁场通过伪标量转换
    print(f"\n【电磁场二重向量】")
    print(f"  F = E + IB (概念表示)")
    print(f"  其中 E 作为向量部分，IB 作为二重向量部分")

    print("\n【物理意义】")
    print("  F 统一了电场和磁场")
    print("  不同参考系下的观察者看到不同的 E 和 B 分量")
    print("  但 F 本身是一个协变的几何对象")


def section_2_field_invariants():
    """第二节：场不变量"""
    print("\n" + "=" * 60)
    print("第二节：场不变量")
    print("=" * 60)

    alg = nblade.Algebra.euclidean(3)
    e1, e2, e3 = alg.basis_vectors()
    I = alg.config.volume_element()

    print("\n电磁场有两个洛伦兹不变量：")
    print("  1. F · F̃ = E² - B² (标量)")
    print("  2. F ∧ F̃ = E · B (伪标量)")
    print("其中 F̃ 是 F 的对偶")

    # 示例1：纯电场
    print("\n【示例1：纯电场】")
    E = e1.scale(2.0)
    B = e1.scale(0.0) + e2.scale(0.0) + e3.scale(0.0)

    E_sq = (E | E).scalar_part()
    B_sq = (B | B).scalar_part()

    print(f"  E = 2e₁")
    print(f"  B = 0")
    print(f"  E² = {E_sq:.6f}")
    print(f"  B² = {B_sq:.6f}")
    print(f"  E² - B² = {E_sq - B_sq:.6f} (不变量)")
    print(f"  E·B = 0 (不变量)")

    # 示例2：纯磁场
    print("\n【示例2：纯磁场】")
    E = e1.scale(0.0)
    B = e2.scale(3.0)

    E_sq = (E | E).scalar_part()
    B_sq = (B | B).scalar_part()

    print(f"  E = 0")
    print(f"  B = 3e₂")
    print(f"  E² = {E_sq:.6f}")
    print(f"  B² = {B_sq:.6f}")
    print(f"  E² - B² = {E_sq - B_sq:.6f} (不变量)")

    # 示例3：正交电场和磁场
    print("\n【示例3：正交电场和磁场】")
    E = e1.scale(4.0)
    B = e2.scale(3.0)

    E_sq = (E | E).scalar_part()
    B_sq = (B | B).scalar_part()
    E_dot_B = (E | B).scalar_part()

    print(f"  E = 4e₁")
    print(f"  B = 3e₂ (与 E 正交)")
    print(f"  E² - B² = {E_sq - B_sq:.6f}")
    print(f"  E·B = {E_dot_B:.6f}")

    # 示例4：平行电场和磁场
    print("\n【示例4：平行电场和磁场】")
    E = e1.scale(4.0)
    B = e1.scale(3.0)

    E_sq = (E | E).scalar_part()
    B_sq = (B | B).scalar_part()
    E_dot_B = (E | B).scalar_part()

    print(f"  E = 4e₁")
    print(f"  B = 3e₁ (与 E 平行)")
    print(f"  E² - B² = {E_sq - B_sq:.6f}")
    print(f"  E·B = {E_dot_B:.6f}")

    print("\n【物理意义】")
    print("  E² - B² = 0 → 电磁波或零场")
    print("  E·B = 0 → 正交场，可找到某参考系使 E 或 B 为零")
    print("  这些不变量决定了场的类型和能量密度")


def section_3_lorentz_force():
    """第三节：洛伦兹力"""
    print("\n" + "=" * 60)
    print("第三节：洛伦兹力")
    print("=" * 60)

    alg = nblade.Algebra.euclidean(3)
    e1, e2, e3 = alg.basis_vectors()
    I = alg.config.volume_element()

    print("\n洛伦兹力描述带电粒子在电磁场中的运动：")
    print("  F = q(E + v × B)")
    print("\n在几何代数中，叉积用对偶和外积表示：")
    print("  v × B = (v ∧ IB)⁻")
    print("  或等价地：v × B = -I(v ∧ B)")

    # 粒子参数
    q = 1.0  # 电荷
    m = 1.0  # 质量

    # 速度
    v = e1.scale(1.0) + e2.scale(0.5)
    print(f"\n【粒子状态】")
    print(f"  电荷 q = {q}")
    print(f"  质量 m = {m}")
    print(f"  速度 v = {v}")

    # 电场
    E = e3.scale(2.0)
    print(f"  电场 E = {E}")

    # 磁场
    B = e3.scale(1.0)
    print(f"  磁场 B = {B}")

    # 电场力
    F_electric = E.scale(q)
    print(f"\n【电场力】")
    print(f"  F_E = qE = {F_electric}")

    # 磁场力 (v × B)
    # 在GA中：v × B = dual(v ∧ B) 在3D中
    v_wedge_B = v ^ B
    v_cross_B = v_wedge_B.dual()

    print(f"\n【磁场力计算】")
    print(f"  v ∧ B = {v_wedge_B}")
    print(f"  v × B = (v ∧ B)* = {v_cross_B}")

    F_magnetic = v_cross_B.scale(q)
    print(f"  F_B = q(v × B) = {F_magnetic}")

    # 总洛伦兹力
    F_total = F_electric + F_magnetic
    print(f"\n【总洛伦兹力】")
    print(f"  F = q(E + v × B) = {F_total}")

    # 加速度
    a = F_total.scale(1.0 / m)
    print(f"  加速度 a = F/m = {a}")

    print("\n【几何代数的优势】")
    print("  - 叉积自然表示为对偶的外积")
    print("  - 统一处理电场和磁场力")
    print("  - 公式在任意维度都成立（除了对偶操作）")


def section_4_maxwell_equations():
    """第四节：麦克斯韦方程组的几何代数形式"""
    print("\n" + "=" * 60)
    print("第四节：麦克斯韦方程组的几何代数形式")
    print("=" * 60)

    print("\n传统麦克斯韦方程组（4个方程）：")
    print("  ∇·E = ρ/ε₀      (高斯定律)")
    print("  ∇·B = 0         (无磁单极)")
    print("  ∇×E = -∂B/∂t    (法拉第定律)")
    print("  ∇×B = μ₀J + μ₀ε₀∂E/∂t  (安培-麦克斯韦定律)")

    print("\n在时空代数(STA)中，合并为一个优雅的方程：")
    print("  ∇F = J")
    print("\n其中：")
    print("  ∇ = γ^μ ∂_μ (时空导数)")
    print("  F = E + IB (电磁场二重向量)")
    print("  J = 时空电流矢量")

    print("\n【方程结构】")
    print("  ∇F = ∇·F + ∇∧F")
    print("  = J (矢量部分) + 0 (三重向量部分)")

    print("\n【分解】")
    print("  矢量部分 (∇·F)：")
    print("    → ∇·E = ρ (高斯定律)")
    print("    → ∇×B - ∂E/∂t = J (安培定律)")
    print("\n  三重向量部分 (∇∧F = 0)：")
    print("    → ∇·B = 0 (无磁单极)")
    print("    → ∇×E + ∂B/∂t = 0 (法拉第定律)")

    print("\n【几何意义】")
    print("  单一方程揭示了电磁学的几何结构：")
    print("  - F 是二重向量（有向面积元）")
    print("  - ∇F = J 表示场与源的关系")
    print("  - 洛伦兹协变性自动满足")

    print("\n【规范不变性】")
    print("  势形式：F = ∇A")
    print("  其中 A 是电磁势（矢量+标量部分）")
    print("  规范变换：A → A + ∇χ 不改变 F")

    # 简单示例：验证连续性方程
    print("\n【连续性方程】")
    print("  从 ∇F = J 和 ∇∧F = 0，可以得到：")
    print("  ∇·J = 0 (电荷守恒)")

    alg = nblade.Algebra.euclidean(3)
    e1, e2, e3 = alg.basis_vectors()

    # 示例电流分布
    J = e1.scale(1.0) + e2.scale(2.0) + e3.scale(3.0)
    print(f"\n  示例电流密度 J = {J}")
    print(f"  ∇·J = 0 表示稳态电流（无源无汇）")


def section_5_plane_wave():
    """第五节：平面波解"""
    print("\n" + "=" * 60)
    print("第五节：平面波解")
    print("=" * 60)

    alg = nblade.Algebra.euclidean(3)
    e1, e2, e3 = alg.basis_vectors()
    I = alg.config.volume_element()

    print("\n电磁波是麦克斯韦方程组的重要解")
    print("在真空中，平面波满足：")
    print("  E ⊥ B ⊥ k (传播方向)")
    print("  |E| = |B| (在自然单位制)")

    # 波参数
    wavelength = 1.0
    k = 2 * math.pi / wavelength  # 波数
    omega = k  # 角频率（光速 c=1）

    print(f"\n【波参数】")
    print(f"  波长 λ = {wavelength}")
    print(f"  波数 k = 2π/λ = {k:.6f}")
    print(f"  角频率 ω = {omega:.6f} (c=1)")

    # 沿 e3 方向传播的波
    # 电场沿 e1，磁场沿 e2
    print(f"\n【场配置】")
    print(f"  传播方向: e₃ (z轴)")
    print(f"  电场方向: e₁ (x轴)")
    print(f"  磁场方向: e₂ (y轴)")

    # 计算不同时刻的场值
    times = [0.0, 0.125, 0.25, 0.375, 0.5]

    print(f"\n【场随时间演化】")
    print(f"  时间 t      E 幅值      B 幅值      相位")
    print("  " + "-" * 50)

    for t in times:
        phase = omega * t
        E_amp = math.cos(phase)
        B_amp = math.cos(phase)

        print(f"  {t:.3f}       {E_amp:+.6f}    {B_amp:+.6f}    {phase / math.pi:.3f}π")

    # 验证 E ⊥ B
    print(f"\n【正交性验证】")
    E_wave = e1  # 电场方向
    B_wave = e2  # 磁场方向

    E_dot_B = (E_wave | B_wave).scalar_part()
    print(f"  E·B = {E_dot_B:.10f} (应为零)")

    # 验证 E ⊥ k 和 B ⊥ k
    k_dir = e3  # 传播方向

    E_dot_k = (E_wave | k_dir).scalar_part()
    B_dot_k = (B_wave | k_dir).scalar_part()

    print(f"  E·k = {E_dot_k:.10f} (应为零)")
    print(f"  B·k = {B_dot_k:.10f} (应为零)")

    print("\n【电磁场的二重向量表示】")
    print("  F(t,z) = E₀ cos(kz - ωt) e₁ + I B₀ cos(kz - ωt) e₂")
    print("  F = E + IB，统一表示电磁波")

    print("\n【能量密度】")
    E_sq = (E_wave | E_wave).scalar_part()
    B_sq = (B_wave | B_wave).scalar_part()
    u = 0.5 * (E_sq + B_sq)
    print(f"  能量密度 u = (E² + B²)/2 = {u:.6f}")
    print(f"  在真空中 |E| = |B|，所以 u = E² = B²")

    print("\n【能流密度 (坡印廷矢量)】")
    print("  S = E × B = (E ∧ B)*")
    E_cross_B = (E_wave ^ B_wave).dual()
    print(f"  S = E × B = {E_cross_B}")
    print(f"  能流方向与传播方向一致")


def section_6_ga_advantages():
    """第六节：几何代数在电磁学中的优势"""
    print("\n" + "=" * 60)
    print("第六节：几何代数在电磁学中的优势")
    print("=" * 60)

    print("\n几何代数为电磁学提供了统一且优雅的框架：")

    print("\n【1. 场的统一表示】")
    print("  传统: E 和 B 作为分离的向量")
    print("  GA:   F = E + IB 作为单一的二重向量")
    print("  → 场变换更简洁，协变性更明显")

    print("\n【2. 麦克斯韦方程组的简化】")
    print("  传统: 4个微分方程")
    print("  GA:   ∇F = J (单一方程)")
    print("  → 结构更清晰，物理意义更明确")

    print("\n【3. 洛伦兹力的自然表示】")
    print("  传统: F = q(E + v × B)")
    print("  GA:   dp/dτ = qF·v (时空形式)")
    print("  → 相对论协变性自动满足")

    print("\n【4. 势的几何意义】")
    print("  传统: A = (φ, A) 作为标量+矢量势")
    print("  GA:   A 作为时空矢量")
    print("  → F = ∇∧A 明确表示场是势的导数")

    print("\n【5. 规范不变性】")
    print("  传统: 需要分离讨论电势和磁势")
    print("  GA:   A → A + ∇χ 自然地统一")
    print("  → 规范变换是简单的几何操作")

    print("\n【6. 电磁波与场不变量】")
    print("  平面波 F² = 0 (零场)")
    print("  驻波 F² > 0 或 F² < 0 (非零不变量)")
    print("  → 场分类变得简单")

    print("\n【7. 相对论变换】")
    print("  传统: 需要复杂的变换矩阵")
    print("  GA:   F' = RFR⁻ (转子变换)")
    print("  → 洛伦兹变换统一为转子运算")

    print("\n【8. 与其他物理理论的统一】")
    print("  GA 同样适用于：")
    print("  - 狭义相对论")
    print("  - 量子力学 (狄拉克方程)")
    print("  - 引力理论")
    print("  → 统一的数学语言")

    # 演示：不同场类型的 F²
    print("\n【场类型示例】")
    alg = nblade.Algebra.euclidean(3)
    e1, e2, e3 = alg.basis_vectors()

    print("\n  场类型         E² vs B²      E·B      特点")
    print("  " + "-" * 55)
    print("  纯电场         E² > B²       0        静电场")
    print("  纯磁场         E² < B²       0        静磁场")
    print("  电磁波         E² = B²       0        辐射场")
    print("  斜交场         任意          0        可找到零场参考系")
    print("  复杂场         任意          ≠0       无法简化为零场")


def main():
    """主函数"""
    print("=" * 60)
    print("nblade 电磁学示例")
    print("几何代数在电磁理论中的应用")
    print("=" * 60)

    section_1_electromagnetic_field()
    section_2_field_invariants()
    section_3_lorentz_force()
    section_4_maxwell_equations()
    section_5_plane_wave()
    section_6_ga_advantages()

    print("\n" + "=" * 60)
    print("示例完成！")
    print("=" * 60)
    print("\n关键点:")
    print("  - 电磁场统一为二重向量 F = E + IB")
    print("  - 场不变量 E² - B² 和 E·B 具有重要物理意义")
    print("  - 洛伦兹力自然表示为几何代数运算")
    print("  - 麦克斯韦方程组简化为 ∇F = J")
    print("  - 平面波是电磁场的基本解")
    print("  - 几何代数提供了电磁学的统一框架")
    print("\n几何代数让电磁学的结构更加清晰优雅！")


if __name__ == "__main__":
    main()
