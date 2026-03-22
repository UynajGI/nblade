"""
nblade 时空代数示例
====================

演示几何代数在狭义相对论中的应用：
- 时空代数 G(1,3,0) 的创建与基本运算
- 四向量：位置、速度、动量
- 洛伦兹变换与转子（推进 + 旋转）
- 相对论速度叠加
- 时空区间与光锥
- 与狭义相对论的关系

作者: nblade 团队
"""

import nblade
import math


def section_1_spacetime_algebra():
    """第一节：创建时空代数 - G(1,3,0)"""
    print("\n" + "=" * 60)
    print("第一节：创建时空代数 G(1,3,0)")
    print("=" * 60)

    # 创建时空代数
    # G(1,3,0) 表示：1 个正平方基向量 (时间)，3 个负平方基向量 (空间)
    # 度量签名: (+, -, -, -)
    alg = nblade.Algebra.spacetime(4)

    print("\n【时空代数基础】")
    print(f"代数: {alg}")
    print(f"维度: {alg.dimension}")
    print(f"签名 (p, q, r): {alg.signature}")
    print("  p=1: 一个正平方基向量 (时间 e₀)")
    print("  q=3: 三个负平方基向量 (空间 e₁, e₂, e₃)")
    print("  r=0: 没有零平方基向量")

    # 获取基向量
    e0, e1, e2, e3 = alg.basis_vectors()

    print("\n【基向量及其平方】")
    # 时间基向量 e₀² = +1
    e0_sq = (e0 | e0).scalar_part()
    print(f"e₀² = {e0_sq}  (时间方向，正平方)")

    # 空间基向量 e₁² = e₂² = e₃² = -1
    e1_sq = (e1 | e1).scalar_part()
    e2_sq = (e2 | e2).scalar_part()
    e3_sq = (e3 | e3).scalar_part()
    print(f"e₁² = {e1_sq}  (空间方向，负平方)")
    print(f"e₂² = {e2_sq}  (空间方向，负平方)")
    print(f"e₃² = {e3_sq}  (空间方向，负平方)")

    # 验证度量张量
    print("\n【度量张量验证】")
    print("时空度规 η = diag(+1, -1, -1, -1)")
    print("任意向量 v = v⁰e₀ + v¹e₁ + v²e₂ + v³e₃")
    print("v² = (v⁰)² - (v¹)² - (v²)² - (v³)²")

    # 伪标量（体积元）
    I = alg.config.volume_element()
    print(f"\n【伪标量（体积元）】")
    print(f"I = e₀ ∧ e₁ ∧ e₂ ∧ e₃ = {I}")
    I_sq = (I | I).scalar_part()
    print(f"I² = {I_sq}")


def section_2_four_vectors():
    """第二节：四向量 - 位置、速度、动量"""
    print("\n" + "=" * 60)
    print("第二节：四向量")
    print("=" * 60)

    alg = nblade.Algebra.spacetime(4)
    e0, e1, e2, e3 = alg.basis_vectors()

    # 光速 (使用自然单位制 c = 1)
    c = 1.0

    print("\n【四位置向量】")
    print("在时空代数中，四位置 X = ct e₀ + x e₁ + y e₂ + z e₃")

    # 创建四位置 (ct, x, y, z) = (1.0, 2.0, 3.0, 0.0)
    ct = 1.0
    x, y, z = 2.0, 3.0, 0.0
    X = e0.scale(ct * c) + e1.scale(x) + e2.scale(y) + e3.scale(z)
    print(f"四位置 X = {ct}e₀ + {x}e₁ + {y}e₂ + {z}e₃")
    print(f"         = ct e₀ + x e₁ + y e₂ + z e₃")
    print(f"         = {X}")

    # 计算时空区间
    spacetime_interval = (X | X).scalar_part()
    print(f"\n时空区间 X² = (ct)² - x² - y² - z²")
    print(f"          = {ct * c}² - {x}² - {y}² - {z}²")
    print(f"          = {ct * c * ct * c - x * x - y * y - z * z}")
    print(f"X² = {spacetime_interval}")

    if spacetime_interval > 0:
        print("  → 类时间隔 (timelike)")
    elif spacetime_interval < 0:
        print("  → 类空间隔 (spacelike)")
    else:
        print("  → 类光间隔 (lightlike/null)")

    print("\n【四速度向量】")
    print("四速度 U = γ(c, vₓ, vᵧ, vᵤ)，其中 γ = 1/√(1-v²/c²)")

    # 三维速度 (以光速为单位)
    v_x, v_y, v_z = 0.6, 0.0, 0.0  # 沿 x 方向 0.6c
    v_sq = v_x**2 + v_y**2 + v_z**2
    gamma = 1.0 / math.sqrt(1 - v_sq)  # 洛伦兹因子

    print(f"三维速度: v = ({v_x}c, {v_y}c, {v_z}c)")
    print(f"|v| = {math.sqrt(v_sq):.4f}c")
    print(f"洛伦兹因子 γ = {gamma:.6f}")

    # 四速度
    U = (
        e0.scale(gamma * c)
        + e1.scale(gamma * v_x * c)
        + e2.scale(gamma * v_y * c)
        + e3.scale(gamma * v_z * c)
    )
    print(f"\n四速度 U = γ(c, vₓ, vᵧ, vᵤ)")
    print(f"         = {gamma:.4f}({c}, {v_x}c, {v_y}c, {v_z}c)")
    print(f"         = {U}")

    # 验证四速度归一化: U² = c²
    U_sq = (U | U).scalar_part()
    print(f"\n验证 U² = {U_sq:.6f} (应等于 c² = {c * c})")

    print("\n【四动量向量】")
    print("四动量 P = m₀U = E/c e₀ + pₓ e₁ + pᵧ e₂ + pᵤ e₃")

    # 静止质量 (自然单位)
    m0 = 1.0

    # 四动量
    E = gamma * m0 * c**2  # 总能量
    p_x = gamma * m0 * v_x * c  # 动量分量
    p_y = gamma * m0 * v_y * c
    p_z = gamma * m0 * v_z * c

    P = e0.scale(E / c) + e1.scale(p_x) + e2.scale(p_y) + e3.scale(p_z)

    print(f"静止质量 m₀ = {m0}")
    print(f"总能量 E = γm₀c² = {E:.6f}")
    print(f"动量 p = ({p_x:.6f}, {p_y:.6f}, {p_z:.6f})")
    print(f"\n四动量 P = E/c e₀ + pₓ e₁ + pᵧ e₂ + pᵤ e₃")
    print(f"         = {P}")

    # 验证质能关系: P² = m₀²c²
    P_sq = (P | P).scalar_part()
    expected = m0**2 * c**2
    print(f"\n验证 P² = {P_sq:.6f} (应等于 m₀²c² = {expected})")
    print(f"质能关系: E² - |p|²c² = m₀²c⁴")


def section_3_lorentz_transformations():
    """第三节：洛伦兹变换（推进和旋转）"""
    print("\n" + "=" * 60)
    print("第三节：洛伦兹变换")
    print("=" * 60)

    alg = nblade.Algebra.spacetime(4)
    e0, e1, e2, e3 = alg.basis_vectors()

    c = 1.0

    print("\n【洛伦兹推进（Boost）】")
    print("推进是沿某方向的洛伦兹变换")
    print("使用转子实现: X' = L X L~")

    # 推进速度 (沿 x 方向)
    beta = 0.5  # v/c
    v = beta * c

    print(f"\n推进速度: v = {beta}c (沿 x 方向)")

    # 洛伦兹因子
    gamma = 1.0 / math.sqrt(1 - beta**2)
    print(f"洛伦兹因子 γ = {gamma:.6f}")

    # 快度 (rapidity)
    # φ = artanh(β)
    rapidity = 0.5 * math.log((1 + beta) / (1 - beta))
    print(f"快度 φ = artanh(β) = {rapidity:.6f} rad")

    # 创建推进转子
    # 在时空代数中，推进转子为 L = exp(-φ e₀e₁/2)
    # 这对应于在 e₀e₁ 平面内的 "旋转"
    boost_plane = e0 ^ e1
    print(f"\n推进平面: {boost_plane}")

    # 使用转子构造推进
    # 注意: 推进角度使用快度
    try:
        boost_rotor = alg.rotor(boost_plane, rapidity)
        print(f"推进转子 L = exp(-φ e₀e₁/2)")
        print(f"           = {boost_rotor}")
    except Exception as ex:
        print(f"注意: 转子创建可能受限于库实现: {ex}")
        print("使用代数方法计算...")

    print("\n【手动计算推进变换】")
    # 洛伦兹变换矩阵形式
    print("洛伦兹变换沿 x 方向:")
    print("  ct' = γ(ct - βx)")
    print("  x'  = γ(x - βct)")
    print("  y'  = y")
    print("  z'  = z")

    # 示例：变换一个四位置
    ct_orig = 2.0
    x_orig = 1.0
    y_orig = 0.0
    z_orig = 0.0

    X_orig = e0.scale(ct_orig) + e1.scale(x_orig) + e2.scale(y_orig) + e3.scale(z_orig)

    print(f"\n原始四位置: X = ({ct_orig}, {x_orig}, {y_orig}, {z_orig})")

    # 应用洛伦兹变换
    ct_new = gamma * (ct_orig - beta * x_orig)
    x_new = gamma * (x_orig - beta * ct_orig)
    y_new = y_orig
    z_new = z_orig

    X_new = e0.scale(ct_new) + e1.scale(x_new) + e2.scale(y_new) + e3.scale(z_new)

    print(f"变换后: X' = ({ct_new:.6f}, {x_new:.6f}, {y_new:.6f}, {z_new:.6f})")

    # 验证时空区间不变
    interval_orig = (X_orig | X_orig).scalar_part()
    interval_new = (X_new | X_new).scalar_part()
    print(f"\n验证时空区间不变性:")
    print(f"  原始 X² = {interval_orig:.6f}")
    print(f"  变换后 X'² = {interval_new:.6f}")
    print(f"  差值 = {abs(interval_orig - interval_new):.2e}")

    print("\n【旋转与推进的结合】")
    print("完整的洛伦兹变换包含推进和旋转")
    print("在时空代数中，可以统一表示为转子变换")

    # 先推进后旋转
    print("\n变换顺序: X' = R L X L~ R~")
    print("其中 L 是推进转子，R 是旋转转子")

    # 创建 90 度旋转转子 (绕 z 轴)
    rotation_plane = e1 ^ e2
    rotation_angle = math.pi / 2
    rotation_rotor = alg.rotor(rotation_plane, rotation_angle)

    print(f"\n旋转转子 (绕 z 轴 90°): R = {rotation_rotor}")


def section_4_velocity_addition():
    """第四节：相对论速度叠加"""
    print("\n" + "=" * 60)
    print("第四节：相对论速度叠加")
    print("=" * 60)

    alg = nblade.Algebra.spacetime(4)
    e0, e1, e2, e3 = alg.basis_vectors()

    c = 1.0

    print("\n【相对论速度叠加公式】")
    print("在狭义相对论中，速度叠加遵循:")
    print("  u' = (u + v) / (1 + uv/c²)")
    print("而非伽利略的简单相加: u' = u + v")

    # 第一个参考系相对于静止系的速度
    beta1 = 0.6  # v1/c
    gamma1 = 1.0 / math.sqrt(1 - beta1**2)

    # 物体相对于第一个参考系的速度
    beta2 = 0.5  # v2/c
    gamma2 = 1.0 / math.sqrt(1 - beta2**2)

    print(f"\n参考系 S' 相对 S 以 v₁ = {beta1}c 运动 (沿 x 方向)")
    print(f"物体相对 S' 以 v₂ = {beta2}c 运动 (沿 x 方向)")

    # 相对论速度叠加
    # u' = (u + v) / (1 + uv/c²)
    beta_combined = (beta1 + beta2) / (1 + beta1 * beta2)
    print(f"\n相对论叠加:")
    print(f"  β₁ = {beta1}, β₂ = {beta2}")
    print(f"  β = (β₁ + β₂) / (1 + β₁β₂)")
    print(f"    = ({beta1} + {beta2}) / (1 + {beta1} × {beta2})")
    print(f"    = {beta_combined:.6f}")

    # 验证不超过光速
    print(f"\n验证: {beta_combined:.6f} < 1 (光速上限)")

    # 伽利略叠加 (错误)
    galilean = beta1 + beta2
    print(f"\n伽利略叠加 (错误): {galilean}")
    print(f"这会超过光速！")

    print("\n【使用几何代数验证】")
    # 计算合成洛伦兹因子
    gamma_combined = 1.0 / math.sqrt(1 - beta_combined**2)
    print(f"合成洛伦兹因子 γ = {gamma_combined:.6f}")

    # 另一种计算方法: γ = γ₁γ₂(1 + β₁β₂)
    gamma_product = gamma1 * gamma2 * (1 + beta1 * beta2)
    print(f"验证: γ₁γ₂(1 + β₁β₂) = {gamma_product:.6f}")

    # 快度叠加
    # 相对论中快度可以简单相加
    rapidity1 = 0.5 * math.log((1 + beta1) / (1 - beta1))
    rapidity2 = 0.5 * math.log((1 + beta2) / (1 - beta2))
    rapidity_sum = rapidity1 + rapidity2

    print(f"\n【快度叠加】")
    print(f"快度 φ₁ = {rapidity1:.6f} rad")
    print(f"快度 φ₂ = {rapidity2:.6f} rad")
    print(f"快度叠加 φ = φ₁ + φ₂ = {rapidity_sum:.6f} rad")

    # 从快度反推速度
    beta_from_rapidity = math.tanh(rapidity_sum)
    print(f"从快度反推: β = tanh(φ) = {beta_from_rapidity:.6f}")
    print(f"与相对论叠加结果一致: {beta_combined:.6f}")


def section_5_spacetime_intervals():
    """第五节：时空区间与光锥"""
    print("\n" + "=" * 60)
    print("第五节：时空区间与光锥")
    print("=" * 60)

    alg = nblade.Algebra.spacetime(4)
    e0, e1, e2, e3 = alg.basis_vectors()

    c = 1.0

    print("\n【时空区间分类】")
    print("时空区间 Δs² = (cΔt)² - Δx² - Δy² - Δz²")
    print("  • Δs² > 0: 类时间隔 (timelike)")
    print("  • Δs² < 0: 类空间隔 (spacelike)")
    print("  • Δs² = 0: 类光间隔 (lightlike/null)")

    # 创建不同类型的事件间隔
    print("\n【类时间隔示例】")
    # 时间间隔大于空间间隔
    dt = 2.0
    dx, dy, dz = 0.5, 0.0, 0.0
    X_timelike = e0.scale(dt * c) + e1.scale(dx) + e2.scale(dy) + e3.scale(dz)
    interval_t = (X_timelike | X_timelike).scalar_part()
    print(f"Δt = {dt}, Δx = {dx}, Δy = {dy}, Δz = {dz}")
    print(f"Δs² = {interval_t:.6f} > 0 → 类时间隔")
    print("存在参考系使两事件发生在同一位置")

    print("\n【类空间隔示例】")
    # 空间间隔大于时间间隔
    dt = 0.5
    dx, dy, dz = 2.0, 0.0, 0.0
    X_spacelike = e0.scale(dt * c) + e1.scale(dx) + e2.scale(dy) + e3.scale(dz)
    interval_s = (X_spacelike | X_spacelike).scalar_part()
    print(f"Δt = {dt}, Δx = {dx}, Δy = {dy}, Δz = {dz}")
    print(f"Δs² = {interval_s:.6f} < 0 → 类空间隔")
    print("存在参考系使两事件同时发生")
    print("因果次序可能颠倒！")

    print("\n【类光间隔示例】")
    # 光线连接的事件
    dt = 1.0
    dx, dy, dz = 1.0, 0.0, 0.0  # 空间距离等于 c*dt
    X_lightlike = e0.scale(dt * c) + e1.scale(dx) + e2.scale(dy) + e3.scale(dz)
    interval_l = (X_lightlike | X_lightlike).scalar_part()
    print(f"Δt = {dt}, Δx = {dx}, Δy = {dy}, Δz = {dz}")
    print(f"Δs² = {interval_l:.6f} ≈ 0 → 类光间隔")
    print("只有光可以连接这两个事件")

    print("\n【光锥结构】")
    print("光锥定义了因果结构:")
    print("  • 未来光锥内: 类时事件，可被影响")
    print("  • 过去光锥内: 类时事件，曾影响")
    print("  • 光锥面上: 类光事件，光信号可达")
    print("  • 光锥外: 类空事件，无因果关联")

    # 光锥边界的法向量
    print("\n【光锥边界】")
    print("光锥边界由类光向量构成")
    print("光锥方程: (ct)² - x² - y² - z² = 0")

    # 沿各方向的光
    light_x_pos = e0.scale(1.0) + e1.scale(1.0)
    light_x_neg = e0.scale(1.0) - e1.scale(1.0)
    light_y_pos = e0.scale(1.0) + e2.scale(1.0)
    light_y_neg = e0.scale(1.0) - e2.scale(1.0)

    print(f"\n沿 +x 方向的光: {light_x_pos}")
    print(f"  平方 = {(light_x_pos | light_x_pos).scalar_part():.6f} (类光)")

    print(f"\n沿 -x 方向的光: {light_x_neg}")
    print(f"  平方 = {(light_x_neg | light_x_neg).scalar_part():.6f} (类光)")


def section_6_special_relativity():
    """第六节：与狭义相对论的关系"""
    print("\n" + "=" * 60)
    print("第六节：时空代数与狭义相对论")
    print("=" * 60)

    alg = nblade.Algebra.spacetime(4)
    e0, e1, e2, e3 = alg.basis_vectors()

    c = 1.0

    print("\n【时空代数的优势】")
    print("1. 几何直观: 向量和张量统一为多向量")
    print("2. 计算简洁: 复杂的矩阵运算变为几何积")
    print("3. 坐标无关: 物理量直接表示，无需坐标")

    print("\n【麦克斯韦方程组的统一】")
    print("传统形式需要 4 个方程，时空代数中只需 1 个:")
    print("  ∇F = J/ε₀c")
    print("其中 F 是电磁场二重向量，J 是四电流")

    print("\n【狄拉克方程的几何解释】")
    print("狄拉克方程可以用时空代数简洁表示")
    print("揭示了电子自旋的几何本质")

    print("\n【时间膨胀验证】")
    # 运动时钟变慢
    beta = 0.8
    gamma = 1.0 / math.sqrt(1 - beta**2)

    print(f"\n以 v = {beta}c 运动的参考系:")
    print(f"洛伦兹因子 γ = {gamma:.6f}")
    print(f"时间膨胀: Δt' = γΔt = {gamma:.6f}Δt")
    print(f"运动的时钟看起来慢了 {gamma:.6f} 倍")

    # 四速度归一化验证
    v_x, v_y, v_z = beta, 0.0, 0.0
    U = (
        e0.scale(gamma * c)
        + e1.scale(gamma * v_x * c)
        + e2.scale(gamma * v_y * c)
        + e3.scale(gamma * v_z * c)
    )
    U_sq = (U | U).scalar_part()

    print(f"\n四速度 U = γ(c, vₓc, vᵧc, vᵤc)")
    print(f"U² = {U_sq:.6f} (恒等于 c² = {c * c})")

    print("\n【长度收缩验证】")
    print(f"运动方向长度收缩: L' = L/γ = L/{gamma:.6f}")
    print(f"静止的 {1.0 / gamma:.4f} 米尺在运动系中看起来有 1 米")

    print("\n【质能关系验证】")
    m0 = 1.0
    E = gamma * m0 * c**2
    p = gamma * m0 * beta * c

    P = e0.scale(E / c) + e1.scale(p)
    P_sq = (P | P).scalar_part()

    print(f"\n静止质量 m₀ = {m0}")
    print(f"总能量 E = γm₀c² = {E:.6f}")
    print(f"动能 T = E - m₀c² = {E - m0 * c * c:.6f}")
    print(f"动量 |p| = {p:.6f}")
    print(f"\nP² = E²/c² - |p|² = {P_sq:.6f}")
    print(f"预期 m₀²c² = {m0**2 * c**2:.6f}")
    print(f"验证 E² - |p|²c² = m₀²c⁴ ✓")

    print("\n【时空代数与传统方法的对比】")
    print("-" * 40)
    print("| 操作         | 传统方法      | 时空代数      |")
    print("|-------------|---------------|---------------|")
    print("| 洛伦兹变换   | 4×4 矩阵乘法  | 转子变换      |")
    print("| 电磁场      | E 和 B 分离   | F 二重向量    |")
    print("| 旋转+推进   | 矩阵分解      | 单一转子      |")
    print("| 不变量      | 多个方程      | 单一多向量    |")

    print("\n" + "=" * 60)
    print("时空代数将狭义相对论统一于几何框架")
    print("使物理关系更加清晰和优雅")
    print("=" * 60)


def main():
    """主函数"""
    print("=" * 60)
    print("nblade 时空代数示例")
    print("几何代数在狭义相对论中的应用")
    print("=" * 60)

    section_1_spacetime_algebra()
    section_2_four_vectors()
    section_3_lorentz_transformations()
    section_4_velocity_addition()
    section_5_spacetime_intervals()
    section_6_special_relativity()

    print("\n" + "=" * 60)
    print("示例完成！")
    print("=" * 60)
    print("\n关键点:")
    print("  - 时空代数 G(1,3,0): e₀²=+1, e₁²=e₂²=e₃²=-1")
    print("  - 四向量: 位置、速度、动量的统一描述")
    print("  - 洛伦兹变换: 用转子统一推进和旋转")
    print("  - 速度叠加: 遵循相对论公式")
    print("  - 时空区间: 分类因果结构")
    print("\n参考:")
    print("  - Hestenes, D. 'Space-Time Algebra' (1966)")
    print("  - Doran & Lasenby 'Geometric Algebra for Physicists'")


if __name__ == "__main__":
    main()
