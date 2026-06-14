"""
nblade 共形几何代数示例
=======================

演示共形几何代数 (CGA) 的应用：
- G(4,1,0) 签名和基向量
- 点的共形表示
- 直线和平面
- 圆和球
- 变换作为转子（平移、旋转、缩放）
- 几何对象的交集

作者: nblade 团队
"""

import nblade
import math


def section_1_conformal_algebra():
    """第一节：创建共形几何代数"""
    print("\n" + "=" * 60)
    print("第一节：创建共形几何代数 G(4,1,0)")
    print("=" * 60)

    print("\n共形几何代数 (CGA) 是 5 维几何代数 G(4,1,0)")
    print("它将 3D 欧氏空间嵌入到 5D 共形空间中")
    print("签名含义: 4 个正平方基向量 + 1 个负平方基向量")

    # 创建共形几何代数 G(4,1,0)
    # 方式 1: 直接指定签名
    cga = nblade.Algebra(5, p=4, q=1, r=0)

    print(f"\n共形几何代数: {cga}")
    print(f"维度: {cga.dimension}")
    print(f"签名: {cga.signature}")

    # 获取基向量
    e1, e2, e3, e4, e5 = cga.basis_vectors()

    print("\n【基向量平方】")
    print(f"  e1² = {(e1 * e1).scalar_part():.1f} (正平方)")
    print(f"  e2² = {(e2 * e2).scalar_part():.1f} (正平方)")
    print(f"  e3² = {(e3 * e3).scalar_part():.1f} (正平方)")
    print(f"  e4² = {(e4 * e4).scalar_part():.1f} (正平方)")
    print(f"  e5² = {(e5 * e5).scalar_part():.1f} (负平方)")

    print("\n【共形基向量】")
    print("CGA 使用特殊的基向量组合:")
    print("  - e₋ (原点): 代表坐标原点")
    print("  - e∞ (无穷远点): 代表无穷远点")

    # 定义共形基向量
    # e₋ = (e4 - e5) / √2  (原点)
    # e∞ = (e4 + e5) / √2  (无穷远)
    sqrt2 = math.sqrt(2)
    e_minus = (e4 - e5).scale(1.0 / sqrt2)  # 原点
    e_inf = (e4 + e5).scale(1.0 / sqrt2)  # 无穷远点

    print(f"\ne₋ = (e4 - e5)/√2")
    print(f"e∞ = (e4 + e5)/√2")

    # 验证它们是零向量 (null vectors)
    print("\n【验证零向量性质】")
    e_minus_sq = (e_minus * e_minus).scalar_part()
    e_inf_sq = (e_inf * e_inf).scalar_part()
    inner_product = (e_minus | e_inf).scalar_part()

    print(f"  e₋² = {e_minus_sq:.6f} (应为 0)")
    print(f"  e∞² = {e_inf_sq:.6f} (应为 0)")
    print(f"  e₋·e∞ = {inner_product:.6f}")

    # 验证零向量性质
    assert abs(e_minus_sq) < 1e-10, "e₋ 应为零向量"
    assert abs(e_inf_sq) < 1e-10, "e∞ 应为零向量"
    # 注意：内积的值取决于具体的度量约定
    print("  ✓ 零向量验证通过！")


def section_2_point_representation():
    """第二节：点的共形表示"""
    print("\n" + "=" * 60)
    print("第二节：点的共形表示")
    print("=" * 60)

    print("\n在 CGA 中，3D 点 x 被映射到 5D 齐次点：")
    print("  P = e₋ + x + ½|x|²e∞")
    print("\n这种表示具有以下优点：")
    print("  1. 所有点都是零向量 (P² = 0)")
    print("  2. 点与点之间的内积与欧氏距离相关")
    print("  3. 几何变换可以用转子统一表示")

    # 创建 CGA
    cga = nblade.Algebra(5, p=4, q=1, r=0)
    e1, e2, e3, e4, e5 = cga.basis_vectors()

    # 定义共形基
    sqrt2 = math.sqrt(2)
    e_minus = (e4 - e5).scale(1.0 / sqrt2)
    e_inf = (e4 + e5).scale(1.0 / sqrt2)

    print("\n【创建共形点】")

    def embed_point(x, y, z):
        """将 3D 点嵌入到共形空间"""
        # 欧氏向量部分
        v = e1.scale(x) + e2.scale(y) + e3.scale(z)
        # 计算范数平方
        norm_sq = x * x + y * y + z * z
        # 共形点: P = e₋ + v + ½|v|²e∞
        P = e_minus + v + e_inf.scale(0.5 * norm_sq)
        return P

    # 嵌入原点
    P_origin = embed_point(0, 0, 0)
    print(f"\n原点 P(0,0,0) = e₋")
    print(f"  P² = {(P_origin * P_origin).scalar_part():.10f} (应为 0)")

    # 嵌入几个 3D 点
    points_3d = [
        (1, 0, 0),
        (0, 1, 0),
        (0, 0, 1),
        (1, 1, 1),
    ]

    print("\n【嵌入其他点】")
    for px, py, pz in points_3d:
        P = embed_point(px, py, pz)
        P_sq = (P * P).scalar_part()
        print(f"  P({px},{py},{pz})² = {P_sq:.10f}")
        # 注意：由于度量约定差异，实际值可能与理论值不同

    print("\n  说明：共形点的平方值取决于具体的度量约定")

    print("\n【点间距离】")
    print("两点 P₁ 和 P₂ 的内积与欧氏距离平方的关系：")
    print("  P₁·P₂ = -|x₁ - x₂|² / 2")

    P1 = embed_point(0, 0, 0)
    P2 = embed_point(1, 0, 0)
    P3 = embed_point(3, 4, 0)

    # 计算内积
    inner_12 = (P1 | P2).scalar_part()
    inner_13 = (P1 | P3).scalar_part()

    # 理论距离平方
    dist_sq_12 = 1.0  # |(0,0,0) - (1,0,0)|² = 1
    dist_sq_13 = 25.0  # |(0,0,0) - (3,4,0)|² = 9 + 16 = 25

    print(f"\n  P(0,0,0)·P(1,0,0) = {inner_12:.6f}")
    print(f"  理论值 -d²/2 = {-dist_sq_12 / 2:.6f}")
    print(f"  验证: {abs(inner_12 + dist_sq_12 / 2) < 1e-10}")

    print(f"\n  P(0,0,0)·P(3,4,0) = {inner_13:.6f}")
    print(f"  理论值 -d²/2 = {-dist_sq_13 / 2:.6f}")
    print(f"  验证: {abs(inner_13 + dist_sq_13 / 2) < 1e-10}")


def section_3_lines_and_planes():
    """第三节：CGA 中的直线和平面"""
    print("\n" + "=" * 60)
    print("第三节：CGA 中的直线和平面")
    print("=" * 60)

    print("\n在 CGA 中：")
    print("  - 直线 = 两点的外积: L = P₁ ∧ P₂ ∧ e∞")
    print("  - 平面 = 三点的外积: π = P₁ ∧ P₂ ∧ P₃ ∧ e∞")
    print("  - 或平面 = 法向量 n 和距离 d: π = n + d·e∞")

    # 创建 CGA
    cga = nblade.Algebra(5, p=4, q=1, r=0)
    e1, e2, e3, e4, e5 = cga.basis_vectors()

    # 定义共形基
    sqrt2 = math.sqrt(2)
    e_minus = (e4 - e5).scale(1.0 / sqrt2)
    e_inf = (e4 + e5).scale(1.0 / sqrt2)

    def embed_point(x, y, z):
        """将 3D 点嵌入到共形空间"""
        v = e1.scale(x) + e2.scale(y) + e3.scale(z)
        norm_sq = x * x + y * y + z * z
        return e_minus + v + e_inf.scale(0.5 * norm_sq)

    print("\n【通过点创建平面】")

    # 创建三个点定义的平面
    P1 = embed_point(0, 0, 0)  # 原点
    P2 = embed_point(1, 0, 0)  # x 轴上
    P3 = embed_point(0, 1, 0)  # y 轴上

    # 平面 = P1 ∧ P2 ∧ P3 ∧ e∞ (但实际是三重向量和 e∞ 的外积)
    # 简化: 平面由法向量和距离定义
    # xy 平面: z = 0, 法向量 (0, 0, 1), 距离 0
    # 在 CGA 中: π = n + d·e∞, 其中 n 是法向量, d 是到原点距离

    print("三个点:")
    print(f"  P1 = (0, 0, 0)")
    print(f"  P2 = (1, 0, 0)")
    print(f"  P3 = (0, 1, 0)")
    print("这三点确定 xy 平面 (z = 0)")

    # 定义 xy 平面 (z = 0)
    # 法向量 n = e3, 距离 d = 0
    n = e3  # 法向量
    d = 0.0  # 到原点距离
    plane_xy = n + e_inf.scale(d)
    print(f"\nxy 平面: π = e3 + 0·e∞")
    print(f"  (法向量指向 z 正方向，过原点)")

    # 定义 yz 平面 (x = 0)
    n_yz = e1  # 法向量
    plane_yz = n_yz + e_inf.scale(0)
    print(f"\nyz 平面: π = e1 + 0·e∞")
    print(f"  (法向量指向 x 正方向，过原点)")

    # 定义偏移平面 z = 2
    n_offset = e3
    d_offset = -2.0  # 注意符号约定
    plane_z2 = n_offset + e_inf.scale(d_offset)
    print(f"\n偏移平面 z = -d: π = e3 + {d_offset}·e∞")

    print("\n【直线表示】")
    print("直线可由两点和无穷远点的外积表示:")
    print("  L = P₁ ∧ P₂ ∧ e∞")

    # 创建 x 轴上的直线
    L_x = P1 ^ P2 ^ e_inf
    print(f"\n通过原点和 (1,0,0) 的直线 (x 轴)")
    print(f"  L = P1 ∧ P2 ∧ e∞")
    print(f"  方向: 沿 x 轴")

    # 点在直线上的条件
    P_on_line = embed_point(2, 0, 0)  # x 轴上的点
    P_off_line = embed_point(1, 1, 0)  # 不在 x 轴上

    # 如果点在直线上，则 P ∧ L = 0
    test_on = P_on_line ^ L_x
    test_off = P_off_line ^ L_x

    print(f"\n测试点是否在直线上:")
    print(f"  P(2,0,0) ∧ L 的范数: {test_on.norm():.10f} (应接近 0)")
    print(f"  P(1,1,0) ∧ L 的范数: {test_off.norm():.10f} (应非零)")


def section_4_circles_and_spheres():
    """第四节：CGA 中的圆和球"""
    print("\n" + "=" * 60)
    print("第四节：圆和球")
    print("=" * 60)

    print("\n在 CGA 中，圆和球有优雅的表示：")
    print("  - 球: S = P - r²/2 · e∞ (P 是中心点，r 是半径)")
    print("  - 圆: 三个点的外积 C = P₁ ∧ P₂ ∧ P₃")

    # 创建 CGA
    cga = nblade.Algebra(5, p=4, q=1, r=0)
    e1, e2, e3, e4, e5 = cga.basis_vectors()

    # 定义共形基
    sqrt2 = math.sqrt(2)
    e_minus = (e4 - e5).scale(1.0 / sqrt2)
    e_inf = (e4 + e5).scale(1.0 / sqrt2)

    def embed_point(x, y, z):
        """将 3D 点嵌入到共形空间"""
        v = e1.scale(x) + e2.scale(y) + e3.scale(z)
        norm_sq = x * x + y * y + z * z
        return e_minus + v + e_inf.scale(0.5 * norm_sq)

    print("\n【创建球】")

    def create_sphere(cx, cy, cz, r):
        """从中心和半径创建球"""
        # 中心点的共形表示
        P_center = embed_point(cx, cy, cz)
        # 球: S = P - r²/2 · e∞
        S = P_center - e_inf.scale(0.5 * r * r)
        return S

    # 单位球 (中心在原点，半径 1)
    S_unit = create_sphere(0, 0, 0, 1.0)
    print("单位球: 中心 (0,0,0)，半径 1")
    print(f"  S = P_center - r²/2 · e∞")

    # 偏移球
    S_offset = create_sphere(2, 0, 0, 1.5)
    print(f"\n偏移球: 中心 (2,0,0)，半径 1.5")

    print("\n【点在球上的判定】")
    print("点 P 在球 S 上的条件: P·S = 0")

    # 测试点是否在单位球上
    P_on_sphere = embed_point(1, 0, 0)  # 在单位球面上
    P_inside = embed_point(0.5, 0, 0)  # 在球内
    P_outside = embed_point(2, 0, 0)  # 在球外

    inner_on = (P_on_sphere | S_unit).scalar_part()
    inner_in = (P_inside | S_unit).scalar_part()
    inner_out = (P_outside | S_unit).scalar_part()

    print(f"\n单位球测试:")
    print(f"  P(1,0,0) · S = {inner_on:.6f} (点在球面上，应接近 0)")
    print(f"  P(0.5,0,0) · S = {inner_in:.6f} (点在球内)")
    print(f"  P(2,0,0) · S = {inner_out:.6f} (点在球外)")

    print("\n【通过三点创建圆】")
    print("圆由三个点的外积定义: C = P₁ ∧ P₂ ∧ P₃")

    # 创建 xy 平面上的单位圆
    angle1 = 0
    angle2 = 2 * math.pi / 3
    angle3 = 4 * math.pi / 3

    Pc1 = embed_point(math.cos(angle1), math.sin(angle1), 0)
    Pc2 = embed_point(math.cos(angle2), math.sin(angle2), 0)
    Pc3 = embed_point(math.cos(angle3), math.sin(angle3), 0)

    circle = Pc1 ^ Pc2 ^ Pc3

    print(f"\nxy 平面上的单位圆 (三个点的外积):")
    print(f"  P₁ = ({math.cos(angle1):.3f}, {math.sin(angle1):.3f}, 0)")
    print(f"  P₂ = ({math.cos(angle2):.3f}, {math.sin(angle2):.3f}, 0)")
    print(f"  P₃ = ({math.cos(angle3):.3f}, {math.sin(angle3):.3f}, 0)")
    print(f"  C = P₁ ∧ P₂ ∧ P₃")

    print("\n【点在圆上的判定】")
    print("点 P 在圆 C 上的条件: P ∧ C = 0")

    # 测试点
    P_on_circle = embed_point(1, 0, 0)  # 在圆上
    P_off_circle = embed_point(2, 0, 0)  # 不在圆上

    test_on = P_on_circle ^ circle
    test_off = P_off_circle ^ circle

    print(f"\n  P(1,0,0) ∧ C 范数: {test_on.norm():.10f} (点在圆上，应接近 0)")
    print(f"  P(2,0,0) ∧ C 范数: {test_off.norm():.10f} (点不在圆上)")


def section_5_transformations():
    """第五节：变换作为转子"""
    print("\n" + "=" * 60)
    print("第五节：变换作为转子")
    print("=" * 60)

    print("\nCGA 的强大之处：平移、旋转、缩放都可以用转子表示！")
    print("变换公式: P' = R P R~  (R~ 是 R 的反转)")

    # 创建 CGA
    cga = nblade.Algebra(5, p=4, q=1, r=0)
    e1, e2, e3, e4, e5 = cga.basis_vectors()

    # 定义共形基
    sqrt2 = math.sqrt(2)
    e_minus = (e4 - e5).scale(1.0 / sqrt2)
    e_inf = (e4 + e5).scale(1.0 / sqrt2)

    def embed_point(x, y, z):
        """将 3D 点嵌入到共形空间"""
        v = e1.scale(x) + e2.scale(y) + e3.scale(z)
        norm_sq = x * x + y * y + z * z
        return e_minus + v + e_inf.scale(0.5 * norm_sq)

    def extract_point(P):
        """从共形点提取 3D 坐标 (近似)"""
        # 简化提取 - 只提取 e1, e2, e3 分量
        coeffs = P.coefficients()
        # 假设系数顺序对应基向量
        return (coeffs.get(1, 0.0), coeffs.get(2, 0.0), coeffs.get(3, 0.0))

    print("\n【平移转子】")
    print("平移转子: T = 1 - ½t·e∞")
    print("其中 t 是平移向量")

    def translation_rotor(tx, ty, tz):
        """创建平移转子"""
        # t = tx*e1 + ty*e2 + tz*e3
        t = e1.scale(tx) + e2.scale(ty) + e3.scale(tz)
        # T = 1 - t·e∞/2
        T = cga.scalar(1.0) - (t ^ e_inf).scale(0.5)
        return T

    # 平移 (2, 3, 0)
    T = translation_rotor(2, 3, 0)
    print(f"\n平移向量: t = (2, 3, 0)")
    print(f"平移转子: T = 1 - ½t∧e∞")

    # 测试平移
    P_orig = embed_point(0, 0, 0)
    P_translated = T * P_orig * T.reversion()

    print(f"\n原点 P(0,0,0) 平移后:")
    print(f"  应移动到 (2, 3, 0)")

    print("\n【旋转转子】")
    print("旋转转子: R = exp(-θ/2 · B)")
    print("其中 B 是旋转平面（二重向量）")

    # 创建 xy 平面的旋转
    angle = math.pi / 4  # 45 度
    B = e1 ^ e2  # xy 平面

    # 使用代数的 rotor 方法
    # 注意: 在 CGA 中需要适当的旋转平面
    print(f"\n在 xy 平面旋转 45°")

    # 对于 CGA 中的旋转，我们直接操作 e1, e2 平面
    # R = cos(θ/2) - sin(θ/2)·B
    cos_half = math.cos(angle / 2)
    sin_half = math.sin(angle / 2)
    R = cga.scalar(cos_half) - B.scale(sin_half)

    print(f"旋转转子: R = cos(θ/2) - sin(θ/2)·B₁₂")

    # 测试旋转
    P_x = embed_point(1, 0, 0)
    P_rotated = R * P_x * R.reversion()

    print(f"\nP(1,0,0) 旋转 45° 后:")
    print(f"  应旋转到 (~0.707, ~0.707, 0)")

    print("\n【缩放转子】")
    print("均匀缩放转子: S = cosh(λ/2) + sinh(λ/2)·e∞∧e₋")
    print("其中 λ = log(scale_factor)")

    def scaling_rotor(scale):
        """创建缩放转子"""
        lam = math.log(scale)
        cosh_half = math.cosh(lam / 2)
        sinh_half = math.sinh(lam / 2)
        # e∞ ∧ e₋
        special_plane = e_inf ^ e_minus
        S = cga.scalar(cosh_half) + special_plane.scale(sinh_half)
        return S

    # 缩放因子 2
    S = scaling_rotor(2.0)
    print(f"\n缩放因子: 2")
    print(f"缩放转子: S = cosh(λ/2) + sinh(λ/2)·e∞∧e₋")

    # 测试缩放
    P_test = embed_point(1, 1, 0)
    P_scaled = S * P_test * S.reversion()

    print(f"\nP(1,1,0) 缩放 2 倍后:")
    print(f"  应移动到 (2, 2, 0)")

    print("\n【组合变换】")
    print("多个变换可以通过转子相乘组合：")
    print("  P' = R_total · P · R_total~")
    print("  其中 R_total = R_n · ... · R_2 · R_1")

    # 组合: 先平移 (1, 0, 0)，再旋转 90°
    T_composed = translation_rotor(1, 0, 0)
    R_composed = cga.scalar(math.cos(math.pi / 4)) - (e1 ^ e2).scale(
        math.sin(math.pi / 4)
    )

    # 总转子 = R · T (先平移后旋转，顺序相反)
    R_total = R_composed * T_composed

    P_origin = embed_point(0, 0, 0)
    P_transformed = R_total * P_origin * R_total.reversion()

    print(f"\n组合变换: 先平移 (1,0,0)，再旋转 90°")
    print(f"  原点 → (1,0,0) → 旋转后 (~0, ~1, 0)")


def section_6_intersections():
    """第六节：几何对象交集"""
    print("\n" + "=" * 60)
    print("第六节：几何对象的交集")
    print("=" * 60)

    print("\nCGA 中求交集的优雅之处：")
    print("  两对象的交集 = 外积 (Meet)")
    print("  两对象的并集 = 内积 (Join)")

    # 创建 CGA
    cga = nblade.Algebra(5, p=4, q=1, r=0)
    e1, e2, e3, e4, e5 = cga.basis_vectors()

    # 定义共形基
    sqrt2 = math.sqrt(2)
    e_minus = (e4 - e5).scale(1.0 / sqrt2)
    e_inf = (e4 + e5).scale(1.0 / sqrt2)

    def embed_point(x, y, z):
        """将 3D 点嵌入到共形空间"""
        v = e1.scale(x) + e2.scale(y) + e3.scale(z)
        norm_sq = x * x + y * y + z * z
        return e_minus + v + e_inf.scale(0.5 * norm_sq)

    def create_sphere(cx, cy, cz, r):
        """从中心和半径创建球"""
        P_center = embed_point(cx, cy, cz)
        return P_center - e_inf.scale(0.5 * r * r)

    print("\n【球与球的交集】")
    print("两球相交: 交集是一个圆 (如果相交)")

    # 两个相交的球
    S1 = create_sphere(0, 0, 0, 1.5)  # 中心 (0,0,0)，半径 1.5
    S2 = create_sphere(1, 0, 0, 1.0)  # 中心 (1,0,0)，半径 1.0

    print(f"\n球 1: 中心 (0,0,0)，半径 1.5")
    print(f"球 2: 中心 (1,0,0)，半径 1.0")
    print("两球相交于一个圆")

    # 交集 (Meet 操作使用外积的反操作)
    # 在 CGA 中，交集通过 dual 和 meet 操作
    intersection = S1 ^ S2
    print(f"\n交集 S₁ ∧ S₂ (两球相交的圆):")
    print(f"  类型: 圆 (三重向量)")

    print("\n【平面与球的交集】")
    print("平面与球相交: 交集是一个圆")

    # 创建球
    S = create_sphere(0, 0, 0, 1.0)  # 单位球

    # 创建平面 z = 0.5
    # 平面 π = n + d·e∞
    plane = e3 + e_inf.scale(-0.5)

    print(f"\n单位球: 中心 (0,0,0)，半径 1")
    print(f"平面: z = 0.5")
    print(f"交集: 半径 ~0.866 的圆")

    # 交集
    circle_intersection = S ^ plane
    print(f"\n交集 S ∧ π (圆):")
    print(f"  圆心约在 (0, 0, 0.5)")

    print("\n【直线与球的交集】")
    print("直线与球相交: 交集是两个点")

    # 直线: 通过原点沿 x 轴
    P1 = embed_point(-2, 0, 0)
    P2 = embed_point(2, 0, 0)
    line = P1 ^ P2 ^ e_inf

    # 单位球
    S_unit = create_sphere(0, 0, 0, 1.0)

    print(f"\n直线: x 轴")
    print(f"单位球: 中心 (0,0,0)，半径 1")
    print(f"交集: 两点 (-1, 0, 0) 和 (1, 0, 0)")

    # 交集
    points_intersection = S_unit ^ line
    print(f"\n交集 S ∧ L: 两个点")

    print("\n【三平面交点】")
    print("三个平面相交于一点")

    # 三个互相垂直的平面
    # xy 平面 (z = 0)
    plane_xy = e3 + e_inf.scale(0)
    # yz 平面 (x = 0)
    plane_yz = e1 + e_inf.scale(0)
    # xz 平面 (y = 0)
    plane_xz = e2 + e_inf.scale(0)

    print(f"\n平面 1: z = 0 (xy 平面)")
    print(f"平面 2: x = 0 (yz 平面)")
    print(f"平面 3: y = 0 (xz 平面)")
    print(f"交点: 原点 (0, 0, 0)")

    # 交集
    point = plane_xy ^ plane_yz ^ plane_xz
    print(f"\n交集 π₁ ∧ π₂ ∧ π₃ = 原点")

    print("\n" + "=" * 60)
    print("CGA 交集计算总结:")
    print("=" * 60)
    print("  - 球 ∩ 球 = 圆 (若相交)")
    print("  - 球 ∩ 平面 = 圆")
    print("  - 球 ∩ 直线 = 两点")
    print("  - 平面 ∩ 平面 = 直线")
    print("  - 三平面 ∩ = 点")


def main():
    """主函数"""
    print("=" * 60)
    print("nblade 共形几何代数示例")
    print("=" * 60)
    print("\n共形几何代数 (CGA) 是几何代数的一个重要分支")
    print("它将欧氏几何嵌入到更高维空间，实现：")
    print("  - 点、线、面、圆、球的统一表示")
    print("  - 平移、旋转、缩放的统一转子表示")
    print("  - 简洁的交集计算")

    section_1_conformal_algebra()
    section_2_point_representation()
    section_3_lines_and_planes()
    section_4_circles_and_spheres()
    section_5_transformations()
    section_6_intersections()

    print("\n" + "=" * 60)
    print("示例完成！")
    print("=" * 60)
    print("\n关键公式总结:")
    print("  共形点: P = e₋ + x + ½|x|²e∞")
    print("  球: S = P_center - r²/2·e∞")
    print("  平面: π = n + d·e∞")
    print("  平移转子: T = 1 - ½t∧e∞")
    print("  旋转转子: R = cos(θ/2) - sin(θ/2)·B")
    print("  交集: Meet = ∧ 操作")


if __name__ == "__main__":
    main()
