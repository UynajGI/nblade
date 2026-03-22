"""
几何代数应用示例

演示几何代数在实际问题中的应用：
- 2D/3D 旋转
- 物理应用
- 几何变换
"""

from nblade import Algebra
import math


def demo_2d_complex():
    """演示 2D 几何代数与复数的等价性"""
    print("=== 2D 几何代数与复数 ===\n")

    # 2D 欧几里得代数 G(2,0,0)
    alg_2d = Algebra.euclidean(2)
    e1, e2 = alg_2d.basis_vectors()

    # 2D 伪标量相当于虚数单位 i
    I = e1 ^ e2  # I² = -1
    print(f"2D 伪标量 I = e1∧e2 = {I}")
    print(f"I² = {I * I}\n")

    # 复数表示
    z1 = alg_2d.from_scalar(3.0) + 4.0 * I  # 3 + 4i
    z2 = alg_2d.from_scalar(1.0) + 2.0 * I  # 1 + 2i

    print(f"复数 z1 = {z1}")
    print(f"复数 z2 = {z2}\n")

    # 复数运算
    print("复数运算:")
    print(f"  z1 + z2 = {z1 + z2}")
    print(f"  z1 * z2 = {z1 * z2}")

    # 共轭
    z1_conj = z1.grade_involution()  # 对于偶子代数，这相当于复共轭
    print(f"  z1* = {z1_conj}")

    # 模长平方
    modulus_sq = z1_conj * z1
    print(f"  |z1|² = z1*z1 = {modulus_sq}\n")

    # 旋转 (使用复数乘法)
    angle = math.pi / 4  # 45 度
    exp_i_theta = alg_2d.from_scalar(math.cos(angle)) + math.sin(angle) * I
    print(f"旋转 45° 算子: {exp_i_theta}")

    # 旋转向量 e1
    rotated_e1 = exp_i_theta * e1
    print(f"  e1 旋转 45°: {rotated_e1}")
    print(
        f"  与 e1*cos(π/4) + e2*sin(π/4) 比较: {math.cos(angle) * e1 + math.sin(angle) * e2}\n"
    )


def demo_3d_quaternions():
    """演示 3D 几何代数与四元数的等价性"""
    print("=== 3D 几何代数与四元数 ===\n")

    # 3D 欧几里得代数 G(3,0,0)
    alg_3d = Algebra.euclidean(3)
    e1, e2, e3 = alg_3d.basis_vectors()

    # 四元数单位 (二重向量)
    i = e2 ^ e3  # e2e3
    j = e3 ^ e1  # e3e1
    k = e1 ^ e2  # e1e2

    print("四元数单位 (作为二重向量):")
    print(f"  i = e2∧e3 = {i}")
    print(f"  j = e3∧e1 = {j}")
    print(f"  k = e1∧e2 = {k}\n")

    # 验证四元数乘法规则
    print("四元数乘法规则验证:")
    print(f"  i² = {i * i}")
    print(f"  j² = {j * j}")
    print(f"  k² = {k * k}")
    print(f"  ij = {i * j}")
    print(f"  jk = {j * k}")
    print(f"  ki = {k * i}")
    print(f"  ijk = {i * j * k}\n")

    # 创建四元数 (标量 + 二重向量)
    q = alg_3d.from_scalar(1.0) + 2 * i + 3 * j + 4 * k
    print(f"四元数 q = {q}")

    # 四元数共轭
    q_conj = q.grade_involution()  # 对于四元数（偶子代数），共轭就是阶次对合
    print(f"四元数共轭 q* = {q_conj}")

    # 范数平方
    norm_sq = q_conj * q
    print(f"范数平方 |q|² = {norm_sq}\n")


def demo_rotations():
    """演示旋转"""
    print("=== 旋转运算 ===\n")

    alg = Algebra.euclidean(3)
    e1, e2, e3 = alg.basis_vectors()

    # 创建旋转平面 (xy 平面)
    plane = e1 ^ e2
    print(f"旋转平面: {plane}")

    # 创建旋转角度的转子
    angle = math.pi / 3  # 60 度
    # 转子 R = exp(-Iθ/2) 其中 I 是平面的单位二重向量
    try:
        rotor = alg.rotor(plane, angle)
        print(f"60° 旋转转子: {rotor}")

        # 验证转子性质
        rotor_inv = rotor.inverse()
        print(f"转子逆: {rotor_inv}")
        print(f"RR† = {(rotor * rotor.reversion()).scalar_part():.6f} (应为 1)\n")

        # 旋转向量
        print("旋转向量:")
        original = e1
        rotated = rotor * original * rotor_inv
        print(f"  e1 旋转 60°: {rotated}")

        # 理论结果：e1 旋转 60° 应该是 cos(60°)e1 + sin(60°)e2
        expected = math.cos(angle) * e1 + math.sin(angle) * e2
        print(f"  理论结果: {expected}")
        error = (rotated - expected).norm()
        print(f"  误差: {error:.10f}\n")

        # 旋转二重向量
        bivector_to_rotate = e2 ^ e3
        rotated_biv = rotor * bivector_to_rotate * rotor_inv
        print(f"  (e2∧e3) 旋转 60°: {rotated_biv}\n")
    except AttributeError:
        print("rotor 方法尚未实现，跳过旋转测试\n")


def demo_physics_application():
    """物理应用：角动量作为二重向量"""
    print("=== 物理应用：角动量 ===\n")

    alg = Algebra.euclidean(3)
    e1, e2, e3 = alg.basis_vectors()

    # 位置向量
    r = alg.vector([3, 4, 0])
    print(f"位置向量 r = {r}")

    # 动量向量
    p = alg.vector([2, -1, 5])
    print(f"动量向量 p = {p}\n")

    # 角动量作为二重向量 L = r∧p
    angular_momentum = r ^ p
    print(f"角动量二重向量 L = r∧p = {angular_momentum}")

    # 解释角动量的分量
    # L = L₁₂(e₁∧e₂) + L₁₃(e₁∧e₃) + L₂₃(e₂∧e₃)
    l12 = angular_momentum.grade(2).project_to(e1 ^ e2).scalar_part()
    l13 = angular_momentum.grade(2).project_to(e1 ^ e3).scalar_part()
    l23 = angular_momentum.grade(2).project_to(e2 ^ e3).scalar_part()

    print(f"  L_xy (e1∧e2) 分量: {l12:.6f}")
    print(f"  L_xz (e1∧e3) 分量: {l13:.6f}")
    print(f"  L_yz (e2∧e3) 分量: {l23:.6f}\n")

    # 验证经典公式 L = r × p (通过 Hodge 对偶)
    # 在 3D 中，L = *(r∧p) 对应于经典叉积
    dual_l = angular_momentum.dual()
    print(f"角动量对偶 (对应经典叉积): {dual_l}\n")


def demo_geometric_transformations():
    """几何变换"""
    print("=== 几何变换 ===\n")

    alg = Algebra.euclidean(3)
    e1, e2, e3 = alg.basis_vectors()

    # 创建一个三角形 (三个点)
    p1 = e1
    p2 = e2
    p3 = e3
    print("原始三角形顶点:")
    print(f"  p1 = {p1}")
    print(f"  p2 = {p2}")
    print(f"  p3 = {p3}\n")

    # 创建三角形的边
    edge12 = p2 - p1
    edge13 = p3 - p1
    triangle_plane = edge12 ^ edge13
    print(f"三角形平面: {triangle_plane}")
    print(f"三角形面积: {triangle_plane.norm():.6f}\n")

    # 将三角形投影到 xy 平面
    xy_plane = e1 ^ e2
    try:
        proj_p1 = p1.project_to(xy_plane)
        proj_p2 = p2.project_to(xy_plane)
        proj_p3 = p3.project_to(xy_plane)

        print("投影到 xy 平面:")
        print(f"  投影 p1 = {proj_p1}")
        print(f"  投影 p2 = {proj_p2}")
        print(f"  投影 p3 = {proj_p3}\n")

        # 计算投影后的面积
        proj_edge12 = proj_p2 - proj_p1
        proj_edge13 = proj_p3 - proj_p1
        proj_triangle = proj_edge12 ^ proj_edge13
        print(f"投影后面积: {proj_triangle.norm():.6f}\n")
    except AttributeError:
        print("project_to 方法尚未实现，跳过投影测试\n")


def demo_cramers_rule():
    """演示 Cramer 法则"""
    print("=== Cramer 法则 ===\n")

    alg = Algebra.euclidean(3)
    e1, e2, e3 = alg.basis_vectors()

    # 求解线性方程组: x = αa + βb
    # 其中 x = [7, 8, 9], a = [1, 2, 3], b = [4, 5, 6]
    x = alg.vector([7, 8, 9])
    a = alg.vector([1, 2, 3])
    b = alg.vector([4, 5, 6])

    print(f"x = {x}")
    print(f"a = {a}")
    print(f"b = {b}\n")

    # 计算系数
    # β = (x∧a)/(b∧a)
    # α = (x∧b)/(a∧b)

    denom = b ^ a
    print(f"分母 b∧a = {denom}")

    if not denom.is_zero():
        beta_num = x ^ a
        alpha_num = x ^ b

        print(f"β 分子 x∧a = {beta_num}")
        print(f"α 分子 x∧b = {alpha_num}")

        # 由于分母是二重向量，我们需要使用逆
        try:
            beta = (beta_num * denom.inverse()).grade(1)  # 取向量部分
            alpha = (alpha_num * (-denom.inverse())).grade(1)  # 注意符号

            print(f"β = {beta}")
            print(f"α = {alpha}")

            # 验证解
            solution = alpha + beta
            print(f"验证: αa + βb = {solution}")
            print(f"误差: {((solution - x).norm_squared()) ** 0.5:.10f}\n")
        except Exception as e:
            print(f"计算逆时出错: {e}\n")
    else:
        print("分母为零，无法使用此方法\n")


if __name__ == "__main__":
    demo_2d_complex()
    demo_3d_quaternions()
    demo_rotations()
    demo_physics_application()
    demo_geometric_transformations()
    demo_cramers_rule()

    print("所有应用示例完成！")
