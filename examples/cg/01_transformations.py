"""
nblade 计算机图形学变换示例
===========================

演示几何代数在计算机图形学中的应用：
- 基本几何变换 (平移、旋转、缩放)
- 反射和投影
- 视图变换
- 齐次坐标与共形几何代数简介

作者: nblade 团队
"""

import nblade
import math


def section_1_reflections():
    """第一节：反射"""
    print("\n" + "=" * 60)
    print("第一节：反射")
    print("=" * 60)

    alg = nblade.Algebra.euclidean(3)
    e1, e2, e3 = alg.basis_vectors()

    print("\n反射是最基本的几何变换")
    print("公式: v' = -n·v·n (其中 n 是单位法向量)")

    # 向量反射
    print("\n【向量反射】")
    v = alg.vector([1.0, 2.0, 3.0])
    n = e1  # 反射面的法向量 (yz 平面)

    reflected = v.reflect_in(n)
    print(f"原始向量 v = {v}")
    print(f"反射面法向量 n = {n} (yz 平面)")
    print(f"反射结果 v' = {reflected}")
    print("  注意: x 分量取反")

    # 点反射
    print("\n【点关于平面的反射】")
    point = alg.vector([2.0, 3.0, 4.0])
    plane_normal = e2  # xz 平面

    reflected_point = point.reflect_in(plane_normal)
    print(f"点 P = {point}")
    print(f"反射面: xz 平面 (法向量 e2)")
    print(f"反射后 P' = {reflected_point}")

    # 多次反射
    print("\n【多次反射 = 旋转】")
    # 两次反射等价于一次旋转
    n1 = e1
    n2 = (e1 + e2).scale(1.0 / math.sqrt(2))  # 45° 方向

    v = e1
    # 先对 n1 反射，再对 n2 反射
    v_ref1 = v.reflect_in(n1)
    v_ref2 = v_ref1.reflect_in(n2)

    print(f"n1 = {n1}")
    print(f"n2 = {n2} (45° 方向)")
    print(f"v = {v}")
    print(f"第一次反射: {v_ref1}")
    print(f"第二次反射: {v_ref2}")
    print("  两次反射 = 绕两平面交线旋转")


def section_2_projections():
    """第二节：投影"""
    print("\n" + "=" * 60)
    print("第二节：投影")
    print("=" * 60)

    alg = nblade.Algebra.euclidean(3)
    e1, e2, e3 = alg.basis_vectors()

    print("\n投影将向量投影到子空间")

    # 向量投影
    print("\n【向量投影】")
    v = alg.vector([3.0, 4.0, 5.0])
    direction = e1  # 投影到 x 轴

    proj = v.project_to(direction)
    reject = v.reject_from(direction)

    print(f"v = {v}")
    print(f"投影方向: {direction}")
    print(f"投影分量 = {proj}")
    print(f"拒绝分量 = {reject}")
    print(f"验证: proj + reject = {proj + reject}")

    # 平面投影
    print("\n【平面投影】")
    plane = e1 ^ e2  # xy 平面
    v3d = alg.vector([2.0, 3.0, 5.0])

    proj_to_plane = v3d.project_to(plane)
    print(f"v = {v3d}")
    print(f"平面 = {plane} (xy 平面)")
    print(f"投影到平面 = {proj_to_plane}")
    print("  注意: z 分量被消除")

    # 斜投影
    print("\n【斜投影】")
    v = alg.vector([1.0, 1.0, 1.0])
    onto = e1 ^ e2  # xy 平面
    along = e3  # 沿 z 方向

    # 正交投影
    ortho_proj = v.project_to(onto)
    print(f"v = {v}")
    print(f"正交投影到 xy 平面 = {ortho_proj}")


def section_3_rotations():
    """第三节：旋转"""
    print("\n" + "=" * 60)
    print("第三节：旋转")
    print("=" * 60)

    alg = nblade.Algebra.euclidean(3)
    e1, e2, e3 = alg.basis_vectors()

    print("\n旋转使用转子实现，比矩阵更高效")

    # 绕坐标轴旋转
    print("\n【绕坐标轴旋转】")

    # 绕 z 轴旋转 90°
    angle = math.pi / 2
    rotor_z = alg.rotor(e1 ^ e2, angle)

    v = e1
    v_rotated = v.rotate_by(rotor_z)

    print(f"绕 z 轴旋转 90°")
    print(f"原始: v = {v}")
    print(f"旋转后: v' = {v_rotated}")

    # 绕任意轴旋转
    print("\n【绕任意轴旋转】")
    axis = alg.vector([1.0, 1.0, 1.0])
    axis_normalized = axis.scale(1.0 / axis.norm())

    # 旋转平面垂直于旋转轴
    # 对于任意轴，旋转平面由轴决定

    angle = math.pi / 3  # 60°

    print(f"旋转轴: {axis_normalized}")
    print(f"旋转角度: {math.degrees(angle):.1f}°")

    # 使用 xy 平面旋转作为示例
    rotor_custom = alg.rotor(e1 ^ e2, angle)
    v = alg.vector([1.0, 0.0, 0.0])
    v_rot = v.rotate_by(rotor_custom)
    print(f"v = {v} 旋转后 = {v_rot}")

    # 欧拉角
    print("\n【欧拉角】")
    # 绕 z 轴 → 绕 y 轴 → 绕 x 轴
    alpha = math.pi / 4  # 绕 z
    beta = math.pi / 6  # 绕 y
    gamma = math.pi / 3  # 绕 x

    R_z = alg.rotor(e1 ^ e2, alpha)
    R_y = alg.rotor(e2 ^ e3, beta)
    R_x = alg.rotor(e3 ^ e1, gamma)

    # 组合旋转 (注意顺序)
    R_total = R_x * R_y * R_z

    print(
        f"欧拉角: α={math.degrees(alpha):.0f}°, β={math.degrees(beta):.0f}°, γ={math.degrees(gamma):.0f}°"
    )
    print(f"组合转子 R = R_x·R_y·R_z")

    v = e1
    v_rotated = v.rotate_by(R_total)
    print(f"v = {v}")
    print(f"旋转后 = {v_rotated}")


def section_4_scaling():
    """第四节：缩放"""
    print("\n" + "=" * 60)
    print("第四节：缩放")
    print("=" * 60)

    alg = nblade.Algebra.euclidean(3)
    e1, e2, e3 = alg.basis_vectors()

    print("\n缩放是最简单的变换")

    # 均匀缩放
    print("\n【均匀缩放】")
    v = alg.vector([1.0, 2.0, 3.0])
    scale_factor = 2.0

    v_scaled = v.scale(scale_factor)
    print(f"v = {v}")
    print(f"缩放因子: {scale_factor}")
    print(f"缩放后: {v_scaled}")

    # 非均匀缩放
    print("\n【非均匀缩放】")
    # 使用分量乘法
    v = alg.vector([1.0, 2.0, 3.0])
    v_scaled = (
        e1.scale(v.coefficients()[1] * 2.0)
        + e2.scale(v.coefficients()[2] * 3.0)
        + e3.scale(v.coefficients()[3] * 4.0)
    )

    print(f"原始: v = {v}")
    print(f"非均匀缩放 (x2, x3, x4): {v_scaled}")


def section_5_combined_transformations():
    """第五节：组合变换"""
    print("\n" + "=" * 60)
    print("第五节：组合变换")
    print("=" * 60)

    alg = nblade.Algebra.euclidean(3)
    e1, e2, e3 = alg.basis_vectors()

    print("\n多个变换可以通过几何积组合")

    # 变换序列
    print("\n【变换序列示例】")

    # 1. 缩放
    scale = 2.0
    v = alg.vector([1.0, 1.0, 0.0])
    v_scaled = v.scale(scale)
    print(f"1. 原始向量: {v}")
    print(f"   缩放 x{scale}: {v_scaled}")

    # 2. 旋转
    angle = math.pi / 4
    rotor = alg.rotor(e1 ^ e2, angle)
    v_rotated = v_scaled.rotate_by(rotor)
    print(f"\n2. 旋转 45°: {v_rotated}")

    # 3. 反射
    reflected = v_rotated.reflect_in(e1)
    print(f"\n3. 反射 (yz 平面): {reflected}")

    # 验证顺序
    print("\n【变换顺序的重要性】")
    v = e1

    # 先旋转后反射
    v_rr = v.rotate_by(rotor).reflect_in(e1)
    print(f"先旋转后反射: {v_rr}")

    # 先反射后旋转
    v_ref = v.reflect_in(e1)
    v_rr2 = v_ref.rotate_by(rotor)
    print(f"先反射后旋转: {v_rr2}")
    print("  顺序不同，结果不同！")


def section_6_triangle_transformations():
    """第六节：三角形变换示例"""
    print("\n" + "=" * 60)
    print("第六节：三角形变换示例")
    print("=" * 60)

    alg = nblade.Algebra.euclidean(3)
    e1, e2, e3 = alg.basis_vectors()

    print("\n对一个三角形应用各种变换")

    # 定义三角形顶点
    p1 = e1
    p2 = e2
    p3 = e3 * 0.5

    print(f"\n原始三角形:")
    print(f"  p1 = {p1}")
    print(f"  p2 = {p2}")
    print(f"  p3 = {p3}")

    # 计算三角形法向量
    edge1 = p2 - p1
    edge2 = p3 - p1
    normal = (edge1 ^ edge2).dual()
    print(f"\n法向量: {normal.scale(1.0 / normal.norm())}")

    # 旋转变换
    print("\n【旋转变换】")
    rotor = alg.rotor(e1 ^ e2, math.pi / 2)

    p1_rot = p1.rotate_by(rotor)
    p2_rot = p2.rotate_by(rotor)
    p3_rot = p3.rotate_by(rotor)

    print(f"旋转 90° 后:")
    print(f"  p1' = {p1_rot}")
    print(f"  p2' = {p2_rot}")
    print(f"  p3' = {p3_rot}")

    # 计算旋转后的法向量
    edge1_rot = p2_rot - p1_rot
    edge2_rot = p3_rot - p1_rot
    normal_rot = (edge1_rot ^ edge2_rot).dual()
    print(f"旋转后法向量: {normal_rot.scale(1.0 / normal_rot.norm())}")

    # 镜像变换
    print("\n【镜像变换】")
    mirror_plane = e1

    p1_mir = p1.reflect_in(mirror_plane)
    p2_mir = p2.reflect_in(mirror_plane)
    p3_mir = p3.reflect_in(mirror_plane)

    print(f"镜像 (yz 平面) 后:")
    print(f"  p1'' = {p1_mir}")
    print(f"  p2'' = {p2_mir}")
    print(f"  p3'' = {p3_mir}")


def section_7_conformal_ga_intro():
    """第七节：共形几何代数简介"""
    print("\n" + "=" * 60)
    print("第七节：共形几何代数简介")
    print("=" * 60)

    print("\n共形几何代数 (CGA) 可以统一表示:")
    print("  - 点、线、平面")
    print("  - 圆、球")
    print("  - 刚体运动 (包括平移)")

    # 创建 CGA
    try:
        alg_cga = nblade.Algebra.cga()  # G(4,1,0)
        print(f"\n共形几何代数: {alg_cga}")
        print("  5 维空间，其中 4 个正平方基向量，1 个负平方")

        basis = alg_cga.basis_vectors()
        print(f"  基向量数量: {len(basis)}")

        print("\nCGA 中的基本概念:")
        print("  - 原点: e₋ (负平方基向量)")
        print("  - 无穷远点: e₊ (正平方基向量之一)")
        print("  - 点: P = e₋ + x + x²e₊/2")
        print("  - 平面: π = n + δe₋ (n 是法向量，δ 是距离)")
        print("  - 球: s = P - r²e₊/2 (P 是中心，r 是半径)")

    except Exception as e:
        print(f"\nCGA 示例需要完整的 5D 支持: {e}")

    print("\nCGA 的优势:")
    print("  1. 所有几何对象都是多向量")
    print("  2. 交运算通过外积实现")
    print("  3. 平移也是转子，与旋转统一")
    print("  4. 公式简洁，避免特殊情况")


def main():
    """主函数"""
    print("=" * 60)
    print("nblade 计算机图形学变换示例")
    print("=" * 60)

    section_1_reflections()
    section_2_projections()
    section_3_rotations()
    section_4_scaling()
    section_5_combined_transformations()
    section_6_triangle_transformations()
    section_7_conformal_ga_intro()

    print("\n" + "=" * 60)
    print("示例完成！")
    print("=" * 60)
    print("\n关键点:")
    print("  - 反射是最基本变换: v' = -nvn")
    print("  - 旋转使用转子: v' = RvR†")
    print("  - 投影和拒绝分解向量")
    print("  - 多个变换通过几何积组合")
    print("  - 共形几何代数统一表示所有几何对象")


if __name__ == "__main__":
    main()
