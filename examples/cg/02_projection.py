"""
nblade 投影操作示例
==================

演示几何代数中的投影操作：
- 向量投影到另一个向量
- 子空间投影（投影到平面）
- 拒绝分量（正交分量）
- 实际应用：阴影计算

投影公式：proj_B(A) = (A⌋B)⌋B⁻¹
拒绝公式：rej_B(A) = A - proj_B(A)

作者: nblade 团队
"""

import nblade
import math


def section_1_vector_projection():
    """第一节：向量投影"""
    print("\n" + "=" * 60)
    print("第一节：向量投影")
    print("=" * 60)

    alg = nblade.Algebra.euclidean(3)
    e1, e2, e3 = alg.basis_vectors()

    print("\n向量投影将一个向量分解到另一个向量的方向上")
    print("公式: proj_u(v) = (v·u / |u|²) · u")
    print("几何代数: proj_u(v) = (v⌋u)⌋u⁻¹ = (v·u)·u / |u|²")

    # 基本向量投影
    print("\n【基本向量投影】")
    v = alg.vector([3.0, 4.0, 0.0])
    u = e1  # 投影到 x 轴方向

    proj = v.project_to(u)

    print(f"向量 v = {v}")
    print(f"投影方向 u = {u}")
    print(f"投影结果 proj_u(v) = {proj}")
    print("  注意：y 分量被移除，只剩下 x 分量")

    # 投影到斜向向量
    print("\n【投影到斜向向量】")
    v = alg.vector([3.0, 4.0, 0.0])
    u = alg.vector([1.0, 1.0, 0.0])  # 45 度方向

    proj = v.project_to(u)

    print(f"向量 v = {v}")
    print(f"投影方向 u = {u}")
    print(f"投影结果 = {proj}")

    # 验证投影长度
    proj_length = proj.norm()
    v_length = v.norm()
    angle = math.atan2(4, 3) - math.pi / 4  # v 与 45° 方向的夹角
    expected_length = v_length * math.cos(angle)
    print(f"投影长度: {proj_length:.4f}")
    print(f"理论值: {expected_length:.4f}")

    # 单位向量投影
    print("\n【投影到单位向量】")
    v = alg.vector([5.0, 0.0, 0.0])
    u = e1  # 已经是单位向量

    proj = v.project_to(u)
    print(f"v = {v}")
    print(f"投影到 e1 = {proj}")
    print(f"投影系数 = {(v | u).scalar_part():.4f}")


def section_2_subspace_projection():
    """第二节：子空间投影"""
    print("\n" + "=" * 60)
    print("第二节：子空间投影")
    print("=" * 60)

    alg = nblade.Algebra.euclidean(3)
    e1, e2, e3 = alg.basis_vectors()

    print("\n子空间投影将向量投影到更高维的子空间（如平面）")
    print("公式: proj_B(v) = (v⌋B)⌋B⁻¹，其中 B 是子空间的伪标量")

    # 投影到平面
    print("\n【投影到 xy 平面】")
    v = alg.vector([2.0, 3.0, 5.0])
    plane = e1 ^ e2  # xy 平面

    proj = v.project_to(plane)

    print(f"向量 v = {v}")
    print(f"平面 B = {plane} (xy 平面)")
    print(f"投影到平面 = {proj}")
    print("  注意：z 分量被消除")

    # 投影到 yz 平面
    print("\n【投影到 yz 平面】")
    v = alg.vector([4.0, 2.0, 3.0])
    plane = e2 ^ e3  # yz 平面

    proj = v.project_to(plane)

    print(f"向量 v = {v}")
    print(f"平面 B = {plane} (yz 平面)")
    print(f"投影到平面 = {proj}")
    print("  注意：x 分量被消除")

    # 投影到任意平面
    print("\n【投影到任意平面】")
    v = alg.vector([1.0, 2.0, 3.0])
    # 由两个向量张成的平面
    a = alg.vector([1.0, 0.0, 1.0])
    b = alg.vector([0.0, 1.0, 1.0])
    plane = a ^ b

    proj = v.project_to(plane)

    print(f"向量 v = {v}")
    print(f"平面由 a = {a} 和 b = {b} 张成")
    print(f"平面二重向量 B = a∧b = {plane}")
    print(f"投影到平面 = {proj}")

    # 验证投影在平面上
    # 投影结果应该与平面法向量正交
    normal = plane.dual()
    orthogonality = (proj | normal).scalar_part()
    print(f"正交验证 (proj·normal): {orthogonality:.2e}")


def section_3_rejection():
    """第三节：拒绝分量"""
    print("\n" + "=" * 60)
    print("第三节：拒绝分量")
    print("=" * 60)

    alg = nblade.Algebra.euclidean(3)
    e1, e2, e3 = alg.basis_vectors()

    print("\n拒绝分量是向量中与投影方向正交的部分")
    print("公式: rej_u(v) = v - proj_u(v)")
    print("几何代数: rej_u(v) = v - (v⌋u)⌋u⁻¹")

    # 基本拒绝操作
    print("\n【向量拒绝】")
    v = alg.vector([3.0, 4.0, 5.0])
    direction = e1  # x 轴方向

    proj = v.project_to(direction)
    reject = v.reject_from(direction)

    print(f"向量 v = {v}")
    print(f"方向 u = {direction}")
    print(f"投影分量 = {proj}")
    print(f"拒绝分量 = {reject}")
    print(f"验证: proj + reject = {proj + reject}")

    # 验证正交性
    orthogonality = (reject | direction).scalar_part()
    print(f"正交验证 (reject·direction): {orthogonality:.2e}")

    # 拒绝出平面
    print("\n【拒绝出平面】")
    v = alg.vector([2.0, 3.0, 4.0])
    plane = e1 ^ e2  # xy 平面

    proj = v.project_to(plane)
    reject = v.reject_from(plane)

    print(f"向量 v = {v}")
    print(f"平面 = {plane} (xy 平面)")
    print(f"投影到平面 = {proj}")
    print(f"拒绝分量（垂直于平面）= {reject}")
    print(f"验证: proj + reject = {proj + reject}")

    # 分解验证
    print("\n【完整分解】")
    v = alg.vector([1.0, 2.0, 3.0])

    # 分解到三个坐标轴
    proj_x = v.project_to(e1)
    proj_y = v.project_to(e2)
    proj_z = v.project_to(e3)

    print(f"向量 v = {v}")
    print(f"x 分量: {proj_x}")
    print(f"y 分量: {proj_y}")
    print(f"z 分量: {proj_z}")
    print(f"重构: {proj_x + proj_y + proj_z}")

    # 拒绝的拒绝
    print("\n【连续拒绝】")
    v = alg.vector([3.0, 4.0, 5.0])

    # 先拒绝 x 方向
    reject_x = v.reject_from(e1)
    print(f"原始向量: {v}")
    print(f"拒绝 x 方向后: {reject_x}")

    # 再拒绝 y 方向（从剩余部分）
    reject_xy = reject_x.reject_from(e2)
    print(f"再拒绝 y 方向后: {reject_xy}")
    print(f"最终结果仅包含 z 分量")


def section_4_practical_application():
    """第四节：实际应用 - 阴影计算"""
    print("\n" + "=" * 60)
    print("第四节：实际应用 - 阴影计算")
    print("=" * 60)

    alg = nblade.Algebra.euclidean(3)
    e1, e2, e3 = alg.basis_vectors()

    print("\n投影在计算机图形学中用于计算阴影和投影效果")
    print("平行光阴影 = 将物体顶点投影到地面平面")

    # 场景设置
    print("\n【场景设置】")
    ground_plane = e1 ^ e2  # 地面是 xy 平面
    light_direction = alg.vector([1.0, -2.0, -1.0])  # 斜向光源

    print(f"地面平面: {ground_plane} (z = 0)")
    print(f"光源方向: {light_direction}")

    # 定义一个简单物体（三角形）
    print("\n【物体顶点】")
    p1 = alg.vector([1.0, 1.0, 3.0])
    p2 = alg.vector([3.0, 1.0, 3.0])
    p3 = alg.vector([2.0, 3.0, 3.0])

    print(f"三角形顶点:")
    print(f"  p1 = {p1}")
    print(f"  p2 = {p2}")
    print(f"  p3 = {p3}")

    # 计算阴影（正交投影到地面）
    print("\n【正交投影阴影】")
    print("假设光源在正上方无限远处")

    shadow_p1 = p1.project_to(ground_plane)
    shadow_p2 = p2.project_to(ground_plane)
    shadow_p3 = p3.project_to(ground_plane)

    print(f"阴影顶点:")
    print(f"  p1' = {shadow_p1}")
    print(f"  p2' = {shadow_p2}")
    print(f"  p3' = {shadow_p3}")

    # 斜投影（考虑光源方向）
    print("\n【斜投影阴影（考虑光源方向）】")
    print("公式: P_shadow = P - ((P·n) / (d·n)) · d")
    print("其中 n 是平面法向量，d 是光源方向")

    # 平面法向量
    n = ground_plane.dual()
    n_unit = n.scale(1.0 / n.norm())

    # 计算斜投影
    def oblique_projection(point, direction, plane_normal):
        """计算点到平面的斜投影"""
        # t = -(point · n) / (direction · n)
        # projection = point + t * direction
        point_dot_n = (point | plane_normal).scalar_part()
        dir_dot_n = (direction | plane_normal).scalar_part()
        t = -point_dot_n / dir_dot_n
        return point + direction.scale(t)

    oblique_p1 = oblique_projection(p1, light_direction, n_unit)
    oblique_p2 = oblique_projection(p2, light_direction, n_unit)
    oblique_p3 = oblique_projection(p3, light_direction, n_unit)

    print(f"斜投影阴影顶点:")
    print(f"  p1'' = {oblique_p1}")
    print(f"  p2'' = {oblique_p2}")
    print(f"  p3'' = {oblique_p3}")

    # 投影到倾斜平面
    print("\n【投影到倾斜墙面】")
    # 倾斜墙面（yz 平面旋转 30 度）
    angle = math.pi / 6  # 30 度
    wall_normal = e1.scale(math.cos(angle)) + e3.scale(math.sin(angle))
    wall_plane = (e2 ^ wall_normal).dual()  # 垂直于法向量的平面

    print(f"墙面法向量: {wall_normal}")
    print(f"墙面平面: {wall_plane}")

    proj_p1 = p1.project_to(wall_plane)
    proj_p2 = p2.project_to(wall_plane)
    proj_p3 = p3.project_to(wall_plane)

    print(f"投影到墙面的顶点:")
    print(f"  p1' = {proj_p1}")
    print(f"  p2' = {proj_p2}")
    print(f"  p3' = {proj_p3}")

    # 计算阴影面积
    print("\n【阴影面积计算】")
    # 使用叉积计算三角形面积
    edge1 = shadow_p2 - shadow_p1
    edge2 = shadow_p3 - shadow_p1
    area = (edge1 ^ edge2).norm() / 2
    print(f"阴影三角形面积: {area:.4f}")

    # 原始三角形面积
    orig_edge1 = p2 - p1
    orig_edge2 = p3 - p1
    orig_area = (orig_edge1 ^ orig_edge2).norm() / 2
    print(f"原始三角形面积: {orig_area:.4f}")


def main():
    """主函数"""
    print("=" * 60)
    print("nblade 投影操作示例")
    print("=" * 60)

    section_1_vector_projection()
    section_2_subspace_projection()
    section_3_rejection()
    section_4_practical_application()

    print("\n" + "=" * 60)
    print("示例完成！")
    print("=" * 60)
    print("\n关键点:")
    print("  - 向量投影: proj_u(v) = (v⌋u)⌋u⁻¹")
    print("  - 子空间投影: proj_B(v) = (v⌋B)⌋B⁻¹")
    print("  - 拒绝分量: rej_u(v) = v - proj_u(v)")
    print("  - 投影 + 拒绝 = 原始向量")
    print("  - 应用: 阴影计算、正交分解")


if __name__ == "__main__":
    main()
