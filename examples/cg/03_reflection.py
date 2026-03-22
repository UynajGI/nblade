"""
nblade 反射操作示例
==================

演示几何代数中的反射操作：
- 向量反射（在超平面/向量法向）
- 平面反射
- 多次反射的组合 = 旋转
- 实际应用：镜面反射模拟

反射公式：v' = -n·v·n （反射在向量 n 上，即反射在垂直于 n 的超平面）
         v' = B·v·B （反射在二重向量 B 上）

作者: nblade 团队
"""

import nblade
import math


def section_1_vector_reflection():
    """第一节：向量反射（超平面反射）"""
    print("\n" + "=" * 60)
    print("第一节：向量反射（超平面反射）")
    print("=" * 60)

    alg = nblade.Algebra.euclidean(3)
    e1, e2, e3 = alg.basis_vectors()

    print("\n反射是最基本的几何变换，所有等距变换都可以由反射组合而成")
    print("在几何代数中，反射公式为: v' = -n·v·n")
    print("其中 n 是反射超平面的单位法向量")

    # 基本反射
    print("\n【基本反射】")
    v = alg.vector([1.0, 2.0, 3.0])
    n = e1  # 反射面法向量（yz 平面）

    reflected = v.reflect_in(n)

    print(f"原始向量 v = {v}")
    print(f"反射面法向量 n = {n}")
    print(f"反射结果 v' = {reflected}")
    print("  注意：x 分量被取反")

    # 反射在 y 轴方向
    print("\n【反射在 xz 平面】")
    v = alg.vector([2.0, 3.0, 4.0])
    n = e2  # y 轴方向，即 xz 平面

    reflected = v.reflect_in(n)

    print(f"原始向量 v = {v}")
    print(f"反射面法向量 n = {n} (xz 平面)")
    print(f"反射结果 v' = {reflected}")
    print("  注意：y 分量被取反")

    # 反射在斜向平面
    print("\n【反射在斜向平面】")
    v = alg.vector([1.0, 0.0, 0.0])
    # 45 度平面
    n = alg.vector([1.0, 1.0, 0.0])
    n_unit = n.scale(1.0 / n.norm())

    reflected = v.reflect_in(n_unit)

    print(f"原始向量 v = {v}")
    print(f"反射面法向量 n = {n_unit} (45 度平面)")
    print(f"反射结果 v' = {reflected}")

    # 验证反射性质
    print("\n【反射性质验证】")
    v = alg.vector([3.0, 4.0, 5.0])
    n = e1

    reflected = v.reflect_in(n)
    reflected_again = reflected.reflect_in(n)

    print(f"原始向量: {v}")
    print(f"第一次反射: {reflected}")
    print(f"第二次反射: {reflected_again}")
    print(f"两次反射 = 恒等变换: {reflected_again == v}")

    # 验证保持长度
    print(f"\n长度保持验证:")
    print(f"原始长度: {v.norm():.4f}")
    print(f"反射后长度: {reflected.norm():.4f}")


def section_2_plane_reflection():
    """第二节：平面反射"""
    print("\n" + "=" * 60)
    print("第二节：平面反射")
    print("=" * 60)

    alg = nblade.Algebra.euclidean(3)
    e1, e2, e3 = alg.basis_vectors()

    print("\n在二重向量（平面）上的反射")
    print("公式: v' = B·v·B，其中 B 是单位二重向量")

    # 在 xy 平面上反射
    print("\n【在 xy 平面上反射】")
    v = alg.vector([1.0, 2.0, 3.0])
    plane = e1 ^ e2  # xy 平面

    reflected = v.reflect_in(plane)

    print(f"原始向量 v = {v}")
    print(f"反射平面 B = {plane} (xy 平面)")
    print(f"反射结果 v' = {reflected}")
    print("  注意：垂直于平面的分量（z）被取反")

    # 在 yz 平面上反射
    print("\n【在 yz 平面上反射】")
    v = alg.vector([3.0, 2.0, 1.0])
    plane = e2 ^ e3  # yz 平面

    reflected = v.reflect_in(plane)

    print(f"原始向量 v = {v}")
    print(f"反射平面 B = {plane} (yz 平面)")
    print(f"反射结果 v' = {reflected}")

    # 在任意平面上反射
    print("\n【在任意平面上反射】")
    v = alg.vector([1.0, 0.0, 0.0])
    # 创建一个倾斜平面
    a = alg.vector([1.0, 1.0, 0.0])
    b = alg.vector([0.0, 1.0, 1.0])
    plane = a ^ b

    # 归一化平面
    plane_norm = plane.norm()
    plane_unit = plane.scale(1.0 / plane_norm)

    reflected = v.reflect_in(plane_unit)

    print(f"原始向量 v = {v}")
    print(f"反射平面 B = {plane} (由 {a} 和 {b} 张成)")
    print(f"反射结果 v' = {reflected}")

    # 比较两种反射方式
    print("\n【向量反射 vs 平面反射】")
    v = alg.vector([2.0, 3.0, 4.0])

    # 向量反射：在 e3 上反射（即 xy 平面）
    reflected_vector = v.reflect_in(e3)
    # 平面反射：在 e1^e2 平面上反射（也是 xy 平面）
    reflected_plane = v.reflect_in(e1 ^ e2)

    print(f"原始向量: {v}")
    print(f"在 e3 上反射（向量反射）: {reflected_vector}")
    print(f"在 e1∧e2 上反射（平面反射）: {reflected_plane}")
    print("  两种方法等价！")


def section_3_multiple_reflections():
    """第三节：多次反射的组合"""
    print("\n" + "=" * 60)
    print("第三节：多次反射的组合 = 旋转")
    print("=" * 60)

    alg = nblade.Algebra.euclidean(3)
    e1, e2, e3 = alg.basis_vectors()

    print("\n两次反射的组合等价于一次旋转")
    print("旋转角度 = 两反射面夹角的两倍")

    # 两次反射 = 旋转
    print("\n【两次反射 = 旋转】")
    v = e1

    # 第一个反射面：xz 平面（法向量 e2）
    n1 = e2
    # 第二个反射面：与 xz 平面成 30 度的平面
    angle_between = math.pi / 6  # 30 度
    n2 = e2.scale(math.cos(angle_between)) + e3.scale(math.sin(angle_between))
    n2 = n2.scale(1.0 / n2.norm())

    print(f"原始向量 v = {v}")
    print(f"第一次反射面法向量 n1 = {n1}")
    print(f"第二次反射面法向量 n2 = {n2}")
    print(f"两平面夹角: {math.degrees(angle_between):.1f}°")

    # 依次反射
    reflected1 = v.reflect_in(n1)
    reflected2 = reflected1.reflect_in(n2)

    print(f"\n第一次反射后: {reflected1}")
    print(f"第二次反射后: {reflected2}")

    # 计算等效旋转
    # 两次反射 = 绕两平面交线旋转 2*angle_between
    rotation_angle = 2 * angle_between
    print(f"\n等效旋转角度: {math.degrees(rotation_angle):.1f}°")

    # 验证与转子旋转等价
    # 交线方向 = n1 × n2
    intersection = (n1 ^ n2).dual()
    print(f"旋转轴方向: {intersection}")

    # 使用转子验证
    rotation_plane = e2 ^ e3  # 旋转平面
    rotor = alg.rotor(rotation_plane, rotation_angle)
    rotated_by_rotor = v.rotate_by(rotor)

    print(f"转子旋转结果: {rotated_by_rotor}")
    print(f"两次反射结果: {reflected2}")

    # 三次反射
    print("\n【三次反射 = 反射 + 旋转】")
    v = alg.vector([1.0, 0.0, 0.0])
    n1 = e1
    n2 = alg.vector([1.0, 1.0, 0.0]).scale(1.0 / math.sqrt(2))
    n3 = e2

    reflected = v.reflect_in(n1).reflect_in(n2).reflect_in(n3)

    print(f"原始向量: {v}")
    print(f"三次反射后: {reflected}")
    print("  奇数次反射 = 反射（改变手性）")

    # 四次反射 = 恒等变换（或双旋转）
    print("\n【四次反射】")
    v = alg.vector([1.0, 2.0, 3.0])
    reflected = v.reflect_in(e1).reflect_in(e2).reflect_in(e1).reflect_in(e2)

    print(f"原始向量: {v}")
    print(f"四次反射后: {reflected}")
    print("  偶数次反射保持手性")


def section_4_practical_example():
    """第四节：实际应用 - 镜面反射模拟"""
    print("\n" + "=" * 60)
    print("第四节：实际应用 - 镜面反射模拟")
    print("=" * 60)

    alg = nblade.Algebra.euclidean(3)
    e1, e2, e3 = alg.basis_vectors()

    print("\n模拟光线在镜面表面的反射")
    print("入射角 = 反射角（相对于法向量）")

    # 场景设置
    print("\n【场景设置】")
    # 镜面平面（xy 平面）
    mirror_plane = e1 ^ e2
    mirror_normal = e3

    print(f"镜面平面: {mirror_plane} (xy 平面)")
    print(f"镜面法向量: {mirror_normal}")

    # 入射光线
    print("\n【入射光线】")
    # 从斜上方射入的光线
    incident = alg.vector([1.0, 0.0, -1.0])
    incident = incident.scale(1.0 / incident.norm())  # 单位化

    print(f"入射光线方向: {incident}")

    # 计算反射光线
    reflected = incident.reflect_in(mirror_normal)

    print(f"反射光线方向: {reflected}")

    # 验证入射角 = 反射角
    incident_angle = math.acos(abs((incident | mirror_normal).scalar_part()))
    reflected_angle = math.acos(abs((reflected | mirror_normal).scalar_part()))

    print(f"\n入射角: {math.degrees(incident_angle):.2f}°")
    print(f"反射角: {math.degrees(reflected_angle):.2f}°")

    # 多面镜反射
    print("\n【多面镜反射】")
    print("模拟光线在两个垂直镜面之间的反射")

    # 两个垂直的镜面
    mirror1_normal = e1  # yz 平面
    mirror2_normal = e2  # xz 平面

    # 入射光线
    ray = alg.vector([1.0, 1.0, -1.0])
    ray = ray.scale(1.0 / ray.norm())

    print(f"初始光线: {ray}")

    # 在第一个镜面反射
    ray_after_m1 = ray.reflect_in(mirror1_normal)
    print(f"在镜面1反射后: {ray_after_m1}")

    # 在第二个镜面反射
    ray_after_m2 = ray_after_m1.reflect_in(mirror2_normal)
    print(f"在镜面2反射后: {ray_after_m2}")

    # 计算等效旋转
    print(f"\n两次垂直反射 = 绕交线旋转 180°")
    print(f"交线方向: e3 = {e3}")

    # 角落反射器（ retroreflector ）
    print("\n【角反射器（回射器）】")
    print("三个互相垂直的镜面，光线总是原路返回")

    # 三个垂直镜面
    ray = alg.vector([1.0, 2.0, 3.0])
    ray = ray.scale(1.0 / ray.norm())

    print(f"入射光线: {ray}")

    # 在三个垂直镜面依次反射
    reflected = ray.reflect_in(e1).reflect_in(e2).reflect_in(e3)
    print(f"三次反射后: {reflected}")
    print(f"方向是否相反: {abs((reflected + ray).norm()) < 1e-10}")

    # 物体反射
    print("\n【物体反射】")
    print("对一个三角形进行镜像变换")

    # 定义三角形
    p1 = alg.vector([1.0, 1.0, 2.0])
    p2 = alg.vector([3.0, 1.0, 2.0])
    p3 = alg.vector([2.0, 3.0, 2.0])

    print(f"原始三角形:")
    print(f"  p1 = {p1}")
    print(f"  p2 = {p2}")
    print(f"  p3 = {p3}")

    # 在 yz 平面反射
    mirror = e1

    p1_mir = p1.reflect_in(mirror)
    p2_mir = p2.reflect_in(mirror)
    p3_mir = p3.reflect_in(mirror)

    print(f"\n镜像后（yz 平面）:")
    print(f"  p1' = {p1_mir}")
    print(f"  p2' = {p2_mir}")
    print(f"  p3' = {p3_mir}")

    # 计算面积（反射保持面积）
    edge1 = p2 - p1
    edge2 = p3 - p1
    area_orig = (edge1 ^ edge2).norm() / 2

    edge1_mir = p2_mir - p1_mir
    edge2_mir = p3_mir - p1_mir
    area_mir = (edge1_mir ^ edge2_mir).norm() / 2

    print(f"\n原始面积: {area_orig:.4f}")
    print(f"镜像面积: {area_mir:.4f}")

    # 法向量方向（反射改变手性）
    normal_orig = (edge1 ^ edge2).dual()
    normal_mir = (edge1_mir ^ edge2_mir).dual()

    print(f"原始法向量: {normal_orig}")
    print(f"镜像法向量: {normal_mir}")
    print("  注意：法向量方向相反（手性改变）")


def main():
    """主函数"""
    print("=" * 60)
    print("nblade 反射操作示例")
    print("=" * 60)

    section_1_vector_reflection()
    section_2_plane_reflection()
    section_3_multiple_reflections()
    section_4_practical_example()

    print("\n" + "=" * 60)
    print("示例完成！")
    print("=" * 60)
    print("\n关键点:")
    print("  - 向量反射: v' = -n·v·n (在垂直于 n 的超平面反射)")
    print("  - 平面反射: v' = B·v·B (在二重向量 B 表示的平面反射)")
    print("  - 两次反射 = 旋转 (角度 = 2×反射面夹角)")
    print("  - 奇数次反射改变手性，偶数次反射保持手性")
    print("  - 应用: 镜面反射、光线追踪、角反射器")


if __name__ == "__main__":
    main()
