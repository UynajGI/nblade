"""
高级几何代数运算示例

演示几何代数的高级运算：
- 投影和拒绝
- 反射和旋转
- 对偶运算
- 逆运算
"""

from nblade import Algebra
import math


def demo_advanced_operations():
    """演示高级运算"""
    print("=== 高级几何代数运算示例 ===\n")

    # 创建 3D 欧几里得代数
    alg = Algebra.euclidean(3)
    e1, e2, e3 = alg.basis_vectors()

    print("基向量:")
    print(f"  e1 = {e1}")
    print(f"  e2 = {e2}")
    print(f"  e3 = {e3}\n")

    # 创建向量和平面
    v = alg.vector([1, 2, 3])
    plane = e1 ^ e2  # xy 平面
    print(f"向量 v = {v}")
    print(f"平面 e1∧e2 = {plane}\n")

    # 正交投影
    print("正交投影:")
    proj_v_to_plane = v.project_to(plane)
    print(f"  v 投影到平面: {proj_v_to_plane}")

    # 正交拒绝
    rej_v_from_plane = v.reject_from(plane)
    print(f"  v 拒绝自平面: {rej_v_from_plane}")

    # 验证 v = proj + rej
    sum_proj_rej = proj_v_to_plane + rej_v_from_plane
    print(f"  投影 + 拒绝: {sum_proj_rej}")
    print(f"  误差: {((v - sum_proj_rej).norm_squared()) ** 0.5:.10f}\n")

    # 反射
    print("反射:")
    reflected_v = v.reflect_in(plane)
    print(f"  v 反射在平面中: {reflected_v}")

    # 反射两次应该回到原向量
    reflected_twice = reflected_v.reflect_in(plane)
    print(f"  v 反射两次: {reflected_twice}")
    error = ((v - reflected_twice).norm_squared()) ** 0.5
    print(f"  误差: {error:.10f}\n")

    # 旋转
    print("旋转:")
    # 创建 90 度旋转的转子
    rotor = alg.rotor(plane, math.pi / 2)  # 90 度
    print(f"  90° 旋转转子: {rotor}")

    # 旋转向量
    rotated_v = v.rotate_by(rotor)
    print(f"  v 旋转 90°: {rotated_v}")

    # 旋转基向量
    rotated_e1 = e1.rotate_by(rotor)
    print(f"  e1 旋转 90°: {rotated_e1}")
    print(f"  误差 (应接近 e2): {((rotated_e1 - e2).norm_squared()) ** 0.5:.10f}\n")

    # 对偶运算
    print("对偶运算:")
    # 在 3D 中，向量的对偶是二重向量
    dual_v = v.dual()
    print(f"  v 的对偶: {dual_v}")

    # 二重向量的对偶是向量
    dual_plane = plane.dual()
    print(f"  e1∧e2 的对偶: {dual_plane}")

    # 双重对偶性质
    dual_dual_v = v.dual().dual()
    print(f"  v 双重对偶: {dual_dual_v}")
    print(f"  误差 (3D 中应为 -v): {((v + dual_dual_v).norm_squared()) ** 0.5:.10f}\n")

    # 逆运算
    print("逆运算:")
    # 创建可逆多向量
    mv = alg.scalar(1.0) + e1 + (e2 ^ e3)
    print(f"  多向量: {mv}")

    try:
        inv_mv = mv.inverse()
        print(f"  逆: {inv_mv}")

        # 验证 A*A⁻¹ = 1
        product = mv * inv_mv
        print(f"  AA⁻¹ = {product}")
        error = ((product - alg.scalar(1.0)).norm_squared()) ** 0.5
        print(f"  误差: {error:.10f}\n")
    except Exception as e:
        print(f"  无法计算逆: {e}\n")

    # 复杂多向量运算
    print("复杂多向量运算:")
    # 创建二重向量
    bivector = e1 ^ e2 + 2 * (e2 ^ e3)
    print(f"  二重向量: {bivector}")

    # 二重向量的平方
    bivector_sq = bivector * bivector
    print(f"  平方: {bivector_sq}")

    # 二重向量的对合
    print(f"  阶次对合: {~bivector}")
    print(f"  反转: {bivector.reversion()}")
    print(f"  Clifford 共轭: {bivector.clifford_conjugate()}\n")

    print("高级运算示例完成！\n")


def demo_higher_dimensions():
    """演示高维空间运算"""
    print("=== 高维几何代数运算示例 ===\n")

    # 创建 4D 时空代数 G(1,3,0)
    alg_4d = Algebra(4, p=1, q=3, r=0)  # 时空代数
    print(f"时空代数: {alg_4d}")
    print(f"签名: {alg_4d.signature}\n")

    # 获取基向量
    e0, e1, e2, e3 = alg_4d.basis_vectors()
    print("时空基向量:")
    print(f"  e0² = {(e0 * e0).scalar_part():.1f} (时间)")  # 应该是 +1
    print(f"  e1² = {(e1 * e1).scalar_part():.1f} (空间)")  # 应该是 -1
    print(f"  e2² = {(e2 * e2).scalar_part():.1f} (空间)")  # 应该是 -1
    print(f"  e3² = {(e3 * e3).scalar_part():.1f} (空间)")  # 应该是 -1\n")

    # 4D 向量
    four_vector = alg_4d.vector([1, 0.5, 0.3, 0.2])  # [t, x, y, z]
    print(f"四维向量 (时空向量): {four_vector}")

    # 计算范数平方 (闵可夫斯基范数)
    norm_sq = four_vector.norm_squared()
    print(f"范数平方: {norm_sq:.6f}\n")

    # 创建 5D 共形几何代数 G(4,1,0)
    cga = Algebra(5, p=4, q=1, r=0)  # 共形几何代数
    print(f"共形几何代数: {cga}")
    print(f"维度: {cga.dimension}, 签名: {cga.signature}\n")

    # 共形空间基向量
    e1, e2, e3, e4, e5 = cga.basis_vectors()
    print("共形基向量:")
    print(f"  e1² = {(e1 * e1).scalar_part():.1f}")
    print(f"  e2² = {(e2 * e2).scalar_part():.1f}")
    print(f"  e3² = {(e3 * e3).scalar_part():.1f}")
    print(f"  e4² = {(e4 * e4).scalar_part():.1f}")
    print(f"  e5² = {(e5 * e5).scalar_part():.1f}\n")

    print("高维运算示例完成！\n")


if __name__ == "__main__":
    demo_advanced_operations()
    demo_higher_dimensions()
