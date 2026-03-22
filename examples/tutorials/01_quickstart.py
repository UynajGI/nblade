"""
nblade 快速入门示例
===================

5 分钟快速了解 nblade 的核心功能：
- 创建几何代数
- 向量和基向量
- 三种基本乘积
- 旋转操作

作者: nblade 团队
"""

import nblade
import math


def main():
    """快速入门主函数"""
    print("=" * 60)
    print("nblade 快速入门 - 5 分钟掌握核心功能")
    print("=" * 60)

    # ========================================
    # 1. 创建几何代数
    # ========================================
    print("\n【1. 创建几何代数】")
    print("-" * 40)

    # 创建 3D 欧几里得几何代数 G(3,0,0)
    # G(p,q,r) 表示有 p 个正平方基向量，q 个负平方，r 个零平方
    alg = nblade.Algebra.euclidean(3)

    print(f"代数类型: {alg}")
    print(f"维度: {alg.dimension}")
    print(f"签名 (p, q, r): {alg.signature}")

    # 其他预定义代数
    # alg_2d = nblade.Algebra.euclidean(2)   # 2D 欧几里得
    # alg_st = nblade.Algebra.spacetime(4)   # 时空代数 G(1,3,0)
    # alg_cga = nblade.Algebra.cga()         # 共形几何代数 G(4,1,0)

    # ========================================
    # 2. 基向量和向量创建
    # ========================================
    print("\n【2. 基向量和向量创建】")
    print("-" * 40)

    # 获取所有基向量
    e1, e2, e3 = alg.basis_vectors()
    print(f"基向量 e1 = {e1}")
    print(f"基向量 e2 = {e2}")
    print(f"基向量 e3 = {e3}")

    # 从列表创建向量
    v = alg.vector([1.0, 2.0, 3.0])
    print(f"\n向量 v = alg.vector([1, 2, 3])")
    print(f"      = {v}")

    # 创建标量
    s = alg.scalar(5.0)
    print(f"\n标量 s = alg.scalar(5.0) = {s}")

    # ========================================
    # 3. 三种基本乘积
    # ========================================
    print("\n【3. 三种基本乘积】")
    print("-" * 40)

    # 创建两个向量
    a = alg.vector([1.0, 0.0, 0.0])  # 沿 e1 方向
    b = alg.vector([1.0, 1.0, 0.0])  # 在 xy 平面内

    print(f"向量 a = {a}")
    print(f"向量 b = {b}")

    # 几何积: ab = a·b + a∧b (内积 + 外积)
    geometric = a * b
    print(f"\n几何积 a * b = {geometric}")
    print("  (包含标量部分和二重向量部分)")

    # 外积 (楔积): 反对称部分，表示张成的子空间
    outer = a ^ b
    print(f"\n外积 a ^ b = {outer}")
    print("  (二重向量，表示 a 和 b 张成的平面)")

    # 内积 (左收缩): 对称部分，表示投影
    inner = a | b
    print(f"\n内积 a | b = {inner}")
    print("  (标量，表示 a 在 b 上的投影)")

    # 验证基本关系式: ab = a·b + a∧b
    print("\n验证 ab = a·b + a∧b:")
    reconstructed = inner + outer
    print(f"  内积 + 外积 = {reconstructed}")
    print(f"  几何积      = {geometric}")

    # ========================================
    # 4. 旋转操作
    # ========================================
    print("\n【4. 旋转操作】")
    print("-" * 40)

    # 创建旋转平面 (xy 平面)
    plane = e1 ^ e2
    print(f"旋转平面: {plane}")

    # 创建转子: 在 xy 平面内旋转 45 度
    angle = math.pi / 4  # 45 度
    rotor = alg.rotor(plane, angle)
    print(f"转子 (旋转 45°): {rotor}")

    # 使用转子旋转向量
    # 旋转公式: v' = R * v * R⁻¹
    rotated = e1.rotate_by(rotor)
    print(f"\ne1 旋转 45° 后: {rotated}")

    # 理论验证: e1 旋转 45° 应该是 cos(45°)e1 + sin(45°)e2
    cos_45 = math.cos(angle)
    sin_45 = math.sin(angle)
    expected = alg.vector([cos_45, sin_45, 0.0])
    print(f"理论预期: {expected}")

    # 计算误差
    diff = rotated - expected
    error = diff.norm()
    print(f"误差: {error:.2e}")

    # ========================================
    # 5. 其他有用操作
    # ========================================
    print("\n【5. 其他有用操作】")
    print("-" * 40)

    # 范数
    v = alg.vector([3.0, 4.0, 0.0])
    print(f"向量 v = {v}")
    print(f"范数 |v| = {v.norm():.6f}  (应为 5)")
    print(f"范数平方 |v|² = {v.norm_squared():.6f}")

    # 对偶: 在 3D 中，向量 a 的对偶是二重向量
    dual_v = v.dual()
    print(f"\n对偶 dual(v) = {dual_v}")

    # 逆向量
    v_inv = v.inverse()
    print(f"逆 v⁻¹ = {v_inv}")
    print(f"验证 v * v⁻¹ = {(v * v_inv).scalar_part():.6f}  (应为 1)")

    # 阶次操作
    # 创建混合多向量 (标量 + 向量 + 二重向量)
    mixed = alg.scalar(1.0) + e1 + (e1 ^ e2)
    print(f"\n混合多向量: {mixed}")
    print(f"  标量部分: {mixed.grade(0)}")
    print(f"  向量部分: {mixed.grade(1)}")
    print(f"  二重向量部分: {mixed.grade(2)}")

    # ========================================
    # 总结
    # ========================================
    print("\n" + "=" * 60)
    print("快速入门完成！")
    print("=" * 60)
    print("\n核心要点:")
    print("  1. 使用 nblade.Algebra.euclidean(n) 创建 n 维几何代数")
    print("  2. 使用 alg.vector([...]) 创建向量")
    print("  3. 三种乘积: * (几何积), ^ (外积), | (内积)")
    print("  4. 使用转子实现旋转: rotor = alg.rotor(plane, angle)")
    print("\n下一步:")
    print("  - 运行 02_basic_operations.py 了解更多运算细节")
    print("  - 运行 05_rotations.py 深入学习旋转")
    print("  - 查看 docs/ 目录获取完整文档")


if __name__ == "__main__":
    main()
