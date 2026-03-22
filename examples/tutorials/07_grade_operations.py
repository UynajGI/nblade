"""
nblade 阶次操作教程
===================

详细介绍多向量阶次的概念和操作：
- 阶次的概念（0阶到n阶）
- 阶次提取
- 偶部和奇部分解
- 阶次对合（Grade Involution）
- 反序（Reversion）

作者: nblade 团队
"""

import nblade
import math


def section_1_what_are_grades():
    """第一节：什么是阶次"""
    print("\n" + "=" * 60)
    print("第一节：什么是阶次 (Grade)")
    print("=" * 60)

    alg = nblade.Algebra.euclidean(3)
    e1, e2, e3 = alg.basis_vectors()

    print("\n在多向量中，阶次表示基向量的数量：")
    print("  - 0阶：标量（纯数值）")
    print("  - 1阶：向量（一个基向量）")
    print("  - 2阶：二重向量（两个基向量楔积）")
    print("  - 3阶：三重向量（三个基向量楔积）")
    print("  - n阶：n重向量（n个基向量楔积）")

    # 0阶：标量
    print("\n【0阶：标量 (Grade 0)】")
    scalar = alg.scalar(5.0)
    print(f"标量: {scalar}")
    print("  不包含任何基向量，是纯数值")

    # 1阶：向量
    print("\n【1阶：向量 (Grade 1)】")
    vector = e1 + 2 * e2
    print(f"向量: {vector}")
    print("  包含一个基向量，表示方向和大小")

    # 2阶：二重向量
    print("\n【2阶：二重向量 (Grade 2)】")
    bivector = e1 ^ e2
    print(f"二重向量: {bivector}")
    print("  两个基向量的楔积，表示平面元素")

    # 3阶：三重向量
    print("\n【3阶：三重向量 (Grade 3)】")
    trivector = e1 ^ e2 ^ e3
    print(f"三重向量: {trivector}")
    print("  三个基向量的楔积，表示体积元素")

    # 伪标量
    print("\n【伪标量 (Pseudoscalar)】")
    I = alg.config.volume_element()
    print(f"3D伪标量 I = e1∧e2∧e3 = {I}")
    print("  在n维空间中，n阶元素称为伪标量")

    # 混合多向量
    print("\n【混合多向量】")
    mixed = alg.scalar(2.0) + e1 + (e1 ^ e2) + (e1 ^ e2 ^ e3)
    print(f"混合多向量: {mixed}")
    print("  同时包含多个阶次的元素")


def section_2_grade_extraction():
    """第二节：阶次提取"""
    print("\n" + "=" * 60)
    print("第二节：阶次提取")
    print("=" * 60)

    alg = nblade.Algebra.euclidean(3)
    e1, e2, e3 = alg.basis_vectors()

    print("\n使用 grade(r) 方法可以提取多向量的 r 阶部分")

    # 创建包含所有阶次的混合多向量
    print("\n【创建混合多向量】")
    mixed = (
        alg.scalar(3.0)  # 0阶
        + e1.scale(2.0)
        + e2.scale(4.0)  # 1阶
        + (e1 ^ e2).scale(5.0)  # 2阶
        + (e1 ^ e2 ^ e3).scale(6.0)
    )  # 3阶

    print(f"混合多向量: {mixed}")
    print("  包含 0阶、1阶、2阶、3阶 部分")

    # 提取各阶部分
    print("\n【提取各阶部分】")

    grade_0 = mixed.grade(0)
    print(f"grade(0) - 标量部分: {grade_0}")

    grade_1 = mixed.grade(1)
    print(f"grade(1) - 向量部分: {grade_1}")

    grade_2 = mixed.grade(2)
    print(f"grade(2) - 二重向量部分: {grade_2}")

    grade_3 = mixed.grade(3)
    print(f"grade(3) - 三重向量部分: {grade_3}")

    # 验证：将各阶部分相加应得到原多向量
    print("\n【验证】")
    reconstructed = grade_0 + grade_1 + grade_2 + grade_3
    print(f"各阶之和: {reconstructed}")
    print(f"原多向量: {mixed}")
    print(f"误差: {(reconstructed - mixed).norm():.2e}")

    # 提取不存在的阶次
    print("\n【提取不存在的阶次】")
    grade_5 = mixed.grade(5)
    print(f"grade(5) (3D中不存在): {grade_5}")
    print("  不存在的阶次返回零多向量")

    # 在更高维度中
    print("\n【在5D中提取阶次】")
    alg_5d = nblade.Algebra.euclidean(5)
    basis_5d = alg_5d.basis_vectors()
    e1_5d, e2_5d = basis_5d[0], basis_5d[1]

    # 创建一个4阶元素
    quadvector = e1_5d ^ e2_5d ^ basis_5d[2] ^ basis_5d[3]
    print(f"4阶元素 (四重向量): {quadvector}")
    print(f"提取4阶: {quadvector.grade(4)}")


def section_3_even_odd_parts():
    """第三节：偶部和奇部"""
    print("\n" + "=" * 60)
    print("第三节：偶部和奇部分解")
    print("=" * 60)

    alg = nblade.Algebra.euclidean(3)
    e1, e2, e3 = alg.basis_vectors()

    print("\n多向量可以分解为偶部和奇部：")
    print("  - 偶部 (even_part): 包含偶数阶（0, 2, 4, ...）")
    print("  - 奇部 (odd_part): 包含奇数阶（1, 3, 5, ...）")

    # 创建混合多向量
    print("\n【创建混合多向量】")
    mixed = (
        alg.scalar(1.0)  # 0阶（偶）
        + e1
        + e2
        + e3  # 1阶（奇）
        + (e1 ^ e2)
        + (e2 ^ e3)  # 2阶（偶）
        + (e1 ^ e2 ^ e3)
    )  # 3阶（奇）

    print(f"混合多向量: {mixed}")

    # 提取偶部和奇部
    print("\n【提取偶部和奇部】")
    even = mixed.even_part()
    odd = mixed.odd_part()

    print(f"偶部 (even_part): {even}")
    print("  包含 0阶（标量）和 2阶（二重向量）")

    print(f"\n奇部 (odd_part): {odd}")
    print("  包含 1阶（向量）和 3阶（三重向量）")

    # 验证
    print("\n【验证】")
    print(f"偶部 + 奇部 = {even + odd}")
    print(f"原多向量 = {mixed}")
    error = (even + odd - mixed).norm()
    print(f"误差: {error:.2e}")

    # 偶部和奇部是正交的
    print("\n【正交性】")
    geometric = even * odd
    print(f"偶部 * 奇部 = {geometric}")
    print("  偶元素和奇元素的几何积产生奇元素")

    # 为什么转子必须是偶部
    print("\n【应用：转子 (Rotor) 必须是偶元素】")
    plane = e1 ^ e2
    rotor = alg.rotor(plane, math.pi / 4)
    print(f"转子: {rotor}")

    # 验证转子是纯偶元素
    rotor_even = rotor.even_part()
    rotor_odd = rotor.odd_part()
    print(f"转子偶部: {rotor_even}")
    print(f"转子奇部范数: {rotor_odd.norm():.2e} (应接近 0)")
    print("  转子只包含偶阶元素（标量 + 二重向量）")

    # 偶部构成子代数
    print("\n【子代数性质】")
    print("偶部元素在几何积下封闭，构成子代数")
    even_element = alg.scalar(1.0) + (e1 ^ e2).scale(0.5)
    product = even_element * even_element
    product_odd = product.odd_part()
    print(f"偶元素 * 偶元素的奇部范数: {product_odd.norm():.2e}")


def section_4_grade_involution():
    """第四节：阶次对合"""
    print("\n" + "=" * 60)
    print("第四节：阶次对合 (Grade Involution)")
    print("=" * 60)

    alg = nblade.Algebra.euclidean(3)
    e1, e2, e3 = alg.basis_vectors()

    print("\n阶次对合（~ 操作符）对每个阶次乘以 (-1)^r:")
    print("  - r 为偶数时：保持不变")
    print("  - r 为奇数时：符号取反")
    print("\n数学定义: A* = (-1)^r A_r")

    # 创建各阶元素
    print("\n【各阶元素的阶次对合】")

    scalar = alg.scalar(5.0)
    print(f"0阶: {scalar}")
    print(f"  ~scalar = {~scalar}  (不变，因为 (-1)^0 = 1)")

    vector = e1
    print(f"\n1阶: {vector}")
    print(f"  ~vector = {~vector}  (变号，因为 (-1)^1 = -1)")

    bivector = e1 ^ e2
    print(f"\n2阶: {bivector}")
    print(f"  ~bivector = {~bivector}  (不变，因为 (-1)^2 = 1)")

    trivector = e1 ^ e2 ^ e3
    print(f"\n3阶: {trivector}")
    print(f"  ~trivector = {~trivector}  (变号，因为 (-1)^3 = -1)")

    # 混合多向量
    print("\n【混合多向量的阶次对合】")
    mixed = scalar + e1.scale(2.0) + bivector.scale(3.0) + trivector.scale(4.0)
    print(f"混合多向量: {mixed}")
    print(f"阶次对合后: {~mixed}")
    print("  奇阶（1, 3）变号，偶阶（0, 2）不变")

    # 验证偶部和奇部的关系
    print("\n【与偶部奇部的关系】")
    even = mixed.even_part()
    odd = mixed.odd_part()
    print(f"原多向量 = 偶部 + 奇部 = {even} + {odd}")
    print(f"阶次对合 = 偶部 - 奇部 = {even} - {odd}")
    print(f"验证: {even - odd}")
    print(f"~mixed = {~mixed}")

    # 双重应用
    print("\n【双重阶次对合】")
    print(f"~~mixed = {~~mixed}")
    print("  两次阶次对合恢复原始多向量")
    error = (~~mixed - mixed).norm()
    print(f"误差: {error:.2e}")

    # 与偶部奇部的关系
    print("\n【阶次对合与偶部奇部的关系】")
    print("  A = A_even + A_odd")
    print("  ~A = A_even - A_odd")
    print("  A + ~A = 2 * A_even")
    print("  A - ~A = 2 * A_odd")

    even_reconstructed = (mixed + ~mixed).scale(0.5)
    odd_reconstructed = (mixed - ~mixed).scale(0.5)
    print(f"\n从混合多向量重建偶部: {even_reconstructed}")
    print(f"直接提取偶部: {even}")
    print(f"从混合多向量重建奇部: {odd_reconstructed}")
    print(f"直接提取奇部: {odd}")


def section_5_reversion():
    """第五节：反序"""
    print("\n" + "=" * 60)
    print("第五节：反序 (Reversion)")
    print("=" * 60)

    alg = nblade.Algebra.euclidean(3)
    e1, e2, e3 = alg.basis_vectors()

    print("\n反序 (A†) 反转基向量的顺序：")
    print("  - 标量 (0阶): 不变")
    print("  - 向量 (1阶): 不变")
    print("  - 二重向量 (2阶): e1∧e2 → e2∧e1 = -e1∧e2")
    print("  - 三重向量 (3阶): e1∧e2∧e3 → e3∧e2∧e1 = e1∧e2∧e3")
    print("\n数学定义: A† = (-1)^(r(r-1)/2) A_r")

    # 各阶的反序
    print("\n【各阶元素的反序】")

    scalar = alg.scalar(5.0)
    print(f"0阶: {scalar}")
    print(f"  reversion = {scalar.reversion()}  (r=0: (-1)^0 = 1)")

    vector = e1.scale(2.0)
    print(f"\n1阶: {vector}")
    print(f"  reversion = {vector.reversion()}  (r=1: (-1)^0 = 1)")

    bivector = e1 ^ e2
    print(f"\n2阶: {bivector}")
    print(f"  reversion = {bivector.reversion()}  (r=2: (-1)^1 = -1)")

    trivector = e1 ^ e2 ^ e3
    print(f"\n3阶: {trivector}")
    print(f"  reversion = {trivector.reversion()}  (r=3: (-1)^3 = -1)")

    # 验证反序公式
    print("\n【反序公式验证】")
    B = e1 ^ e2
    print(f"二重向量 B = {B}")
    print(f"B† = {B.reversion()}")

    # e1*e2 = e1∧e2 + e1·e2 = e1∧e2 (正交)
    # (e1*e2)† = e2*e1 = -e1*e2
    product = e1 * e2
    product_rev = product.reversion()
    print(f"\n(e1*e2) = {product}")
    print(f"(e1*e2)† = {product_rev}")
    print(f"e2*e1 = {e2 * e1}")

    # 转子的反序
    print("\n【转子的反序】")
    plane = e1 ^ e2
    angle = math.pi / 4
    rotor = alg.rotor(plane, angle)
    rotor_rev = rotor.reversion()

    print(f"转子 R = {rotor}")
    print(f"R† = {rotor_rev}")

    # 验证 R·R† = 1
    print("\n【验证转子归一化】")
    product = rotor * rotor_rev
    print(f"R·R† = {product}")
    print(f"标量部分: {product.scalar_part():.6f} (应接近 1)")

    # 反序的应用：旋转公式
    print("\n【反序的应用：旋转公式】")
    v = e1
    rotated = v.rotate_by(rotor)
    print(f"使用 rotate_by: v' = {rotated}")

    rotated_manual = rotor * v * rotor.reversion()
    print(f"手动计算: R·v·R† = {rotated_manual}")
    print(f"误差: {(rotated - rotated_manual).norm():.2e}")

    # 双重反序
    print("\n【双重反序】")
    print(f"(B†)† = {(B.reversion().reversion())}")
    print(f"原二重向量 B = {B}")
    print("两次反序恢复原始多向量")


def section_6_summary():
    """第六节：总结"""
    print("\n" + "=" * 60)
    print("第六节：总结与对比")
    print("=" * 60)

    alg = nblade.Algebra.euclidean(3)
    e1, e2, e3 = alg.basis_vectors()

    # 创建一个包含所有阶次的混合多向量
    mixed = (
        alg.scalar(1.0)
        + e1.scale(2.0)
        + e2.scale(3.0)
        + (e1 ^ e2).scale(4.0)
        + (e2 ^ e3).scale(5.0)
        + (e1 ^ e2 ^ e3).scale(6.0)
    )

    print("\n【多向量操作对比】")
    print(f"原始多向量 A = {mixed}")

    print(f"\n阶次提取:")
    print(f"  grade(0): {mixed.grade(0)}")
    print(f"  grade(1): {mixed.grade(1)}")
    print(f"  grade(2): {mixed.grade(2)}")
    print(f"  grade(3): {mixed.grade(3)}")

    print(f"\n偶部和奇部:")
    print(f"  even_part(): {mixed.even_part()}")
    print(f"  odd_part(): {mixed.odd_part()}")

    print(f"\n对合操作:")
    print(f"  grade_involution (~A): {~mixed}")
    print(f"  reversion (A†): {mixed.reversion()}")

    # 表格对比
    print("\n【阶次对合和反序的符号变化】")
    print("-" * 50)
    print(f"{'阶次 r':<8} {'(-1)^r':<12} {'(-1)^(r(r-1)/2)':<20} {'名称'}")
    print("-" * 50)
    print(f"{'0':<8} {'+1':<12} {'+1':<20} {'标量'}")
    print(f"{'1':<8} {'-1':<12} {'+1':<20} {'向量'}")
    print(f"{'2':<8} {'+1':<12} {'-1':<20} {'二重向量'}")
    print(f"{'3':<8} {'-1':<12} {'-1':<20} {'三重向量'}")
    print(f"{'4':<8} {'+1':<12} {'+1':<20} {'四重向量'}")
    print("-" * 50)

    # 实际应用
    print("\n【实际应用】")
    print("1. 阶次提取：分离多向量的不同几何成分")
    print("2. 偶部/奇部：转子必须是偶元素；旋量可以是奇元素")
    print("3. 阶次对合：用于计算Clifford共轭")
    print("4. 反序：用于计算范数和逆；旋转公式 R·v·R†")


def main():
    """主函数"""
    print("=" * 60)
    print("nblade 阶次操作教程")
    print("理解多向量的阶次结构")
    print("=" * 60)

    section_1_what_are_grades()
    section_2_grade_extraction()
    section_3_even_odd_parts()
    section_4_grade_involution()
    section_5_reversion()
    section_6_summary()

    print("\n" + "=" * 60)
    print("教程完成！")
    print("=" * 60)
    print("\n核心要点:")
    print("  1. 阶次 = 基向量数量（0阶标量，1阶向量，2阶二重向量...）")
    print("  2. grade(r) 提取 r 阶部分")
    print("  3. even_part() 和 odd_part() 提取偶部和奇部")
    print("  4. ~A (阶次对合): 奇阶变号，偶阶不变")
    print("  5. A† (反序): 反转基向量顺序")
    print("\n下一步:")
    print("  - 运行 05_rotations.py 学习转子旋转")
    print("  - 运行 advanced/02_advanced_operations.py 学习高级操作")


if __name__ == "__main__":
    main()
