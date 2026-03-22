"""
几何代数基础运算示例

演示几何代数的基本运算：
- 向量创建和基本操作
- 几何积、外积、内积
- 对合运算
"""

from nblade import Algebra


def demo_basic_operations():
    """演示基本运算"""
    print("=== 基础几何代数运算示例 ===\n")

    # 创建 3D 欧几里得代数 G(3,0,0)
    alg = Algebra.euclidean(3)
    print(f"代数: {alg}")
    print(f"维度: {alg.dimension}")
    print(f"签名: {alg.signature}\n")

    # 获取基向量
    e1, e2, e3 = alg.basis_vectors()
    print("基向量:")
    print(f"  e1 = {e1}")
    print(f"  e2 = {e2}")
    print(f"  e3 = {e3}\n")

    # 创建一般向量
    v1 = alg.vector([1, 2, 3])
    v2 = alg.vector([4, 5, 6])
    print("一般向量:")
    print(f"  v1 = {v1}")
    print(f"  v2 = {v2}\n")

    # 基本运算
    print("基本运算:")

    # 几何积
    geom_prod = v1 * v2
    print(f"  几何积 v1*v2 = {geom_prod}")

    # 外积
    wedge_prod = v1 ^ v2
    print(f"  外积 v1^v2 = {wedge_prod}")

    # 左内积
    left_inner = v1 | v2
    print(f"  左内积 v1|v2 = {left_inner}")

    # 右内积
    right_inner = v1.right_inner(v2)
    print(f"  右内积 v1⌊v2 = {right_inner}\n")

    # 基向量的特殊性质
    print("基向量性质:")
    print(f"  e1² = {e1 * e1}")  # 应该是标量 1
    print(f"  e1·e2 = {e1 | e2}")  # 应该是标量 0 (正交)
    print(f"  e1^e2 = {e1 ^ e2}")  # 二重向量
    print(f"  (e1^e2)² = {(e1 ^ e2) * (e1 ^ e2)}\n")  # 二重向量的平方

    # 阶次投影
    print("阶次投影:")
    mv = alg.scalar(1.0) + v1 + (e1 ^ e2) + (e1 ^ e2 ^ e3)
    print(f"  多向量: {mv}")
    print(f"  标量部分: {mv.grade(0)}")
    print(f"  向量部分: {mv.grade(1)}")
    print(f"  二重向量部分: {mv.grade(2)}")
    print(f"  三重向量部分: {mv.grade(3)}\n")

    # 对合运算
    print("对合运算:")
    print(f"  原向量: {v1}")
    print(f"  阶次对合: {~v1}")  # ~ 表示阶次对合
    print(f"  反转: {v1.reversion()}")  # A†
    print(f"  Clifford 共轭: {v1.clifford_conjugate()}")  # A‡
    print(f"  偶部: {v1.even_part()}")
    print(f"  奇部: {v1.odd_part()}\n")

    print("基础运算示例完成！\n")


if __name__ == "__main__":
    demo_basic_operations()
