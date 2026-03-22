"""
从列表或数组创建几何代数向量的示例
"""

from nblade import Algebra
import numpy as np

# 创建 3D 欧几里得代数
alg = Algebra.euclidean(3)

# 方法 1: 使用 vector() 方法从列表创建向量
v1 = alg.vector([1, 2, 3])
print(f"从列表创建向量: {v1}")

# 方法 2: 使用 vector() 方法从 NumPy 数组创建向量
arr = np.array([4, 5, 6])
v2 = alg.vector(arr)
print(f"从 NumPy 数组创建向量: {v2}")

# 方法 3: 从列表创建更高维度的向量
alg_5d = Algebra.euclidean(5)
v5 = alg_5d.vector([1, 2, 3, 4, 5])
print(f"5D 向量: {v5}")

# 验证创建的向量确实是 1-向量（向量部分）
print(f"v1 的向量部分: {v1.grade(1)}")
print(f"v1 的标量部分: {v1.grade(0)}")  # 应该是 0

# 可以对创建的向量执行几何代数运算
result = v1 * v2  # 几何积
print(f"v1 * v2 = {result}")

wedge_result = v1 ^ v2  # 外积
print(f"v1 ^ v2 = {wedge_result}")

inner_result = v1 | v2  # 内积
print(f"v1 | v2 = {inner_result}")

# 也可以创建其他类型的多向量
# 从系数数组创建一般多向量
coeffs = [1.0, 2.0, 3.0, 0.0, 0.0, 4.0, 5.0, 6.0]  # 对应 3D 空间的所有部分
mv = alg.from_coefficients(coeffs)
print(f"从系数数组创建多向量: {mv}")

# 演示运算
print("\n运算示例:")
print(f"向量 v1 = {v1}")
print(f"向量 v2 = {v2}")
print(f"几何积 v1*v2 = {v1 * v2}")
print(f"外积 v1^v2 = {v1 ^ v2}")
print(f"左内积 v1|v2 = {v1 | v2}")
print(f"v1 的范数 = {v1.norm():.6f}")
print(f"v1 的平方范数 = {v1.norm_squared():.6f}")
