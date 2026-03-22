# 对偶运算教程

本教程介绍几何代数中对偶的概念及其应用。

## 目录

1. [什么是对偶？](#1-什么是对偶)
2. [数学定义](#2-数学定义)
3. [对偶与叉积](#3-对偶与叉积)
4. [逆对偶](#4-逆对偶)
5. [应用](#5-应用)

---

## 1. 什么是对偶？

对偶是几何代数中的一个基本运算，它在多向量的不同阶次之间建立映射。在 3D 欧几里得空间中：

| 原始元素 | 对偶 |
|----------|------|
| 标量 | 三重向量（伪标量）|
| 向量 | 二重向量 |
| 二重向量 | 向量 |
| 三重向量 | 标量 |

对偶本质上给出了一个多向量的"正交补"。

### 直观理解

在 3D 中：
- **向量**的对偶是**二重向量**，表示垂直于该向量的平面
- **二重向量**的对偶是**向量**，垂直于该平面
- **标量**的对偶是**体积元素**（伪标量）

## 2. 数学定义

### 伪标量

伪标量 `I` 是代数中最高阶的元素：

```python
import nblade

alg = nblade.Algebra.euclidean(3)
e1, e2, e3 = alg.basis_vectors()

# 伪标量 I = e1 ∧ e2 ∧ e3
I = e1 ^ e2 ^ e3
print(f"I = {I}")  # e1∧e2∧e3
print(f"I² = {(I * I).scalar_part()}")  # 在 3D 欧几里得空间中为 -1
```

### 对偶运算

多向量 `A` 的对偶定义为：

```
A* = A · I⁻¹  (右收缩)
```

或等价地：

```
A* = A I⁻¹   (几何积)
```

### 在 nblade 中使用对偶

```python
import nblade

alg = nblade.Algebra.euclidean(3)
e1, e2, e3 = alg.basis_vectors()

# 向量的对偶
v = e1
v_dual = v.dual()
print(f"e1 的对偶: {v_dual}")  # e2∧e3

# 二重向量的对偶
B = e1 ^ e2
B_dual = B.dual()
print(f"e1∧e2 的对偶: {B_dual}")  # e3
```

## 3. 对偶与叉积

在 3D 中，叉积可以用对偶表示：

```
a × b = (a ∧ b)*
```

这就是为什么叉积只在 3D 中有效——它需要用对偶把二重向量映射回向量。

### 示例：通过对偶计算叉积

```python
import nblade

alg = nblade.Algebra.euclidean(3)
e1, e2, e3 = alg.basis_vectors()

a = e1
b = e2

# 通过对偶计算叉积: a × b = (a ∧ b)*
wedge = a ^ b        # a ∧ b = e1∧e2
cross = wedge.dual() # 对偶 = e3

print(f"a ∧ b = {wedge}")  # e1∧e2
print(f"(a ∧ b)* = {cross}")  # e3

# 这等价于叉积！
# e1 × e2 = e3
```

### 几何代数方法的优势

| 叉积 | 几何代数对偶 |
|---------------|------------------------|
| 只在 3D 中有效 | 任意维度都有效 |
| 返回向量 | 返回二重向量（更自然）|
| 非结合律 | 外积满足结合律 |

## 4. 逆对偶

逆对偶反转对偶运算：

```
(A*)* = A（在某些度量签名中可能差一个符号）
```

### 在 nblade 中使用逆对偶

```python
import nblade

alg = nblade.Algebra.euclidean(3)
e1, e2, e3 = alg.basis_vectors()

v = e1 + 2*e2 + 3*e3

# 对偶
v_dual = v.dual()

# 逆对偶（应该恢复原始向量）
v_recovered = v_dual.inverse_dual()

print(f"原始向量: {v}")
print(f"对偶: {v_dual}")
print(f"逆对偶: {v_recovered}")
```

## 5. 应用

### 从平面获取法向量

由二重向量 B 定义的平面，其法向量为 n = B*：

```python
import nblade

alg = nblade.Algebra.euclidean(3)
e1, e2, e3 = alg.basis_vectors()

# 由 e1∧e2 定义的平面（xy 平面）
plane = e1 ^ e2

# 法向量（指向 z 方向）
normal = plane.dual()
print(f"xy 平面的法向量: {normal}")  # e3
```

### 面积与体积

对偶自然地关联了面积和体积：

```python
import nblade

alg = nblade.Algebra.euclidean(3)
e1, e2, e3 = alg.basis_vectors()

# 面积元素
area = e1 ^ e2  # 表示面积的二重向量

# 面积的"体积"（大小）
# |B*| = |B| 在 3D 欧几里得空间中
area_magnitude = area.norm()
print(f"面积大小: {area_magnitude}")  # 1.0
```

### 电磁场

在物理学中，电磁场 F = E + iB 使用对偶来关联 E 场和 B 场：

```python
import nblade

alg = nblade.Algebra.euclidean(3)
e1, e2, e3 = alg.basis_vectors()
I = e1 ^ e2 ^ e3  # 伪标量

# 电场
E = e1 + e2

# 磁场作为二重向量（通过对偶）
# B_bivector = I * B_vector = dual(B_vector)
B_vector = e3
B_bivector = B_vector.dual()  # e1∧e2

print(f"磁场二重向量: {B_bivector}")
```

---

## 总结

| 概念 | 公式 | nblade 方法 |
|---------|---------|---------------|
| 对偶 | A* = AI⁻¹ | `A.dual()` |
| 逆对偶 | (A*)* | `A.inverse_dual()` |
| 叉积（3D）| a × b = (a ∧ b)* | `(a ^ b).dual()` |
| 伪标量 | I = e₁∧e₂∧...∧eₙ | `alg.config.volume_element()` |

## 延伸阅读

- 运行示例：`python examples/tutorials/04_dual_operations.py`
- 参见：[阶次操作](./03_grade.md)
- API 参考：[MultiVector.dual()](../api/multivector.md#dual)