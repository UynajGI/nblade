# 计算机图形学示例概览

本文档介绍 nblade 在计算机图形学中的应用示例。

## 示例列表

### 01_transformations.py - 几何变换

**难度**: ⭐⭐⭐ 高级

**内容**:
- 反射变换
- 投影操作
- 旋转操作
- 缩放变换
- 组合变换
- 三角形变换示例
- 共形几何代数简介

**运行**:
```bash
python examples/cg/01_transformations.py
```

**预期输出**:
```
============================================================
nblade 计算机图形学变换示例
============================================================

============================================================
第1节：反射
============================================================

反射是最基础的几何变换
公式: v' = -n·v·n (n 为单位法向量)

【向量反射】
原始向量 v = 1.000000e1 + 2.000000e2 + 3.000000e3
反射平面法向量 n = 1.000000e1 (yz 平面)
反射结果 v' = 1.000000e1 + -2.000000e2 + -3.000000e3
...
```

---

### 02_projection.py - 投影操作

**难度**: ⭐⭐⭐ 高级

**内容**:
- 向量投影到另一个向量
- 子空间投影（投影到平面）
- 拒绝分量（正交分量）
- 实际应用：阴影计算

**运行**:
```bash
python examples/cg/02_projection.py
```

**预期输出**:
```
============================================================
nblade 投影操作示例
============================================================

============================================================
第一节：向量投影
============================================================

向量投影将一个向量分解到另一个向量的方向上
公式: proj_u(v) = (v·u / |u|²) · u

【基本向量投影】
向量 v = 3.000000*e1 + 4.000000*e2
投影方向 u = 1.000000*e1
投影结果 proj_u(v) = 3.000000*e1
  注意：y 分量被移除，只剩下 x 分量
...
```

---

### 03_reflection.py - 反射操作

**难度**: ⭐⭐⭐ 高级

**内容**:
- 向量反射（在超平面/向量法向）
- 平面反射
- 多次反射的组合 = 旋转
- 实际应用：镜面反射模拟

**运行**:
```bash
python examples/cg/03_reflection.py
```

**预期输出**:
```
============================================================
nblade 反射操作示例
============================================================

============================================================
第一节：向量反射（超平面反射）
============================================================

反射是最基本的几何变换，所有等距变换都可以由反射组合而成
在几何代数中，反射公式为: v' = -n·v·n
其中 n 是反射超平面的单位法向量

【基本反射】
原始向量 v = 1.000000*e1 + 2.000000*e2 + 3.000000*e3
反射面法向量 n = 1.000000*e1
反射结果 v' = -1.000000*e1 + 2.000000*e2 + 3.000000*e3
  注意：x 分量被取反
...
```

---

## 核心概念

### 反射

反射是最基础的几何变换。其他变换都可以通过反射的组合得到。

**公式**: `v' = -n v n`

其中 `n` 是反射平面的单位法向量。

```python
v = alg.vector([1, 2, 3])
n = e1  # yz 平面的法向量

reflected = v.reflect_in(n)
```

### 旋转

旋转使用转子实现。

**公式**: `v' = R v R†`

```python
import math

# 创建转子
plane = e1 ^ e2  # xy 平面
rotor = alg.rotor(plane, math.pi / 4)  # 45°

# 旋转
rotated = v.rotate_by(rotor)
```

### 投影

将向量投影到子空间。

**公式**: `proj_B(v) = (v⌋B)⌋B⁻¹`

```python
v = alg.vector([1, 2, 3])

# 投影到 x 轴
proj = v.project_to(e1)

# 正交分量
reject = v.reject_from(e1)

# v = proj + reject
```

### 拒绝分量

拒绝分量是向量中与投影方向正交的部分。

**公式**: `rej_B(v) = v - proj_B(v)`

```python
v = alg.vector([3.0, 4.0, 5.0])
direction = e1

proj = v.project_to(direction)
reject = v.reject_from(direction)

# proj 和 reject 正交
```

### 组合变换

多个变换通过几何积组合。

```python
# 先旋转，再反射
result = v.rotate_by(rotor).reflect_in(n)

# 变换顺序很重要！
result2 = v.reflect_in(n).rotate_by(rotor)  # 不同结果
```

---

## 与传统方法的对比

| 操作 | 传统方法 | 几何代数方法 |
|------|----------|--------------|
| 旋转 | 3×3 矩阵 | 转子（4 个参数） |
| 反射 | 3×3 矩阵 | 向量乘法 |
| 组合 | 矩阵乘法 | 几何积 |
| 插值 | 欧拉角/四元数 | 转子 SLERP |

**优势**:
- 参数更少（转子 vs 旋转矩阵）
- 数值稳定性更好
- 避免万向节锁
- 统一的运算规则

---

## 共形几何代数（CGA）

共形几何代数 G(4,1,0) 扩展了 3D 欧几里得代数，可以表示：

- 点、线、平面
- 圆、球体
- 刚体运动（包括平移）

```python
# 创建 CGA
cga = nblade.Algebra.cga()

# 在 CGA 中，平移也是转子
# 统一了旋转和平移
```

### CGA 的优势

1. **所有几何对象都是多向量**
2. **通过外积进行交集运算**
3. **平移和旋转统一为转子**
4. **公式简洁，避免特殊情况**

---

## 实用技巧

### 批量变换点集

```python
# 创建多个点
points = [e1, e2, e3, e1 + e2]

# 创建转子
rotor = alg.rotor(e1 ^ e2, math.pi / 4)

# 批量旋转
rotated_points = [p.rotate_by(rotor) for p in points]
```

### 从两点确定反射平面

```python
# 两点确定一条直线，该直线的垂直平面可用作反射面
p1 = alg.vector([1, 0, 0])
p2 = alg.vector([0, 1, 0])

line = p2 - p1
# 反射平面法向量需要外部指定
```

### 确定点在平面的哪一侧

```python
point = alg.vector([1, 2, 3])
plane_normal = e3  # xy 平面

# 点积的符号表示在哪一侧
side = (point | plane_normal).scalar_part()
if side > 0:
    print("在平面的正侧")
elif side < 0:
    print("在平面的负侧")
else:
    print("在平面上")
```

### 阴影计算

```python
# 正交投影到地面
ground_plane = e1 ^ e2  # xy 平面
shadow_point = point.project_to(ground_plane)

# 斜投影（考虑光源方向）
def oblique_projection(point, direction, plane_normal):
    point_dot_n = (point | plane_normal).scalar_part()
    dir_dot_n = (direction | plane_normal).scalar_part()
    t = -point_dot_n / dir_dot_n
    return point + direction.scale(t)
```

---

## 常见问题

### Q: 转子与四元数有什么区别？

A: 四元数是几何代数的偶子代数。转子是更一般的概念，可以在任意维度中表示旋转。在 3D 中，它们是等价的。

### Q: 如何实现平移？

A: 在标准几何代数中，平移通过向量加法实现。在共形几何代数中，平移也是转子，与旋转统一。

### Q: 如何进行旋转插值？

A: 使用 SLERP（球面线性插值）对转子进行插值：

```python
def slerp(R1, R2, t):
    # 简化版本
    R = R1.scale(1-t) + R2.scale(t)
    # 归一化...
    return R
```

---

## 参考

- [Geometric Algebra for Computer Science](https://www.elsevier.com/books/geometric-algebra-for-computer-science/dorst/978-0-12-374942-0) - Dorst, Fontijne, Mann
- [GA Wiki - 共形几何代数](https://en.wikipedia.org/wiki/Conformal_geometric_algebra)

---

## 下一步

- [教程示例](tutorials.md) - 基础概念
- [API 参考](../api/algebra.md) - 详细 API
