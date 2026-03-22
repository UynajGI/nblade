# 物理示例概览

本页介绍 nblade 在物理学中的应用示例。

## 示例列表

### 01_rigid_body.py - 刚体动力学

**难度**: ⭐⭐⭐ 进阶

**内容**:
- 角动量作为二重向量
- 惯性张量
- 刚体旋转动力学
- 欧拉方程的几何代数形式
- 转动动能
- 简单刚体模拟

**运行**:
```bash
python examples/physics/01_rigid_body.py
```

**预期输出**:
```
============================================================
nblade 刚体物理示例
几何代数在经典力学中的应用
============================================================

============================================================
第一节：角动量作为二重向量
============================================================

在传统向量分析中，角动量 L = r × p
在几何代数中，角动量自然表示为二重向量 L = r ∧ p

位置向量 r = 3.000000e1 + 4.000000e2
动量向量 p = 2.000000e1 + -1.000000e2 + 5.000000e3

角动量二重向量 L = r∧p = -11.000000e1∧e2 + 15.000000e1∧e3 + 20.000000e2∧e3
...
```

---

## 核心概念

### 角动量的几何表示

在传统向量分析中，角动量定义为叉积：

$$\mathbf{L} = \mathbf{r} \times \mathbf{p}$$

在几何代数中，角动量自然表示为二重向量：

$$L = r \wedge p$$

这种表示的优势：
- 适用于任意维度（叉积只在 3D 有效）
- 明确的几何意义（张成的平面）
- 与其他几何代数运算一致

### 与传统叉积的关系

在 3D 中，角动量的对偶对应于传统叉积：

```python
L_bivector = r ^ p        # 几何代数表示
L_vector = L_bivector.dual()  # 对应传统叉积结果
```

### 刚体旋转

刚体姿态用转子描述：

```python
# 角速度平面
omega_plane = omega.dual()

# 更新转子
delta_rotor = alg.rotor(omega_plane, omega_magnitude * dt)
rotor = delta_rotor * rotor
```

---

## 理论背景

### 欧拉方程

无外力矩刚体的欧拉方程：

$$\frac{d\mathbf{L}}{dt} + \boldsymbol{\omega} \times \mathbf{L} = 0$$

在几何代数中：

$$\frac{dL}{dt} + \omega \wedge L = 0$$

### 转动动能

$$T = \frac{1}{2} \omega \cdot L = \frac{1}{2} \sum_i I_i \omega_i^2$$

在 nblade 中：

```python
T = 0.5 * (omega | L).scalar_part()
```

---

## 常见问题

### Q: 为什么用二重向量表示角动量？

A: 二重向量明确表示旋转平面，比传统向量（表示旋转轴）更具几何直观性。此外，二重向量在任意维度都有效。

### Q: 如何从角速度计算转子？

A: 角速度的对偶是旋转平面，然后使用该平面创建转子：

```python
omega = alg.vector([0, 0, 1])  # 绕 z 轴
omega_plane = omega.dual()      # e1∧e2
rotor = alg.rotor(omega_plane, angle)
```

---

## 参考资料

- [Geometric Algebra for Physicists](https://www.cambridge.org/core/books/geometric-algebra-for-physicists/) - Doran & Lasenby
- [New Foundations for Classical Mechanics](https://www.springer.com/gp/book/9780792355141) - David Hestenes

---

## 下一步

- [图形学示例](cg.md) - 几何变换
- [教程示例](tutorials.md) - 基础概念