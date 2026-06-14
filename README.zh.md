# nblade - N 维 Blade

[![PyPI](https://img.shields.io/pypi/v/nblade.svg)](https://pypi.org/project/nblade/)
[![Docs](https://readthedocs.org/projects/nblade/badge/?version=latest)](https://nblade.readthedocs.io/zh/latest/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**[English](../README.md)** | **[中文](README.zh.md)** | **[文档](https://nblade.readthedocs.io/zh/latest/)**

---

## 概述

**nblade**（N 维 Blade）是一个基于 Rust 实现的高性能几何代数库，提供 Python 绑定。支持任意维度（最高 64 维）和任意度量签名 G(p, q, r)。

## 特性

- **任意维度**：支持最高 64 维向量空间
- **任意签名**：支持 G(p, q, r) 度量签名（欧几里得、时空、共形等）
- **高性能**：Rust 后端，支持并行计算和 SIMD 优化
- **双表示**：自动选择密集或稀疏表示
- **完整运算**：所有标准几何代数运算
- **NumPy 集成**：直接从 NumPy 数组创建向量

## 安装

### pip 安装

```bash
pip install nblade
```

### 从源码安装

```bash
# 安装 maturin
pip install maturin

# 构建并安装
maturin develop --release
```

## 快速开始

```python
import nblade

# 创建 3D 欧几里得几何代数 G(3,0,0)
algebra = nblade.Algebra.euclidean(3)

# 创建基向量
e1 = algebra.basis_vector(0)
e2 = algebra.basis_vector(1)
e3 = algebra.basis_vector(2)

# 从列表创建向量
v = algebra.vector([1, 2, 3])

# 几何积
product = e1 * e2  # 结果为二重向量 e12

# 外积（楔积）
wedge = e1 ^ e2  # 结果为二重向量 e12

# 内积
inner = e1 | e2  # 结果为 0（正交向量）

# 对偶
I = algebra.config.volume_element()  # 伪标量
v_dual = v.dual()

# 使用转子旋转
import math
plane = e1 ^ e2
rotor = algebra.rotor(plane, math.pi / 2)  # 90 度旋转
rotated = e1.rotate_by(rotor)
```

## API 参考

### Algebra 类

```python
algebra = nblade.Algebra(dimension, p=0, q=0, r=0)

# 工厂方法
algebra = nblade.Algebra.euclidean(dimension)  # G(n, 0, 0)
algebra = nblade.Algebra.spacetime(dimension)  # G(1, n-1, 0)
algebra = nblade.Algebra.cga()                 # G(4, 1, 0)

# 属性
algebra.dimension   # 向量空间维度
algebra.signature   # (p, q, r) 元组

# 方法
algebra.basis_vector(i)      # 创建 e_i
algebra.vector([x, y, z])    # 从列表创建向量
algebra.scalar(value)        # 创建标量多向量
algebra.zeros()              # 创建零多向量
algebra.one()                # 创建单位标量
algebra.rotor(plane, angle)  # 创建转子
```

### MultiVector 类

```python
# 构造
mv = nblade.MultiVector.basis_vector(config, i)
mv = nblade.MultiVector.from_scalar(config, value)
mv = nblade.MultiVector.from_coefficients(config, coeffs)
mv = nblade.MultiVector.zeros(config)
mv = nblade.MultiVector.one(config)

# 乘积
mv.geometric_product(other)  # 几何积
mv.outer_product(other)      # 外积/楔积
mv.left_inner(other)         # 左收缩
mv.right_inner(other)        # 右收缩

# 对合
mv.grade_involution()        # 阶次对合 (A*)
mv.reversion()               # 反转 (A†)
mv.clifford_conjugate()      # Clifford 共轭 (A‡)

# 其他运算
mv.dual()                    # 对偶 (A⊥)
mv.inverse_dual()            # 逆对偶 (A⁻⊥)
mv.inverse()                 # 乘法逆
mv.norm()                    # 范数 |A|
mv.norm_squared()            # |A|²

# 阶次运算
mv.grade(r)                  # r 阶部分
mv.even_part()               # 偶部
mv.odd_part()                # 奇部

# 几何运算
mv.project_to(blade)         # 投影到子空间
mv.reject_from(blade)        # 拒绝分量
mv.reflect_in(blade)         # 反射
mv.rotate_by(rotor)          # 转子旋转

# 运算符重载
mv1 + mv2   # 加法
mv1 - mv2   # 减法
mv1 * mv2   # 几何积
mv1 ^ mv2   # 外积
mv1 | mv2   # 左内积
~mv        # 阶次对合
-mv        # 取负
```

### 模块函数

```python
# 创建旋转平面上的转子
rotor = nblade.create_rotor(plane, angle)

# 互逆标架计算
reciprocal = nblade.reciprocal_frame(vectors)

# 基展开
coeffs = nblade.basis_expansion(multivector)
```

## 支持的运算

| 运算 | 符号 | 公式 |
|-----------|--------|---------|
| 几何积 | `AB` | `AB = A·B + A∧B` |
| 外积 | `A∧B` | 反对称部分 |
| 左内积 | `A⌋B` | 左收缩 |
| 右内积 | `A⌊B` | 右收缩 |
| 阶对合 | `A*` | `(-1)^r A_r` |
| 反转 | `A†` | `(-1)^(r(r-1)/2) A_r` |
| Clifford 共轭 | `A‡` | `(A*)†` |
| 对偶 | `A⊥` | `A·I` 或 `AI` |
| 逆对偶 | `A⁻⊥` | `A⌋I` |

## 性能

nblade 针对性能进行了优化：

- **SIMD 加速**：2D-4D 运算（需 `simd` 特性）
- **自适应并行**：高维运算（≥6D）
- **内存池**：减少分配开销（需 `pool` 特性）
- **密集/稀疏自动选择**：基于系数密度

### 基准测试

| 运算 | 维度 | 时间 |
|-----------|-----------|------|
| 几何积 | 3D | ~5ms |
| 几何积 | 5D | ~15ms |
| 几何积 (SIMD) | 3D | ~3ms |

## 示例

### 2D 旋转

```python
import nblade
import math

# 创建 2D 欧几里得代数
algebra = nblade.Algebra.euclidean(2)

# 创建基向量
e1, e2 = algebra.basis_vectors()

# 创建 45 度旋转转子
plane = e1 ^ e2
rotor = algebra.rotor(plane, math.pi / 4)

# 旋转 e1
rotated = e1.rotate_by(rotor)
print(rotated)  # ~ 0.707*e1 + 0.707*e2
```

### 3D 叉积通过对偶

```python
import nblade

algebra = nblade.Algebra.euclidean(3)
e1, e2, e3 = algebra.basis_vectors()

# 叉积：a × b = (a ∧ b)*
# 其中 * 是 3D 对偶

a = e1 + 2*e2
b = 2*e1 + e3

# 计算叉积
a_cross_b = (a ^ b).dual()
print(a_cross_b)
```

## 贡献

欢迎贡献！请参阅 [CONTRIBUTING.md](CONTRIBUTING.md) 了解贡献指南。

### 开发环境设置

```bash
# 克隆仓库
git clone https://github.com/UynajGI/nblade.git
cd nblade

# 安装 Rust（如果尚未安装）
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh

# 构建
cargo build --release

# 运行测试
cargo test --all-features

# 构建 Python 绑定
pip install maturin
maturin develop --release

# 运行 Python 测试
pytest tests/test_python_api.py
```

### 代码风格

- Rust：遵循 `cargo fmt` 和 `cargo clippy`
- Python：遵循 PEP 8，使用 `ruff format`

## 许可证

MIT License
