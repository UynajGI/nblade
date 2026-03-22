# MultiVector 类 API 参考

`MultiVector`（多重向量）是几何代数中的核心数据类型，表示多重向量（标量、向量、二重向量等的线性组合）。

## 类定义

```python
class MultiVector:
    """
    多重向量类，表示几何代数中的多重向量元素。
    可以是标量、向量、二重向量、三重向量等的线性组合。
    """
```

---

## 构造方法

### `zeros(config)`

创建一个零多重向量。

```python
@classmethod
def zeros(cls, config: AlgebraConfig) -> MultiVector
```

**参数**：`config` - 代数配置对象

**返回值**：值为零的多重向量

---

### `one(config)`

创建单位多重向量（标量 1）。

```python
@classmethod
def one(cls, config: AlgebraConfig) -> MultiVector
```

**返回值**：单位标量多重向量

---

### `from_scalar(config, scalar)`

从标量值创建多重向量。

```python
@classmethod
def from_scalar(cls, config: AlgebraConfig, scalar: float) -> MultiVector
```

**参数**：
- `config` - 代数配置
- `scalar` - 标量值

**返回值**：标量多重向量

**示例**：
```python
alg = nblade.Algebra.euclidean(3)
s = nblade.MultiVector.from_scalar(alg.config, 5.0)
# 或者使用便捷方法
s = alg.scalar(5.0)
```

---

### `basis_vector(config, i)`

创建基向量 e_i。

```python
@classmethod
def basis_vector(cls, config: AlgebraConfig, i: int) -> MultiVector
```

**参数**：
- `config` - 代数配置
- `i` - 基向量索引（从 0 开始）

**返回值**：基向量多重向量

---

### `from_coefficients(config, coefficients)`

从系数数组创建多重向量。

```python
@classmethod
def from_coefficients(cls, config: AlgebraConfig, coefficients: List[float]) -> MultiVector
```

**参数**：
- `config` - 代数配置
- `coefficients` - 所有基元素的系数数组

**返回值**：多重向量

**示例**：
```python
# 3D 中的系数顺序：[标量, e1, e2, e3, e12, e13, e23, e123]
coeffs = [1.0, 2.0, 3.0, 4.0, 0.5, 0.0, 0.0, 0.1]
mv = nblade.MultiVector.from_coefficients(alg.config, coeffs)
```

---

### `from_numpy(config, array)`

从 NumPy 数组创建多重向量。

```python
@classmethod
def from_numpy(cls, config: AlgebraConfig, array: numpy.ndarray) -> MultiVector
```

**参数**：
- `config` - 代数配置
- `array` - NumPy 数组

**返回值**：多重向量

---

## 属性

### `config`

```python
@property
def config(self) -> AlgebraConfig
```

返回关联的代数配置。

---

## 基本运算

### 几何积 `geometric_product()`

```python
def geometric_product(self, other: MultiVector) -> MultiVector
```

计算与另一个多重向量的几何积。

**运算符**：`*`

**示例**：
```python
a = alg.vector([1, 0, 0])
b = alg.vector([0, 1, 0])

# 方法调用
product = a.geometric_product(b)

# 或者使用运算符
product = a * b
```

---

### 外积 `outer_product()`

```python
def outer_product(self, other: MultiVector) -> MultiVector
```

计算外积（楔积）。

**运算符**：`^`

**示例**：
```python
a = alg.vector([1, 0, 0])
b = alg.vector([0, 1, 0])

wedge = a ^ b  # 二重向量
```

---

### 左内积 `left_inner()`

```python
def left_inner(self, other: MultiVector) -> MultiVector
```

计算左内积（左缩并）。

**运算符**：`|`

**示例**：
```python
a = alg.vector([1, 2, 3])
b = alg.vector([4, 5, 6])

inner = a | b  # 标量
```

---

### 右内积 `right_inner()`

```python
def right_inner(self, other: MultiVector) -> MultiVector
```

计算右内积（右缩并）。

**示例**：
```python
a = alg.vector([1, 2, 3])
b = alg.vector([4, 5, 6])

right = a.right_inner(b)
```

---

## 阶运算

### `grade(r)`

提取 r 阶部分。

```python
def grade(self, r: int) -> MultiVector
```

**参数**：`r` - 阶数（0=标量，1=向量，2=二重向量，...）

**返回值**：指定阶的多重向量

**示例**：
```python
# 混合多重向量
mv = alg.scalar(1.0) + e1 + (e1 ^ e2)

scalar = mv.grade(0)    # 标量部分
vector = mv.grade(1)    # 向量部分
bivector = mv.grade(2)  # 二重向量部分
```

---

### `even_part()`

提取偶数阶部分。

```python
def even_part(self) -> MultiVector
```

**返回值**：偶数阶分量（标量 + 二重向量 + 四重向量 + ...）

---

### `odd_part()`

提取奇数阶部分。

```python
def odd_part(self) -> MultiVector
```

**返回值**：奇数阶分量（向量 + 三重向量 + ...）

---

### `scalar_part()`

提取标量部分。

```python
def scalar_part(self) -> float
```

**返回值**：标量值

**示例**：
```python
a = alg.vector([3, 4, 0])
b = alg.vector([3, 4, 0])

inner = a | b
print(inner.scalar_part())  # 25
```

---

### `coefficients()`

获取所有基元素的系数。

```python
def coefficients(self) -> List[float]
```

**返回值**：系数列表

---

## 对合运算

### `grade_involution()`

阶对合。

```python
def grade_involution(self) -> MultiVector
```

将 r 阶元素乘以 (-1)^r。

**运算符**：`~`（注意：在 Python 中 ~ 是按位非运算符，因此请使用方法调用）

**示例**：
```python
v = alg.vector([1, 2, 3])
involution = v.grade_involution()  # -v（向量是 1 阶）
```

---

### `reversion()`

反转。

```python
def reversion(self) -> MultiVector
```

将 r 阶元素乘以 (-1)^(r(r-1)/2)。

**示例**：
```python
v = alg.vector([1, 2, 3])
rev = v.reversion()  # v（向量保持不变）

B = e1 ^ e2
B_rev = B.reversion()  # -B（二重向量符号翻转）
```

---

### `clifford_conjugate()`

克利福德共轭。

```python
def clifford_conjugate(self) -> MultiVector
```

等于阶对合加反转。

**示例**：
```python
mv = alg.scalar(1) + e1 + (e1 ^ e2)
conj = mv.clifford_conjugate()
```

---

## 对偶运算

### `dual()`

霍奇对偶。

```python
def dual(self) -> MultiVector
```

将 k 重向量映射到 (n-k) 重向量。

**示例**：
```python
v = alg.vector([1, 2, 3])
dual_v = v.dual()  # 二重向量

# 在 3D 中，对偶建立了叉积与外积之间的联系
# a x b = (a ^ b).dual()
```

---

### `inverse_dual()`

逆对偶。

```python
def inverse_dual(self) -> MultiVector
```

**示例**：
```python
v = alg.vector([1, 2, 3])
original = v.dual().inverse_dual()  # 等于 v
```

---

## 范数与逆元

### `norm()`

计算范数。

```python
def norm(self) -> float
```

**返回值**：|A| = sqrt(<A†A>₀)

**示例**：
```python
v = alg.vector([3, 4, 0])
print(v.norm())  # 5.0
```

---

### `norm_squared()`

计算范数的平方。

```python
def norm_squared(self) -> float
```

**返回值**：|A|²

---

### `inverse()`

计算乘法逆元。

```python
def inverse(self) -> MultiVector
```

**返回值**：A⁻¹ 满足 A * A⁻¹ = 1

**异常**：如果不可逆则抛出异常

**示例**：
```python
v = alg.vector([1, 2, 3])
v_inv = v.inverse()

# 验证
product = v * v_inv
print(product.scalar_part())  # 约等于 1
```

---

### `is_invertible()`

检查是否可逆。

```python
def is_invertible(self) -> bool
```

---

## 几何运算

### `project_to(blade)`

投影到由 blade 表示的子空间。

```python
def project_to(self, blade: MultiVector) -> MultiVector
```

**参数**：`blade` - 表示目标子空间的 blade

**返回值**：投影后的多重向量

**示例**：
```python
v = alg.vector([1, 2, 3])
e1 = alg.basis_vector(0)

proj = v.project_to(e1)  # 投影到 x 轴
print(proj)  # e1
```

---

### `reject_from(blade)`

从由 blade 表示的子空间中拒绝（即求垂直分量）。

```python
def reject_from(self, blade: MultiVector) -> MultiVector
```

**返回值**：正交分量

**示例**：
```python
v = alg.vector([1, 2, 3])
e1 = alg.basis_vector(0)

reject = v.reject_from(e1)  # 垂直于 x 轴的分量
# v = proj + reject
```

---

### `reflect_in(blade)`

在由 blade 表示的超平面中反射。

```python
def reflect_in(self, blade: MultiVector) -> MultiVector
```

**公式**：v' = -n v n（其中 n 是法向量）

**示例**：
```python
v = alg.vector([1, 2, 3])
e1 = alg.basis_vector(0)

reflected = v.reflect_in(e1)  # 在 yz 平面中反射
# x 分量被取反
```

---

### `rotate_by(rotor)`

使用旋量进行旋转。

```python
def rotate_by(self, rotor: MultiVector) -> MultiVector
```

**参数**：`rotor` - 旋转旋量

**返回值**：旋转后的多重向量

**公式**：v' = R v R†

**示例**：
```python
import math

alg = nblade.Algebra.euclidean(3)
e1, e2, e3 = alg.basis_vectors()

# 创建旋量
plane = e1 ^ e2
rotor = alg.rotor(plane, math.pi / 4)

# 旋转
rotated = e1.rotate_by(rotor)
```

---

## 其他运算

### `scale(factor)`

标量乘法。

```python
def scale(self, factor: float) -> MultiVector
```

**参数**：`factor` - 缩放因子

**返回值**：缩放后的多重向量

---

### `add(other)`

加法。

```python
def add(self, other: MultiVector) -> MultiVector
```

**运算符**：`+`

---

### `commutator(other)`

交换子 [A, B] = (AB - BA) / 2。

```python
def commutator(self, other: MultiVector) -> MultiVector
```

---

### `scalar_product(other)`

标量积 <AB>₀。

```python
def scalar_product(self, other: MultiVector) -> float
```

---

### `is_zero()`

检查是否为零。

```python
def is_zero(self) -> bool
```

---

## 运算符重载

| 运算符 | 方法 | 描述 |
|--------|------|------|
| `a + b` | `add` | 加法 |
| `a - b` | `__sub__` | 减法 |
| `a * b` | `geometric_product` | 几何积 |
| `a ^ b` | `outer_product` | 外积 |
| `a \| b` | `left_inner` | 左内积 |
| `-a` | `__neg__` | 取负 |
| `~a` | `grade_involution` | 阶对合 |

---

## 完整示例

```python
import nblade
import math

# 创建代数
alg = nblade.Algebra.euclidean(3)
e1, e2, e3 = alg.basis_vectors()

# 创建向量
v = alg.vector([1.0, 2.0, 3.0])
w = alg.vector([4.0, 5.0, 6.0])

# 基本运算
print("几何积:", v * w)
print("外积:", v ^ w)
print("内积:", v | w)

# 阶运算
mixed = alg.scalar(1.0) + v + (e1 ^ e2)
print("标量部分:", mixed.grade(0))
print("向量部分:", mixed.grade(1))

# 旋转
rotor = alg.rotor(e1 ^ e2, math.pi / 4)
rotated = v.rotate_by(rotor)
print("旋转后:", rotated)

# 范数与逆元
print("范数:", v.norm())
print("逆元:", v.inverse())

# 投影与拒绝
proj = v.project_to(e1)
reject = v.reject_from(e1)
print("投影:", proj)
print("拒绝:", reject)
```
