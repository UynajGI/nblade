# nblade 示例目录

本目录包含几何代数库的各种示例，按功能分类组织：

## 目录结构

```
examples/
├── tutorials/       # 教程示例（推荐入门）
│   ├── 01_quickstart.py        # 快速入门
│   ├── 02_basic_operations.py  # 基本运算
│   ├── 03_vectors.py           # 向量操作
│   ├── 04_dual_operations.py   # 对偶运算
│   ├── 05_rotations.py         # 旋转操作
│   ├── 06_reciprocal_frame.py  # 互逆标架
│   └── 07_grade_operations.py  # 阶次操作
├── basic/           # 基础运算示例
│   ├── 01_basic_operations.py  # 基本运算演示
│   ├── 02_vector_creation.py   # 向量创建演示
│   └── vector_from_array.py    # 从数组创建向量
├── advanced/        # 高级运算示例
│   ├── 02_advanced_operations.py
│   ├── 03_spacetime.py         # 时空代数
│   └── 04_conformal.py         # 共形几何代数
├── applications/    # 应用示例
│   └── 03_applications.py      # 几何应用演示
├── physics/         # 物理应用示例
│   ├── 01_rigid_body.py        # 刚体动力学
│   ├── 02_kepler_problem.py    # 开普勒问题
│   └── 03_electromagnetism.py  # 电磁学
└── cg/              # 计算机图形学示例
    ├── 01_transformations.py   # 几何变换
    ├── 02_projection.py        # 投影操作
    ├── 03_reflection.py        # 反射操作
    └── 04_interpolation.py     # 转子插值
```

## 示例说明

### tutorials/ (推荐入门)
- **01_quickstart.py**: 快速入门，5 分钟上手 nblade
- **02_basic_operations.py**: 几何积、外积、内积等基本运算
- **03_vectors.py**: 向量创建和操作
- **04_dual_operations.py**: Hodge 对偶、与叉积关系
- **05_rotations.py**: 使用转子进行旋转
- **06_reciprocal_frame.py**: 互逆标架、度量张量、向量分解
- **07_grade_operations.py**: 阶次提取、偶部奇部

### basic/
- **01_basic_operations.py**: 演示基本几何代数运算
- **02_vector_creation.py**: 演示从 Python 列表/数组创建向量

### advanced/
- **02_advanced_operations.py**: 演示高级运算（对合、逆、对偶、交换子等）
- **03_spacetime.py**: 时空代数 G(1,3)，四向量、洛伦兹变换
- **04_conformal.py**: 共形几何代数 G(4,1)，点、圆、球、变换

### applications/
- **03_applications.py**: 演示实际应用（旋转、反射、投影等几何变换）

### physics/
- **01_rigid_body.py**: 刚体动力学模拟示例
- **02_kepler_problem.py**: 开普勒问题、角动量、轨道方程
- **03_electromagnetism.py**: 电磁场作为二重向量、麦克斯韦方程

### cg/
- **01_transformations.py**: 计算机图形学中的几何变换
- **02_projection.py**: 向量投影、子空间投影、阴影计算
- **03_reflection.py**: 向量反射、平面反射、镜面模拟
- **04_interpolation.py**: 转子 SLERP 插值、平滑动画

## 运行示例

```bash
# 运行教程示例
python examples/tutorials/01_quickstart.py

# 运行基础示例
python examples/basic/01_basic_operations.py

# 运行高级示例
python examples/advanced/03_spacetime.py
python examples/advanced/04_conformal.py

# 运行物理示例
python examples/physics/01_rigid_body.py
python examples/physics/02_kepler_problem.py

# 运行 CG 示例
python examples/cg/04_interpolation.py
```

## 依赖

所有示例都需要安装 nblade 库：

```bash
pip install maturin
maturin develop  # 在项目根目录运行
```