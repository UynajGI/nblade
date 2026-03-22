//! PyO3 Python 绑定模块 / PyO3 Python bindings module
//!
//! 提供几何代数的高性能 Python 接口
//! Provides high-performance Python interface for geometric algebra

use numpy::{PyArray1, PyArrayMethods};
use pyo3::prelude::*;
use pyo3::types::{PyList, PyModule, PyType};
use std::sync::Arc;

use crate::basis::frame as frame_ops;
use crate::geometry::rotation;
use crate::{AlgebraConfig, MultiVector, Signature};

/// Python 代数配置类 / Python algebra configuration class
///
/// 封装代数配置，用于定义几何代数的维度、签名等参数
/// Wraps algebra configuration to define dimension, signature, etc. for geometric algebra
#[pyclass(name = "AlgebraConfig")]
#[derive(Clone)]
pub struct PyAlgebraConfig {
    inner: Arc<AlgebraConfig>,
}

#[pymethods]
impl PyAlgebraConfig {
    /// 创建新的代数配置 / Create new algebra configuration
    ///
    /// # 参数 / Arguments
    /// * `dimension` - 向量空间维度 / Vector space dimension
    /// * `p` - 正平方基向量数量 (e_i^2 = +1) / Number of positive-square basis vectors
    /// * `q` - 负平方基向量数量 (e_i^2 = -1) / Number of negative-square basis vectors
    /// * `r` - 零平方基向量数量 (e_i^2 = 0) / Number of null (zero-square) basis vectors
    #[pyo3(text_signature = "(dimension, p, q, r)")]
    #[new]
    fn new(dimension: u32, p: u32, q: u32, r: u32) -> PyResult<Self> {
        let signature = Signature::new(p, q, r);
        let config = AlgebraConfig::new(dimension, signature);
        Ok(Self {
            inner: Arc::new(config),
        })
    }

    /// 创建欧几里得代数配置 / Create Euclidean algebra configuration
    ///
    /// 使用正交度量 (e_i^2 = +1) 创建 n 维欧几里得空间的几何代数
    /// Creates geometric algebra for n-dimensional Euclidean space with orthogonal metric (e_i^2 = +1)
    #[pyo3(text_signature = "(dimension)")]
    #[classmethod]
    fn euclidean(_cls: &Bound<'_, PyType>, dimension: u32) -> PyResult<Self> {
        let config = AlgebraConfig::euclidean(dimension);
        Ok(Self {
            inner: Arc::new(config),
        })
    }

    /// 获取维度 / Get dimension
    #[getter]
    fn dimension(&self) -> u32 {
        self.inner.dimension()
    }

    /// 获取基向量数量 (2^n) / Get number of basis vectors (2^n)
    #[getter]
    fn basis_count(&self) -> usize {
        self.inner.basis_count()
    }

    /// 获取度量签名 (p, q, r) / Get metric signature (p, q, r)
    ///
    /// p: 正平方基向量数, q: 负平方基向量数, r: 零平方基向量数
    /// p: positive-square, q: negative-square, r: null-square basis vectors
    #[getter]
    fn signature(&self) -> (u32, u32, u32) {
        let sig = &self.inner.signature;
        (sig.positive, sig.negative, sig.null)
    }

    fn __repr__(&self) -> String {
        format!(
            "AlgebraConfig(dimension={}, signature=G({},{},{}))",
            self.dimension(),
            self.inner.signature.positive,
            self.inner.signature.negative,
            self.inner.signature.null
        )
    }

    fn __str__(&self) -> String {
        format!(
            "G({},{},{}) Algebra ({}D, {} basis vectors)",
            self.inner.signature.positive,
            self.inner.signature.negative,
            self.inner.signature.null,
            self.dimension(),
            self.basis_count()
        )
    }

    /// 创建体积元素（伪标量）I = e₁∧e₂∧...∧eₙ
    ///
    /// Create volume element (pseudoscalar) I = e₁∧e₂∧...∧eₙ
    ///
    /// # 返回 / Returns
    /// 体积元素多向量 / Volume element multivector
    fn volume_element(&self) -> PyMultiVector {
        PyMultiVector {
            inner: self.inner.volume_element(),
        }
    }

    /// 计算体积元素的平方 I²
    ///
    /// Compute volume element squared I²
    ///
    /// # 返回 / Returns
    /// - `Some(value)` 如果签名非退化 / if signature is non-degenerate
    /// - `None` 如果签名退化（存在零向量）/ if signature is degenerate (null vectors exist)
    fn volume_element_squared(&self) -> Option<f64> {
        self.inner.volume_element_squared()
    }

    /// 创建体积元素的逆 I⁻¹
    ///
    /// Create volume element inverse I⁻¹
    ///
    /// # 返回 / Returns
    /// - `Some(mv)` 如果签名非退化 / if signature is non-degenerate
    /// - `None` 如果签名退化（存在零向量）/ if signature is degenerate (null vectors exist)
    fn volume_element_inverse(&self) -> Option<PyMultiVector> {
        self.inner
            .volume_element_inverse()
            .map(|mv| PyMultiVector { inner: mv })
    }
}

/// Python 多向量类 / Python multivector class
///
/// 几何代数的核心对象，可表示标量、向量、双向量及更高阶元素
/// Core object of geometric algebra, represents scalars, vectors, bivectors, and higher-grade elements
#[pyclass(name = "MultiVector")]
#[derive(Clone)]
pub struct PyMultiVector {
    inner: MultiVector,
}

/// Helper to extract Vec<PyMultiVector> from a Python list
fn extract_multivector_list(list: &Bound<'_, PyList>) -> PyResult<Vec<PyMultiVector>> {
    list.iter()
        .map(|item| item.extract::<PyMultiVector>())
        .collect()
}

#[pymethods]
impl PyMultiVector {
    /// 创建零多向量 / Create zero multivector
    ///
    /// 所有系数为零的多向量，作为加法的单位元
    /// Multivector with all zero coefficients, serves as additive identity
    #[pyo3(text_signature = "(config)")]
    #[classmethod]
    fn zeros(_cls: &Bound<'_, PyType>, config: &PyAlgebraConfig) -> PyResult<Self> {
        Ok(Self {
            inner: MultiVector::zeros(config.inner.clone()),
        })
    }

    /// 创建单位多向量 / Create unit multivector
    ///
    /// 标量部分为 1 的多向量，作为几何积的单位元
    /// Multivector with scalar part = 1, serves as multiplicative identity for geometric product
    #[pyo3(text_signature = "(config)")]
    #[classmethod]
    fn one(_cls: &Bound<'_, PyType>, config: &PyAlgebraConfig) -> PyResult<Self> {
        Ok(Self {
            inner: MultiVector::one(config.inner.clone()),
        })
    }

    /// 创建标量多向量 / Create scalar multivector
    ///
    /// 仅包含标量部分的多向量 (grade 0)
    /// Multivector containing only scalar part (grade 0)
    #[pyo3(text_signature = "(config, scalar)")]
    #[classmethod]
    fn from_scalar(
        _cls: &Bound<'_, PyType>,
        config: &PyAlgebraConfig,
        scalar: f64,
    ) -> PyResult<Self> {
        Ok(Self {
            inner: MultiVector::from_scalar(config.inner.clone(), scalar),
        })
    }

    /// 创建基向量 / Create basis vector
    ///
    /// 创建第 i 个标准基向量 e_i (i 从 0 开始)
    /// Creates the i-th standard basis vector e_i (0-indexed)
    #[pyo3(text_signature = "(config, i)")]
    #[classmethod]
    fn basis_vector(_cls: &Bound<'_, PyType>, config: &PyAlgebraConfig, i: u32) -> PyResult<Self> {
        Ok(Self {
            inner: MultiVector::basis_vector(config.inner.clone(), i),
        })
    }

    /// 从系数数组创建 / Create from coefficient array
    ///
    /// 所有系数显式指定，数组长度必须等于 2^n (n 为维度)
    /// All coefficients explicitly specified, array length must equal 2^n (n is dimension)
    #[pyo3(text_signature = "(config, coefficients)")]
    #[classmethod]
    fn from_coefficients(
        _cls: &Bound<'_, PyType>,
        config: &PyAlgebraConfig,
        coefficients: Vec<f64>,
    ) -> PyResult<Self> {
        Ok(Self {
            inner: MultiVector::from_coefficients(config.inner.clone(), coefficients),
        })
    }

    /// 获取代数配置 / Get algebra configuration
    #[getter]
    fn config(&self) -> PyAlgebraConfig {
        PyAlgebraConfig {
            inner: self.inner.config().clone(),
        }
    }

    /// 获取标量部分 / Get scalar part
    ///
    /// 返回 grade 0 部分的系数
    /// Returns the coefficient of grade 0 part
    fn scalar_part(&self) -> f64 {
        self.inner.scalar_part()
    }

    /// 获取系数数组 / Get coefficient array
    ///
    /// 对于密集多向量，使用零拷贝返回 NumPy 数组视图
    /// 对于稀疏多向量，返回新分配的 NumPy 数组
    /// For dense multivectors, returns zero-copy NumPy array view
    /// For sparse multivectors, returns newly allocated NumPy array
    fn coefficients<'py>(this: Bound<'py, Self>) -> Bound<'py, PyArray1<f64>> {
        let borrowed = this.borrow();
        let inner = &borrowed.inner;

        match inner {
            MultiVector::Dense(dense) => unsafe {
                PyArray1::borrow_from_array_bound(dense.as_array(), this.clone().into_any())
            },
            MultiVector::Sparse(_) => {
                let size = inner.config().basis_count();
                let mut coeffs = vec![0.0f64; size];
                for i in 0..size {
                    coeffs[i] = inner.get_coefficient(i);
                }
                PyArray1::from_vec_bound(this.py(), coeffs)
            }
        }
    }

    /// 从 NumPy 数组创建多向量 / Create multivector from NumPy array
    #[classmethod]
    fn from_numpy(
        _cls: &Bound<'_, PyType>,
        config: &PyAlgebraConfig,
        array: Bound<'_, PyArray1<f64>>,
    ) -> PyResult<Self> {
        let readonly = array.readonly();
        let slice = readonly.as_slice().map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyValueError, _>(format!(
                "Array is not contiguous: {:?}",
                e
            ))
        })?;
        let coefficients = slice.to_vec();
        Ok(Self {
            inner: MultiVector::from_coefficients(config.inner.clone(), coefficients),
        })
    }

    /// 获取阶次部分
    fn grade(&self, r: u32) -> Self {
        Self {
            inner: self.inner.grade_projection(r),
        }
    }

    /// 获取偶部
    fn even_part(&self) -> Self {
        Self {
            inner: self.inner.even_part(),
        }
    }

    /// 获取奇部
    fn odd_part(&self) -> Self {
        Self {
            inner: self.inner.odd_part(),
        }
    }

    /// 几何积
    fn geometric_product(&self, other: &PyMultiVector) -> PyResult<Self> {
        Ok(Self {
            inner: self.inner.geometric_product(&other.inner),
        })
    }

    /// 外积
    fn outer_product(&self, other: &PyMultiVector) -> PyResult<Self> {
        Ok(Self {
            inner: self.inner.outer_product(&other.inner),
        })
    }

    /// 左内积 A⌋B
    fn left_inner(&self, other: &PyMultiVector) -> PyResult<Self> {
        Ok(Self {
            inner: self.inner.left_inner(&other.inner),
        })
    }

    /// 右内积 A⌊B
    fn right_inner(&self, other: &PyMultiVector) -> PyResult<Self> {
        Ok(Self {
            inner: self.inner.right_inner(&other.inner),
        })
    }

    /// 阶次对合 A*
    fn grade_involution(&self) -> Self {
        Self {
            inner: self.inner.grade_involution(),
        }
    }

    /// 反转 A†
    fn reversion(&self) -> Self {
        Self {
            inner: self.inner.reversion(),
        }
    }

    /// Clifford 共轭 A‡
    fn clifford_conjugate(&self) -> Self {
        Self {
            inner: self.inner.clifford_conjugate(),
        }
    }

    /// 对偶 A⊥
    fn dual(&self) -> Self {
        Self {
            inner: self.inner.dual(),
        }
    }

    /// 逆对偶 A⁻⊥ / Inverse dual A⁻⊥
    ///
    /// # 公式 / Formula
    /// A⁻⊥ = A⌋I (左内积与体积元素) / A⁻⊥ = A⌋I (left inner product with volume element)
    fn inverse_dual(&self) -> Self {
        Self {
            inner: self.inner.inverse_dual(),
        }
    }

    /// 逆 A⁻¹
    fn inverse(&self) -> PyResult<Self> {
        self.inner
            .inverse()
            .map(|mv| Self { inner: mv })
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(e))
    }

    /// 标量积
    fn scalar_product(&self, other: &PyMultiVector) -> f64 {
        self.inner.scalar_product(&other.inner)
    }

    /// 范数平方
    fn norm_squared(&self) -> f64 {
        self.inner.norm_squared()
    }

    /// 范数
    fn norm(&self) -> f64 {
        self.inner.norm()
    }

    /// 交换子
    fn commutator(&self, other: &PyMultiVector) -> Self {
        Self {
            inner: self.inner.commutator(&other.inner),
        }
    }

    /// 正交投影
    fn project_to(&self, blade: &PyMultiVector) -> PyResult<Self> {
        self.inner
            .project_to(&blade.inner)
            .map(|mv| Self { inner: mv })
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(e))
    }

    /// 正交拒绝
    fn reject_from(&self, blade: &PyMultiVector) -> PyResult<Self> {
        self.inner
            .reject_from(&blade.inner)
            .map(|mv| Self { inner: mv })
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(e))
    }

    /// 反射
    fn reflect_in(&self, blade: &PyMultiVector) -> PyResult<Self> {
        self.inner
            .reflect_in(&blade.inner)
            .map(|mv| Self { inner: mv })
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(e))
    }

    /// 旋转
    fn rotate_by(&self, rotor: &PyMultiVector) -> PyResult<Self> {
        self.inner
            .rotate_by(&rotor.inner)
            .map(|mv| Self { inner: mv })
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(e))
    }

    /// 检查是否为零 / Check if zero
    ///
    /// 所有系数绝对值小于 EPSILON (1e-15) 时返回 True
    /// Returns True when all coefficients have absolute value less than EPSILON (1e-15)
    fn is_zero(&self) -> bool {
        self.inner.is_zero()
    }

    /// 检查是否可逆 / Check if invertible
    ///
    /// 范数平方大于 EPSILON (1e-15) 时返回 True
    /// Returns True when norm squared is greater than EPSILON (1e-15)
    fn is_invertible(&self) -> bool {
        self.inner.is_invertible()
    }

    /// 缩放 / Scale
    ///
    /// 将多向量的所有系数乘以标量
    /// Scale all coefficients by a scalar factor
    fn scale(&self, scalar: f64) -> Self {
        Self {
            inner: self.inner.scale(scalar),
        }
    }

    /// 加法 / Addition
    ///
    /// 将两个多向量相加
    /// Add two multivectors together
    fn add(&self, other: &PyMultiVector) -> Self {
        Self {
            inner: self.inner.add(&other.inner),
        }
    }

    // Python 魔术方法 / Python magic methods

    /// 加法 + / Addition +
    fn __add__(&self, other: &PyMultiVector) -> Self {
        Self {
            inner: self.inner.add(&other.inner),
        }
    }

    fn __sub__(&self, other: &PyMultiVector) -> Self {
        Self {
            inner: self.inner.sub(&other.inner),
        }
    }

    /// 几何积 * / Geometric product *
    fn __mul__(&self, other: &PyMultiVector) -> Self {
        Self {
            inner: self.inner.geometric_product(&other.inner),
        }
    }

    /// 反射乘法（标量 * 多向量）/ Reflected multiplication (scalar * multivector)
    fn __rmul__(&self, scalar: f64) -> Self {
        Self {
            inner: self.inner.scale(scalar),
        }
    }

    /// 外积 ^ / Outer product ^
    fn __xor__(&self, other: &PyMultiVector) -> Self {
        Self {
            inner: self.inner.outer_product(&other.inner),
        }
    }

    /// 左内积 | / Left inner product |
    fn __or__(&self, other: &PyMultiVector) -> Self {
        Self {
            inner: self.inner.left_inner(&other.inner),
        }
    }

    fn __invert__(&self) -> Self {
        Self {
            inner: self.inner.grade_involution(),
        }
    }

    fn __neg__(&self) -> Self {
        Self {
            inner: self.inner.neg(),
        }
    }

    fn __repr__(&self) -> String {
        format!("{}", self.inner)
    }

    fn __str__(&self) -> String {
        format!("{}", self.inner)
    }
}

// ============================================================================
// Module-level functions / 模块级函数
// ============================================================================

/// 创建转子 / Create a rotor
#[pyfunction]
fn create_rotor(plane: &PyMultiVector, angle: f64) -> PyResult<PyMultiVector> {
    Ok(PyMultiVector {
        inner: rotation::create_rotor(&plane.inner, angle),
    })
}

/// 计算互逆标架 / Compute reciprocal frame
///
/// 给定一个标架 {a₁, a₂, ..., aₙ}，计算互逆标架 {a¹, a², ..., aⁿ}
/// 使得 aⁱ⌋aⱼ = δⁱⱼ
///
/// Given a frame {a₁, a₂, ..., aₙ}, compute reciprocal frame {a¹, a², ..., aⁿ}
/// such that aⁱ⌋aⱼ = δⁱⱼ
///
/// # 参数 / Arguments
/// * `vectors` - 标架向量列表 / Frame vector list
///
/// # 返回 / Returns
/// 互逆标架向量列表 / Reciprocal frame vector list
///
/// # 异常 / Raises
/// ValueError - 标架为空、向量数量不匹配维度、体积元素不可逆等
#[pyfunction]
#[pyo3(name = "reciprocal_frame")]
fn reciprocal_frame_py(vectors: Bound<'_, PyList>) -> PyResult<Vec<PyMultiVector>> {
    let inner_vectors: Vec<MultiVector> = extract_multivector_list(&vectors)?
        .into_iter()
        .map(|v| v.inner)
        .collect();

    frame_ops::reciprocal_frame(&inner_vectors)
        .map(|recip| {
            recip
                .into_iter()
                .map(|mv| PyMultiVector { inner: mv })
                .collect()
        })
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(e))
}

/// 验证互逆标架条件 / Verify reciprocal frame condition
///
/// 检查 aⁱ⌋aⱼ = δⁱⱼ 是否成立
/// Verify that aⁱ⌋aⱼ = δⁱⱼ holds
///
/// # 参数 / Arguments
/// * `frame` - 原始标架 / Original frame
/// * `reciprocal` - 互逆标架 / Reciprocal frame
/// * `tolerance` - 容差（默认 1e-10）/ Tolerance (default 1e-10)
///
/// # 异常 / Raises
/// ValueError - 验证失败，包含错误详情
#[pyfunction]
#[pyo3(name = "verify_reciprocal_frame")]
#[pyo3(signature = (frame, reciprocal, tolerance=1e-10))]
fn verify_reciprocal_frame_py(
    frame: Bound<'_, PyList>,
    reciprocal: Bound<'_, PyList>,
    tolerance: f64,
) -> PyResult<()> {
    let inner_frame: Vec<MultiVector> = extract_multivector_list(&frame)?
        .into_iter()
        .map(|v| v.inner)
        .collect();
    let inner_reciprocal: Vec<MultiVector> = extract_multivector_list(&reciprocal)?
        .into_iter()
        .map(|v| v.inner)
        .collect();

    frame_ops::verify_reciprocal_frame(&inner_frame, &inner_reciprocal, tolerance)
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(e))
}

/// 计算标架的度量张量 / Compute metric tensor of frame
///
/// 度量张量 g_ij = a_i · a_j（标量积）
/// Metric tensor g_ij = a_i · a_j (scalar product)
///
/// # 参数 / Arguments
/// * `vectors` - 标架向量列表 / Frame vector list
///
/// # 返回 / Returns
/// 度量张量矩阵（n × n）/ Metric tensor matrix (n × n)
#[pyfunction]
#[pyo3(name = "metric_tensor")]
fn metric_tensor_py(vectors: Bound<'_, PyList>) -> PyResult<Vec<Vec<f64>>> {
    let inner_vectors: Vec<MultiVector> = extract_multivector_list(&vectors)?
        .into_iter()
        .map(|v| v.inner)
        .collect();

    frame_ops::metric_tensor(&inner_vectors)
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(e))
}

/// 多向量基展开 / Multivector basis expansion
///
/// 将多向量展开为基刃的线性组合：A = Σᴵ Aᴵ aᴵ
/// Expand multivector as linear combination of basis blades: A = Σᴵ Aᴵ aᴵ
///
/// # 参数 / Arguments
/// * `mv` - 待展开的多向量 / Multivector to expand
///
/// # 返回 / Returns
/// 基索引到系数的列表 [(index, coefficient), ...]
/// List of (basis index, coefficient) pairs
#[pyfunction]
#[pyo3(name = "basis_expansion")]
fn basis_expansion_py(mv: &PyMultiVector) -> Vec<(u64, f64)> {
    frame_ops::basis_expansion(&mv.inner)
}

/// 从展开系数重构多向量 / Reconstruct multivector from expansion
///
/// 给定基展开系数，重构原始多向量
/// Reconstruct original multivector from basis expansion coefficients
///
/// # 参数 / Arguments
/// * `config` - 代数配置 / Algebra configuration
/// * `coefficients` - 基索引到系数的列表 [(index, coefficient), ...]
///                    List of (basis index, coefficient) pairs
///
/// # 返回 / Returns
/// 重构的多向量 / Reconstructed multivector
#[pyfunction]
#[pyo3(name = "basis_reconstruction")]
fn basis_reconstruction_py(
    config: &PyAlgebraConfig,
    coefficients: Vec<(u64, f64)>,
) -> PyMultiVector {
    PyMultiVector {
        inner: frame_ops::basis_reconstruction(config.inner.clone(), &coefficients),
    }
}

/// Python 模块定义 / Python module definition
#[pymodule]
fn _core(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyAlgebraConfig>()?;
    m.add_class::<PyMultiVector>()?;

    // 添加模块级函数 / Add module-level functions
    m.add_function(wrap_pyfunction!(create_rotor, m)?)?;
    m.add_function(wrap_pyfunction!(reciprocal_frame_py, m)?)?;
    m.add_function(wrap_pyfunction!(verify_reciprocal_frame_py, m)?)?;
    m.add_function(wrap_pyfunction!(metric_tensor_py, m)?)?;
    m.add_function(wrap_pyfunction!(basis_expansion_py, m)?)?;
    m.add_function(wrap_pyfunction!(basis_reconstruction_py, m)?)?;

    // 添加常量 / Add constants
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;

    Ok(())
}
