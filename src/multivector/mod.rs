//! 多向量模块
//!
//! 提供密集和稀疏两种多向量表示

pub mod dense;
pub mod sparse;

#[cfg(feature = "pool")]
pub mod pool;

pub use dense::DenseMultiVector;
pub use sparse::SparseMultiVector;

#[cfg(feature = "pool")]
pub use pool::BufferPool;

use crate::algebra_config::AlgebraConfigRef;

/// 多向量枚举（自动选择密集或稀疏表示）
#[derive(Debug, Clone)]
pub enum MultiVector {
    /// 密集表示
    Dense(DenseMultiVector),
    /// 稀疏表示
    Sparse(SparseMultiVector),
}

impl MultiVector {
    /// 创建新的零多向量
    pub fn zeros(config: AlgebraConfigRef) -> Self {
        MultiVector::Dense(DenseMultiVector::zeros(config))
    }

    /// 创建单位多向量
    pub fn one(config: AlgebraConfigRef) -> Self {
        MultiVector::Dense(DenseMultiVector::one(config))
    }

    /// 从系数创建多向量（自动选择表示）
    pub fn from_coefficients(config: AlgebraConfigRef, coefficients: Vec<f64>) -> Self {
        let density = coefficients.iter().filter(|&&c| c.abs() > 1e-15).count() as f64
            / coefficients.len() as f64;

        if density < sparse::SPARSE_THRESHOLD {
            // 转换为稀疏表示
            let mut sparse = SparseMultiVector::zeros(config);
            for (i, &c) in coefficients.iter().enumerate() {
                if c.abs() > 1e-15 {
                    sparse.coefficients.insert(i as u64, c);
                }
            }
            MultiVector::Sparse(sparse)
        } else {
            MultiVector::Dense(DenseMultiVector::from_coefficients(config, coefficients))
        }
    }

    /// 从标量创建
    pub fn from_scalar(config: AlgebraConfigRef, scalar: f64) -> Self {
        MultiVector::Dense(DenseMultiVector::from_scalar(config, scalar))
    }

    /// 创建基向量
    pub fn basis_vector(config: AlgebraConfigRef, i: u32) -> Self {
        MultiVector::Dense(DenseMultiVector::basis_vector(config, i))
    }

    /// 获取标量部分
    pub fn scalar_part(&self) -> f64 {
        match self {
            MultiVector::Dense(d) => d.scalar_part(),
            MultiVector::Sparse(s) => s.scalar_part(),
        }
    }

    /// 获取系数
    pub fn get_coefficient(&self, index: usize) -> f64 {
        match self {
            MultiVector::Dense(d) => d.get_coefficient(index),
            MultiVector::Sparse(s) => s.get_coefficient(index as u64),
        }
    }

    /// 获取代数配置
    pub fn config(&self) -> &AlgebraConfigRef {
        match self {
            MultiVector::Dense(d) => &d.config,
            MultiVector::Sparse(s) => &s.config,
        }
    }

    /// 获取阶次部分
    pub fn grade_projection(&self, r: u32) -> Self {
        match self {
            MultiVector::Dense(d) => MultiVector::Dense(d.grade_projection(r)),
            MultiVector::Sparse(s) => MultiVector::Sparse(s.grade_projection(r)),
        }
    }

    /// 获取偶部
    pub fn even_part(&self) -> Self {
        match self {
            MultiVector::Dense(d) => MultiVector::Dense(d.even_part()),
            MultiVector::Sparse(s) => MultiVector::Sparse(s.even_part()),
        }
    }

    /// 获取奇部
    pub fn odd_part(&self) -> Self {
        match self {
            MultiVector::Dense(d) => MultiVector::Dense(d.odd_part()),
            MultiVector::Sparse(s) => MultiVector::Sparse(s.odd_part()),
        }
    }

    /// 加法
    pub fn add(&self, other: &Self) -> Self {
        match (self, other) {
            (MultiVector::Dense(a), MultiVector::Dense(b)) => MultiVector::Dense(a.add(b)),
            (MultiVector::Sparse(a), MultiVector::Sparse(b)) => MultiVector::Sparse(a.add(b)),
            (MultiVector::Dense(a), MultiVector::Sparse(b)) => {
                MultiVector::Dense(a.add(&b.to_dense()))
            }
            (MultiVector::Sparse(a), MultiVector::Dense(b)) => {
                MultiVector::Dense(a.to_dense().add(b))
            }
        }
    }

    /// 减法
    pub fn sub(&self, other: &Self) -> Self {
        match (self, other) {
            (MultiVector::Dense(a), MultiVector::Dense(b)) => MultiVector::Dense(a.sub(b)),
            (MultiVector::Sparse(a), MultiVector::Sparse(b)) => MultiVector::Sparse(a.sub(b)),
            (MultiVector::Dense(a), MultiVector::Sparse(b)) => {
                MultiVector::Dense(a.sub(&b.to_dense()))
            }
            (MultiVector::Sparse(a), MultiVector::Dense(b)) => {
                MultiVector::Dense(a.to_dense().sub(b))
            }
        }
    }

    /// 标量乘法
    pub fn scale(&self, scalar: f64) -> Self {
        match self {
            MultiVector::Dense(d) => MultiVector::Dense(d.scale(scalar)),
            MultiVector::Sparse(s) => MultiVector::Sparse(s.scale(scalar)),
        }
    }

    /// 取负
    pub fn neg(&self) -> Self {
        match self {
            MultiVector::Dense(d) => MultiVector::Dense(d.neg()),
            MultiVector::Sparse(s) => MultiVector::Sparse(s.neg()),
        }
    }

    /// 检查是否为零
    pub fn is_zero(&self) -> bool {
        match self {
            MultiVector::Dense(d) => d.is_zero(),
            MultiVector::Sparse(s) => s.is_zero(),
        }
    }

    /// 获取系数切片（仅密集表示有效）
    ///
    /// 对于密集多向量，直接返回内部数组的切片。
    /// 对于稀疏多向量，返回 None。
    #[inline]
    pub fn as_slice(&self) -> Option<&[f64]> {
        match self {
            MultiVector::Dense(d) => Some(d.as_slice()),
            MultiVector::Sparse(_) => None,
        }
    }

    /// 检查是否为密集表示
    #[inline]
    pub fn is_dense(&self) -> bool {
        matches!(self, MultiVector::Dense(_))
    }

    /// 转换为密集表示（如果需要则转换）
    pub fn to_dense(&self) -> Self {
        match self {
            MultiVector::Dense(_) => self.clone(),
            MultiVector::Sparse(s) => MultiVector::Dense(s.to_dense()),
        }
    }

    /// 获取内部 ndarray 数组引用（仅密集表示有效）
    pub fn as_dense_array(&self) -> Option<&ndarray::Array1<f64>> {
        match self {
            MultiVector::Dense(d) => Some(d.as_array()),
            MultiVector::Sparse(_) => None,
        }
    }
}

impl std::fmt::Display for MultiVector {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            MultiVector::Dense(d) => write!(f, "{}", d),
            MultiVector::Sparse(s) => write!(f, "{}", s),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::AlgebraConfig;
    use std::sync::Arc;

    fn test_config() -> AlgebraConfigRef {
        Arc::new(AlgebraConfig::euclidean(3))
    }

    fn test_config_4d() -> AlgebraConfigRef {
        Arc::new(AlgebraConfig::euclidean(4))
    }

    #[test]
    fn test_zeros() {
        let config = test_config();
        let mv = MultiVector::zeros(config.clone());
        assert!(mv.is_zero());
        assert!(mv.is_dense());
    }

    #[test]
    fn test_one() {
        let config = test_config();
        let mv = MultiVector::one(config.clone());
        assert!((mv.scalar_part() - 1.0).abs() < 1e-10);
        assert!(mv.grade_projection(1).is_zero());
    }

    #[test]
    fn test_from_scalar() {
        let config = test_config();
        let mv = MultiVector::from_scalar(config.clone(), 5.0);
        assert!((mv.scalar_part() - 5.0).abs() < 1e-10);
        assert!(mv.is_dense());
    }

    #[test]
    fn test_basis_vector() {
        let config = test_config();
        let e1 = MultiVector::basis_vector(config.clone(), 0);
        assert!((e1.get_coefficient(1) - 1.0).abs() < 1e-10);
        assert!(e1.is_dense());
    }

    #[test]
    fn test_from_coefficients_dense() {
        let config = test_config();
        let coeffs = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0];
        let mv = MultiVector::from_coefficients(config.clone(), coeffs);
        assert!(mv.is_dense());
    }

    #[test]
    fn test_from_coefficients_sparse() {
        let config = test_config_4d();
        let coeffs = vec![
            1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        ];
        let mv = MultiVector::from_coefficients(config.clone(), coeffs);
        assert!(!mv.is_dense());
    }

    #[test]
    fn test_scalar_part() {
        let config = test_config();
        let mv = MultiVector::from_scalar(config.clone(), 3.5);
        assert!((mv.scalar_part() - 3.5).abs() < 1e-10);

        let config_4d = test_config_4d();
        let sparse_coeffs = vec![
            2.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        ];
        let sparse_mv = MultiVector::from_coefficients(config_4d.clone(), sparse_coeffs);
        assert!((sparse_mv.scalar_part() - 2.5).abs() < 1e-10);
    }

    #[test]
    fn test_get_coefficient() {
        let config = test_config();
        let e1 = MultiVector::basis_vector(config.clone(), 0);
        assert!((e1.get_coefficient(1) - 1.0).abs() < 1e-10);
        assert!((e1.get_coefficient(0) - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_config_method() {
        let config = test_config();
        let mv = MultiVector::zeros(config.clone());
        assert_eq!(mv.config().dimension(), 3);
    }

    #[test]
    fn test_grade_projection_dense() {
        let config = test_config();
        let coeffs = vec![1.0, 2.0, 0.0, 3.0, 0.0, 0.0, 0.0, 0.0];
        let mv = MultiVector::from_coefficients(config.clone(), coeffs);

        let grade0 = mv.grade_projection(0);
        assert!((grade0.scalar_part() - 1.0).abs() < 1e-10);

        let grade1 = mv.grade_projection(1);
        assert!((grade1.get_coefficient(1) - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_grade_projection_sparse() {
        let config = test_config_4d();
        let coeffs = vec![
            1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        ];
        let mv = MultiVector::from_coefficients(config.clone(), coeffs);

        let grade0 = mv.grade_projection(0);
        assert!((grade0.scalar_part() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_even_part() {
        let config = test_config();
        let coeffs = vec![1.0, 2.0, 0.0, 3.0, 0.0, 0.0, 0.0, 0.0];
        let mv = MultiVector::from_coefficients(config.clone(), coeffs);

        let even = mv.even_part();
        assert!((even.get_coefficient(0) - 1.0).abs() < 1e-10);
        assert!((even.get_coefficient(3) - 3.0).abs() < 1e-10);
        assert!((even.get_coefficient(1) - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_odd_part() {
        let config = test_config();
        let coeffs = vec![1.0, 2.0, 0.0, 3.0, 0.0, 0.0, 0.0, 0.0];
        let mv = MultiVector::from_coefficients(config.clone(), coeffs);

        let odd = mv.odd_part();
        assert!((odd.get_coefficient(1) - 2.0).abs() < 1e-10);
        assert!((odd.get_coefficient(0) - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_add_dense_dense() {
        let config = test_config();
        let a = MultiVector::basis_vector(config.clone(), 0);
        let b = MultiVector::basis_vector(config.clone(), 1);

        let sum = a.add(&b);
        assert!((sum.get_coefficient(1) - 1.0).abs() < 1e-10);
        assert!((sum.get_coefficient(2) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_add_sparse_sparse() {
        let config = test_config_4d();
        let sparse_a = vec![
            1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        ];
        let sparse_b = vec![
            0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        ];
        let a = MultiVector::from_coefficients(config.clone(), sparse_a);
        let b = MultiVector::from_coefficients(config.clone(), sparse_b);

        let sum = a.add(&b);
        assert!((sum.get_coefficient(0) - 1.0).abs() < 1e-10);
        assert!((sum.get_coefficient(1) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_add_dense_sparse() {
        let config = test_config_4d();
        let dense = MultiVector::basis_vector(config.clone(), 0);
        let sparse_coeffs = vec![
            0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        ];
        let sparse = MultiVector::from_coefficients(config.clone(), sparse_coeffs);

        let sum = dense.add(&sparse);
        assert!((sum.get_coefficient(1) - 1.0).abs() < 1e-10);
        assert!((sum.get_coefficient(2) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_add_sparse_dense() {
        let config = test_config_4d();
        let sparse_coeffs = vec![
            0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        ];
        let sparse = MultiVector::from_coefficients(config.clone(), sparse_coeffs);
        let dense = MultiVector::basis_vector(config.clone(), 0);

        let sum = sparse.add(&dense);
        assert!((sum.get_coefficient(1) - 1.0).abs() < 1e-10);
        assert!((sum.get_coefficient(2) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_sub_dense_dense() {
        let config = test_config();
        let a = MultiVector::basis_vector(config.clone(), 0);
        let b = MultiVector::basis_vector(config.clone(), 1);

        let diff = a.sub(&b);
        assert!((diff.get_coefficient(1) - 1.0).abs() < 1e-10);
        assert!((diff.get_coefficient(2) - (-1.0)).abs() < 1e-10);
    }

    #[test]
    fn test_sub_sparse_sparse() {
        let config = test_config_4d();
        let sparse_a = vec![
            1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        ];
        let sparse_b = vec![
            0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        ];
        let a = MultiVector::from_coefficients(config.clone(), sparse_a);
        let b = MultiVector::from_coefficients(config.clone(), sparse_b);

        let diff = a.sub(&b);
        assert!((diff.get_coefficient(0) - 1.0).abs() < 1e-10);
        assert!((diff.get_coefficient(1) - (-1.0)).abs() < 1e-10);
    }

    #[test]
    fn test_sub_dense_sparse() {
        let config = test_config_4d();
        let dense = MultiVector::basis_vector(config.clone(), 0);
        let sparse_coeffs = vec![
            0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        ];
        let sparse = MultiVector::from_coefficients(config.clone(), sparse_coeffs);

        let diff = dense.sub(&sparse);
        assert!((diff.get_coefficient(1) - 1.0).abs() < 1e-10);
        assert!((diff.get_coefficient(2) - (-1.0)).abs() < 1e-10);
    }

    #[test]
    fn test_sub_sparse_dense() {
        let config = test_config_4d();
        let sparse_coeffs = vec![
            0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        ];
        let sparse = MultiVector::from_coefficients(config.clone(), sparse_coeffs);
        let dense = MultiVector::basis_vector(config.clone(), 0);

        let diff = sparse.sub(&dense);
        assert!((diff.get_coefficient(1) - (-1.0)).abs() < 1e-10);
        assert!((diff.get_coefficient(2) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_scale_dense() {
        let config = test_config();
        let mv = MultiVector::basis_vector(config.clone(), 0);
        let scaled = mv.scale(2.5);
        assert!((scaled.get_coefficient(1) - 2.5).abs() < 1e-10);
    }

    #[test]
    fn test_scale_sparse() {
        let config = test_config_4d();
        let sparse_coeffs = vec![
            1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        ];
        let mv = MultiVector::from_coefficients(config.clone(), sparse_coeffs);
        let scaled = mv.scale(3.0);
        assert!((scaled.scalar_part() - 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_neg_dense() {
        let config = test_config();
        let mv = MultiVector::from_scalar(config.clone(), 5.0);
        let negated = mv.neg();
        assert!((negated.scalar_part() - (-5.0)).abs() < 1e-10);
    }

    #[test]
    fn test_neg_sparse() {
        let config = test_config_4d();
        let sparse_coeffs = vec![
            2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        ];
        let mv = MultiVector::from_coefficients(config.clone(), sparse_coeffs);
        let negated = mv.neg();
        assert!((negated.scalar_part() - (-2.0)).abs() < 1e-10);
    }

    #[test]
    fn test_is_zero() {
        let config = test_config();
        let zero = MultiVector::zeros(config.clone());
        assert!(zero.is_zero());

        let non_zero = MultiVector::from_scalar(config.clone(), 1.0);
        assert!(!non_zero.is_zero());
    }

    #[test]
    fn test_as_slice_dense() {
        let config = test_config();
        let mv = MultiVector::from_scalar(config.clone(), 5.0);
        let slice = mv.as_slice();
        assert!(slice.is_some());
        let slice = slice.unwrap();
        assert_eq!(slice.len(), 8);
        assert!((slice[0] - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_as_slice_sparse() {
        let config = test_config_4d();
        let sparse_coeffs = vec![
            1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        ];
        let mv = MultiVector::from_coefficients(config.clone(), sparse_coeffs);
        let slice = mv.as_slice();
        assert!(slice.is_none());
    }

    #[test]
    fn test_is_dense() {
        let config = test_config();
        let dense = MultiVector::from_scalar(config.clone(), 1.0);
        assert!(dense.is_dense());

        let config_4d = test_config_4d();
        let sparse_coeffs = vec![
            1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        ];
        let sparse = MultiVector::from_coefficients(config_4d.clone(), sparse_coeffs);
        assert!(!sparse.is_dense());
    }

    #[test]
    fn test_to_dense_from_dense() {
        let config = test_config();
        let dense = MultiVector::from_scalar(config.clone(), 1.0);
        let still_dense = dense.to_dense();
        assert!(still_dense.is_dense());
        assert!((still_dense.scalar_part() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_to_dense_from_sparse() {
        let config = test_config_4d();
        let sparse_coeffs = vec![
            2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        ];
        let sparse = MultiVector::from_coefficients(config.clone(), sparse_coeffs);
        let dense = sparse.to_dense();
        assert!(dense.is_dense());
        assert!((dense.scalar_part() - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_as_dense_array_dense() {
        let config = test_config();
        let mv = MultiVector::from_scalar(config.clone(), 3.0);
        let array = mv.as_dense_array();
        assert!(array.is_some());
        let array = array.unwrap();
        assert_eq!(array.len(), 8);
        assert!((array[0] - 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_as_dense_array_sparse() {
        let config = test_config_4d();
        let sparse_coeffs = vec![
            1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        ];
        let mv = MultiVector::from_coefficients(config.clone(), sparse_coeffs);
        let array = mv.as_dense_array();
        assert!(array.is_none());
    }

    #[test]
    fn test_display_dense() {
        let config = test_config();
        let mv = MultiVector::from_scalar(config.clone(), 5.0);
        let display = format!("{}", mv);
        assert!(display.contains("5"));
    }

    #[test]
    fn test_display_sparse() {
        let config = test_config_4d();
        let sparse_coeffs = vec![
            1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        ];
        let mv = MultiVector::from_coefficients(config.clone(), sparse_coeffs);
        let display = format!("{}", mv);
        assert!(display.contains("1"));
    }

    #[test]
    fn test_display_zero() {
        let config = test_config();
        let mv = MultiVector::zeros(config.clone());
        let display = format!("{}", mv);
        assert_eq!(display, "0");
    }

    #[test]
    fn test_display_basis_vector() {
        let config = test_config();
        let e1 = MultiVector::basis_vector(config.clone(), 0);
        let display = format!("{}", e1);
        assert!(display.contains("e1"));
    }
}
