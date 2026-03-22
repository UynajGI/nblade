//! 内积模块
//!
//! 实现左内积 (A⌋B) 和右内积 (A⌊B)
//!
//! 左内积公式：A⌋B = Σ⟨ArBs⟩_{s-r}
//! 右内积公式：A⌊B = Σ⟨ArBs⟩_{r-s}

use crate::basis::index::{basis_left_inner, basis_right_inner, BasisIndex};
use crate::multivector::{DenseMultiVector, MultiVector, SparseMultiVector};
use crate::products::should_use_parallel;
use rayon::prelude::*;

/// 密集多向量的左内积（自适应并行调度）
pub fn left_inner_dense(a: &DenseMultiVector, b: &DenseMultiVector) -> DenseMultiVector {
    let config = a.config.clone();
    let size = config.basis_count();

    let non_zero_a: Vec<_> = a
        .coefficients
        .iter()
        .enumerate()
        .filter(|(_, &c)| c.abs() > 1e-15)
        .collect();

    let non_zero_b: Vec<_> = b
        .coefficients
        .iter()
        .enumerate()
        .filter(|(_, &c)| c.abs() > 1e-15)
        .collect();

    if should_use_parallel(size) {
        left_inner_dense_parallel(&non_zero_a, &non_zero_b, config, size)
    } else {
        left_inner_dense_sequential(&non_zero_a, &non_zero_b, config, size)
    }
}

fn left_inner_dense_sequential(
    non_zero_a: &[(usize, &f64)],
    non_zero_b: &[(usize, &f64)],
    config: crate::algebra_config::AlgebraConfigRef,
    size: usize,
) -> DenseMultiVector {
    let mut result = vec![0.0f64; size];

    for (i, &coef_a) in non_zero_a {
        for (j, &coef_b) in non_zero_b {
            if let Some((result_idx, sign)) = basis_left_inner(*i as BasisIndex, *j as BasisIndex) {
                result[result_idx as usize] += sign as f64 * coef_a * coef_b;
            }
        }
    }

    DenseMultiVector::from_coefficients(config, result)
}

fn left_inner_dense_parallel(
    non_zero_a: &[(usize, &f64)],
    non_zero_b: &[(usize, &f64)],
    config: crate::algebra_config::AlgebraConfigRef,
    size: usize,
) -> DenseMultiVector {
    let result: Vec<f64> = non_zero_a
        .par_iter()
        .fold(
            || vec![0.0f64; size],
            |mut local_result, (i, &coef_a)| {
                for (j, &coef_b) in non_zero_b {
                    if let Some((result_idx, sign)) =
                        basis_left_inner(*i as BasisIndex, *j as BasisIndex)
                    {
                        local_result[result_idx as usize] += sign as f64 * coef_a * coef_b;
                    }
                }
                local_result
            },
        )
        .reduce(
            || vec![0.0f64; size],
            |mut a, b| {
                for i in 0..size {
                    a[i] += b[i];
                }
                a
            },
        );

    DenseMultiVector::from_coefficients(config, result)
}

/// 稀疏多向量的左内积
pub fn left_inner_sparse(a: &SparseMultiVector, b: &SparseMultiVector) -> SparseMultiVector {
    let mut result_coeffs = std::collections::HashMap::new();

    for (&i, &coef_a) in &a.coefficients {
        for (&j, &coef_b) in &b.coefficients {
            if let Some((result, sign)) = basis_left_inner(i, j) {
                let product = sign as f64 * coef_a * coef_b;
                *result_coeffs.entry(result).or_insert(0.0) += product;
            }
        }
    }

    result_coeffs.retain(|_, &mut v| v.abs() > 1e-15);

    SparseMultiVector::from_map(a.config.clone(), result_coeffs)
}

/// 密集多向量的右内积（自适应并行调度）
pub fn right_inner_dense(a: &DenseMultiVector, b: &DenseMultiVector) -> DenseMultiVector {
    let config = a.config.clone();
    let size = config.basis_count();

    let non_zero_a: Vec<_> = a
        .coefficients
        .iter()
        .enumerate()
        .filter(|(_, &c)| c.abs() > 1e-15)
        .collect();

    let non_zero_b: Vec<_> = b
        .coefficients
        .iter()
        .enumerate()
        .filter(|(_, &c)| c.abs() > 1e-15)
        .collect();

    if should_use_parallel(size) {
        right_inner_dense_parallel(&non_zero_a, &non_zero_b, config, size)
    } else {
        right_inner_dense_sequential(&non_zero_a, &non_zero_b, config, size)
    }
}

fn right_inner_dense_sequential(
    non_zero_a: &[(usize, &f64)],
    non_zero_b: &[(usize, &f64)],
    config: crate::algebra_config::AlgebraConfigRef,
    size: usize,
) -> DenseMultiVector {
    let mut result = vec![0.0f64; size];

    for (i, &coef_a) in non_zero_a {
        for (j, &coef_b) in non_zero_b {
            if let Some((result_idx, sign)) = basis_right_inner(*i as BasisIndex, *j as BasisIndex)
            {
                result[result_idx as usize] += sign as f64 * coef_a * coef_b;
            }
        }
    }

    DenseMultiVector::from_coefficients(config, result)
}

fn right_inner_dense_parallel(
    non_zero_a: &[(usize, &f64)],
    non_zero_b: &[(usize, &f64)],
    config: crate::algebra_config::AlgebraConfigRef,
    size: usize,
) -> DenseMultiVector {
    let result: Vec<f64> = non_zero_a
        .par_iter()
        .fold(
            || vec![0.0f64; size],
            |mut local_result, (i, &coef_a)| {
                for (j, &coef_b) in non_zero_b {
                    if let Some((result_idx, sign)) =
                        basis_right_inner(*i as BasisIndex, *j as BasisIndex)
                    {
                        local_result[result_idx as usize] += sign as f64 * coef_a * coef_b;
                    }
                }
                local_result
            },
        )
        .reduce(
            || vec![0.0f64; size],
            |mut a, b| {
                for i in 0..size {
                    a[i] += b[i];
                }
                a
            },
        );

    DenseMultiVector::from_coefficients(config, result)
}

/// 稀疏多向量的右内积
pub fn right_inner_sparse(a: &SparseMultiVector, b: &SparseMultiVector) -> SparseMultiVector {
    let mut result_coeffs = std::collections::HashMap::new();

    for (&i, &coef_a) in &a.coefficients {
        for (&j, &coef_b) in &b.coefficients {
            if let Some((result, sign)) = basis_right_inner(i, j) {
                let product = sign as f64 * coef_a * coef_b;
                *result_coeffs.entry(result).or_insert(0.0) += product;
            }
        }
    }

    result_coeffs.retain(|_, &mut v| v.abs() > 1e-15);

    SparseMultiVector::from_map(a.config.clone(), result_coeffs)
}

impl MultiVector {
    /// 左内积 A⌋B
    pub fn left_inner(&self, other: &Self) -> Self {
        match (self, other) {
            (MultiVector::Dense(a), MultiVector::Dense(b)) => {
                MultiVector::Dense(left_inner_dense(a, b))
            }
            (MultiVector::Sparse(a), MultiVector::Sparse(b)) => {
                MultiVector::Sparse(left_inner_sparse(a, b))
            }
            (MultiVector::Dense(a), MultiVector::Sparse(b)) => {
                MultiVector::Dense(left_inner_dense(a, &b.to_dense()))
            }
            (MultiVector::Sparse(a), MultiVector::Dense(b)) => {
                MultiVector::Dense(left_inner_dense(&a.to_dense(), b))
            }
        }
    }

    /// 右内积 A⌊B
    pub fn right_inner(&self, other: &Self) -> Self {
        match (self, other) {
            (MultiVector::Dense(a), MultiVector::Dense(b)) => {
                MultiVector::Dense(right_inner_dense(a, b))
            }
            (MultiVector::Sparse(a), MultiVector::Sparse(b)) => {
                MultiVector::Sparse(right_inner_sparse(a, b))
            }
            (MultiVector::Dense(a), MultiVector::Sparse(b)) => {
                MultiVector::Dense(right_inner_dense(a, &b.to_dense()))
            }
            (MultiVector::Sparse(a), MultiVector::Dense(b)) => {
                MultiVector::Dense(right_inner_dense(&a.to_dense(), b))
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{AlgebraConfig, MultiVector};
    use std::collections::HashMap;
    use std::sync::Arc;

    fn test_config() -> crate::AlgebraConfigRef {
        Arc::new(AlgebraConfig::euclidean(3))
    }

    // ==================== Left Inner Product Tests ====================

    #[test]
    fn test_left_inner_vectors() {
        let config = test_config();
        let e1 = MultiVector::basis_vector(config.clone(), 0);
        let e2 = MultiVector::basis_vector(config.clone(), 1);

        let result = e1.left_inner(&e1);
        assert!((result.get_coefficient(0) - 1.0).abs() < 1e-10);

        let result = e1.left_inner(&e2);
        assert!(result.is_zero());
    }

    #[test]
    fn test_left_inner_vector_bivector() {
        let config = test_config();
        let e1 = MultiVector::basis_vector(config.clone(), 0);
        let e2 = MultiVector::basis_vector(config.clone(), 1);
        let e12 = e1.geometric_product(&e2);

        let result = e1.left_inner(&e12);
        assert!((result.get_coefficient(2) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_left_inner_scalar() {
        let config = test_config();
        let scalar = MultiVector::from_scalar(config.clone(), 3.0);
        let e1 = MultiVector::basis_vector(config.clone(), 0);

        let result = scalar.left_inner(&e1);
        assert!((result.get_coefficient(1) - 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_left_inner_bivector_vector() {
        let config = test_config();
        let e1 = MultiVector::basis_vector(config.clone(), 0);
        let e2 = MultiVector::basis_vector(config.clone(), 1);
        let e12 = e1.geometric_product(&e2);

        let result = e12.left_inner(&e1);
        assert!(result.is_zero());
    }

    #[test]
    fn test_left_inner_trivector() {
        let config = test_config();
        let e1 = MultiVector::basis_vector(config.clone(), 0);
        let e2 = MultiVector::basis_vector(config.clone(), 1);
        let e3 = MultiVector::basis_vector(config.clone(), 2);
        let e12 = e1.geometric_product(&e2);
        let e123 = e12.geometric_product(&e3);

        let result = e1.left_inner(&e123);
        assert!((result.get_coefficient(6) - 1.0).abs() < 1e-10);

        let result = e12.left_inner(&e123);
        assert!((result.get_coefficient(4) + 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_left_inner_zero() {
        let config = test_config();
        let zero = MultiVector::zeros(config.clone());
        let e1 = MultiVector::basis_vector(config.clone(), 0);

        let result = zero.left_inner(&e1);
        assert!(result.is_zero());

        let result = e1.left_inner(&zero);
        assert!(result.is_zero());
    }

    // ==================== Right Inner Product Tests ====================

    #[test]
    fn test_right_inner() {
        let config = test_config();
        let e1 = MultiVector::basis_vector(config.clone(), 0);
        let _e2 = MultiVector::basis_vector(config.clone(), 1);

        let result = e1.right_inner(&e1);
        assert!((result.get_coefficient(0) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_right_inner_bivector_vector() {
        let config = test_config();
        let e1 = MultiVector::basis_vector(config.clone(), 0);
        let e2 = MultiVector::basis_vector(config.clone(), 1);
        let e12 = e1.geometric_product(&e2);

        let result = e12.right_inner(&e2);
        assert!((result.get_coefficient(1) + 1.0).abs() < 1e-10);

        let result = e12.right_inner(&e1);
        assert!((result.get_coefficient(2) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_right_inner_vector_bivector() {
        let config = test_config();
        let e1 = MultiVector::basis_vector(config.clone(), 0);
        let e2 = MultiVector::basis_vector(config.clone(), 1);
        let e12 = e1.geometric_product(&e2);

        let result = e1.right_inner(&e12);
        assert!(result.is_zero());
    }

    #[test]
    fn test_right_inner_scalar() {
        let config = test_config();
        let scalar = MultiVector::from_scalar(config.clone(), 5.0);
        let e1 = MultiVector::basis_vector(config.clone(), 0);

        let result = e1.right_inner(&scalar);
        assert!((result.get_coefficient(1) - 5.0).abs() < 1e-10);

        let result = scalar.right_inner(&scalar);
        assert!((result.get_coefficient(0) - 25.0).abs() < 1e-10);
    }

    #[test]
    fn test_right_inner_trivector() {
        let config = test_config();
        let e1 = MultiVector::basis_vector(config.clone(), 0);
        let e2 = MultiVector::basis_vector(config.clone(), 1);
        let e3 = MultiVector::basis_vector(config.clone(), 2);
        let e12 = e1.geometric_product(&e2);
        let e123 = e12.geometric_product(&e3);

        let result = e123.right_inner(&e1);
        assert!((result.get_coefficient(6) - 1.0).abs() < 1e-10);

        let result = e123.right_inner(&e3);
        assert!((result.get_coefficient(3) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_right_inner_zero() {
        let config = test_config();
        let zero = MultiVector::zeros(config.clone());
        let e1 = MultiVector::basis_vector(config.clone(), 0);

        let result = zero.right_inner(&e1);
        assert!(result.is_zero());

        let result = e1.right_inner(&zero);
        assert!(result.is_zero());
    }

    // ==================== Sparse Multivector Tests ====================

    #[test]
    fn test_left_inner_sparse_vectors() {
        let config = test_config();
        let e1_sparse = MultiVector::Sparse(SparseMultiVector::basis_vector(config.clone(), 0));
        let e2_sparse = MultiVector::Sparse(SparseMultiVector::basis_vector(config.clone(), 1));

        let result = e1_sparse.left_inner(&e1_sparse);
        assert!((result.get_coefficient(0) - 1.0).abs() < 1e-10);

        let result = e1_sparse.left_inner(&e2_sparse);
        assert!(result.is_zero());
    }

    #[test]
    fn test_left_inner_sparse_bivector() {
        let config = test_config();
        let e1_sparse = MultiVector::Sparse(SparseMultiVector::basis_vector(config.clone(), 0));

        let mut coeffs: HashMap<BasisIndex, f64> = HashMap::new();
        coeffs.insert(3, 1.0);
        let e12_sparse = MultiVector::Sparse(SparseMultiVector::from_map(config.clone(), coeffs));

        let result = e1_sparse.left_inner(&e12_sparse);
        assert!((result.get_coefficient(2) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_right_inner_sparse_vectors() {
        let config = test_config();
        let e1_sparse = MultiVector::Sparse(SparseMultiVector::basis_vector(config.clone(), 0));
        let e2_sparse = MultiVector::Sparse(SparseMultiVector::basis_vector(config.clone(), 1));

        let result = e1_sparse.right_inner(&e1_sparse);
        assert!((result.get_coefficient(0) - 1.0).abs() < 1e-10);

        let result = e1_sparse.right_inner(&e2_sparse);
        assert!(result.is_zero());
    }

    #[test]
    fn test_right_inner_sparse_bivector() {
        let config = test_config();
        let e2_sparse = MultiVector::Sparse(SparseMultiVector::basis_vector(config.clone(), 1));

        let mut coeffs: HashMap<BasisIndex, f64> = HashMap::new();
        coeffs.insert(3, 1.0);
        let e12_sparse = MultiVector::Sparse(SparseMultiVector::from_map(config.clone(), coeffs));

        let result = e12_sparse.right_inner(&e2_sparse);
        assert!((result.get_coefficient(1) + 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_left_inner_sparse_direct() {
        let config = test_config();
        let mut coeffs_a: HashMap<BasisIndex, f64> = HashMap::new();
        coeffs_a.insert(1, 2.0);

        let mut coeffs_b: HashMap<BasisIndex, f64> = HashMap::new();
        coeffs_b.insert(1, 3.0);

        let a = SparseMultiVector::from_map(config.clone(), coeffs_a);
        let b = SparseMultiVector::from_map(config.clone(), coeffs_b);

        let result = left_inner_sparse(&a, &b);
        assert!((result.scalar_part() - 6.0).abs() < 1e-10);
    }

    #[test]
    fn test_right_inner_sparse_direct() {
        let config = test_config();
        let mut coeffs_a: HashMap<BasisIndex, f64> = HashMap::new();
        coeffs_a.insert(1, 2.0);

        let mut coeffs_b: HashMap<BasisIndex, f64> = HashMap::new();
        coeffs_b.insert(1, 3.0);

        let a = SparseMultiVector::from_map(config.clone(), coeffs_a);
        let b = SparseMultiVector::from_map(config.clone(), coeffs_b);

        let result = right_inner_sparse(&a, &b);
        assert!((result.scalar_part() - 6.0).abs() < 1e-10);
    }

    // ==================== Mixed Dense/Sparse Tests ====================

    #[test]
    fn test_left_inner_dense_sparse() {
        let config = test_config();
        let e1 = MultiVector::basis_vector(config.clone(), 0);
        let e2_sparse = MultiVector::Sparse(SparseMultiVector::basis_vector(config.clone(), 1));

        let result = e1.left_inner(&e2_sparse);
        assert!(result.is_zero());

        let e1_sparse = MultiVector::Sparse(SparseMultiVector::basis_vector(config.clone(), 0));
        let result = e1.left_inner(&e1_sparse);
        assert!((result.get_coefficient(0) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_left_inner_sparse_dense() {
        let config = test_config();
        let e1_sparse = MultiVector::Sparse(SparseMultiVector::basis_vector(config.clone(), 0));
        let e2 = MultiVector::basis_vector(config.clone(), 1);

        let result = e1_sparse.left_inner(&e2);
        assert!(result.is_zero());

        let e1 = MultiVector::basis_vector(config.clone(), 0);
        let result = e1_sparse.left_inner(&e1);
        assert!((result.get_coefficient(0) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_right_inner_dense_sparse() {
        let config = test_config();
        let e1 = MultiVector::basis_vector(config.clone(), 0);
        let e2_sparse = MultiVector::Sparse(SparseMultiVector::basis_vector(config.clone(), 1));

        let result = e1.right_inner(&e2_sparse);
        assert!(result.is_zero());

        let e1_sparse = MultiVector::Sparse(SparseMultiVector::basis_vector(config.clone(), 0));
        let result = e1.right_inner(&e1_sparse);
        assert!((result.get_coefficient(0) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_right_inner_sparse_dense() {
        let config = test_config();
        let e1_sparse = MultiVector::Sparse(SparseMultiVector::basis_vector(config.clone(), 0));
        let e2 = MultiVector::basis_vector(config.clone(), 1);

        let result = e1_sparse.right_inner(&e2);
        assert!(result.is_zero());

        let e1 = MultiVector::basis_vector(config.clone(), 0);
        let result = e1_sparse.right_inner(&e1);
        assert!((result.get_coefficient(0) - 1.0).abs() < 1e-10);
    }

    // ==================== Parallel Path Tests (6D+) ====================

    #[test]
    fn test_left_inner_parallel() {
        let config = Arc::new(AlgebraConfig::euclidean(6));
        let e1 = MultiVector::basis_vector(config.clone(), 0);
        let e2 = MultiVector::basis_vector(config.clone(), 1);

        let result = e1.left_inner(&e1);
        assert!((result.get_coefficient(0) - 1.0).abs() < 1e-10);

        let result = e1.left_inner(&e2);
        assert!(result.is_zero());

        let e12 = e1.geometric_product(&e2);
        let result = e1.left_inner(&e12);
        assert!((result.get_coefficient(2) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_right_inner_parallel() {
        let config = Arc::new(AlgebraConfig::euclidean(6));
        let e1 = MultiVector::basis_vector(config.clone(), 0);
        let e2 = MultiVector::basis_vector(config.clone(), 1);

        let result = e1.right_inner(&e1);
        assert!((result.get_coefficient(0) - 1.0).abs() < 1e-10);

        let result = e1.right_inner(&e2);
        assert!(result.is_zero());

        let e12 = e1.geometric_product(&e2);
        let result = e12.right_inner(&e2);
        assert!((result.get_coefficient(1) + 1.0).abs() < 1e-10);
    }

    // ==================== Edge Cases ====================

    #[test]
    fn test_left_inner_orthogonal() {
        let config = test_config();
        let e1 = MultiVector::basis_vector(config.clone(), 0);
        let e2 = MultiVector::basis_vector(config.clone(), 1);
        let e3 = MultiVector::basis_vector(config.clone(), 2);

        assert!(e1.left_inner(&e2).is_zero());
        assert!(e1.left_inner(&e3).is_zero());
        assert!(e2.left_inner(&e1).is_zero());
        assert!(e2.left_inner(&e3).is_zero());
        assert!(e3.left_inner(&e1).is_zero());
        assert!(e3.left_inner(&e2).is_zero());
    }

    #[test]
    fn test_right_inner_orthogonal() {
        let config = test_config();
        let e1 = MultiVector::basis_vector(config.clone(), 0);
        let e2 = MultiVector::basis_vector(config.clone(), 1);
        let e3 = MultiVector::basis_vector(config.clone(), 2);

        assert!(e1.right_inner(&e2).is_zero());
        assert!(e1.right_inner(&e3).is_zero());
        assert!(e2.right_inner(&e1).is_zero());
        assert!(e2.right_inner(&e3).is_zero());
        assert!(e3.right_inner(&e1).is_zero());
        assert!(e3.right_inner(&e2).is_zero());
    }

    #[test]
    fn test_left_inner_scalar_scalar() {
        let config = test_config();
        let a = MultiVector::from_scalar(config.clone(), 3.0);
        let b = MultiVector::from_scalar(config.clone(), 4.0);

        let result = a.left_inner(&b);
        assert!((result.get_coefficient(0) - 12.0).abs() < 1e-10);
    }

    #[test]
    fn test_inner_with_coefficients() {
        let config = test_config();
        let mut coeffs = vec![0.0; 8];
        coeffs[1] = 2.0;
        let a = MultiVector::from_coefficients(config.clone(), coeffs);

        let mut coeffs = vec![0.0; 8];
        coeffs[1] = 3.0;
        let b = MultiVector::from_coefficients(config.clone(), coeffs);

        let result = a.left_inner(&b);
        assert!((result.get_coefficient(0) - 6.0).abs() < 1e-10);

        let result = a.right_inner(&b);
        assert!((result.get_coefficient(0) - 6.0).abs() < 1e-10);
    }

    #[test]
    fn test_left_inner_different_grades() {
        let config = test_config();
        let e1 = MultiVector::basis_vector(config.clone(), 0);

        let mut coeffs = vec![0.0; 8];
        coeffs[0] = 1.0;
        coeffs[1] = 2.0;
        coeffs[3] = 3.0;
        let mv = MultiVector::from_coefficients(config.clone(), coeffs);

        let result = e1.left_inner(&mv);
        assert!((result.get_coefficient(0) - 2.0).abs() < 1e-10);
        assert!((result.get_coefficient(2) - 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_right_inner_different_grades() {
        let config = test_config();
        let e2 = MultiVector::basis_vector(config.clone(), 1);

        let mut coeffs = vec![0.0; 8];
        coeffs[0] = 1.0;
        coeffs[1] = 2.0;
        coeffs[3] = 3.0;
        let mv = MultiVector::from_coefficients(config.clone(), coeffs);

        let result = mv.right_inner(&e2);
        assert!((result.get_coefficient(1) + 3.0).abs() < 1e-10);
    }

    // ==================== Direct Dense Function Tests ====================

    #[test]
    fn test_left_inner_dense_direct() {
        let config = test_config();
        let a = DenseMultiVector::from_coefficients(
            config.clone(),
            vec![0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        );
        let b = DenseMultiVector::from_coefficients(
            config.clone(),
            vec![0.0, 3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        );

        let result = left_inner_dense(&a, &b);
        assert!((result.get_coefficient(0) - 6.0).abs() < 1e-10);
    }

    #[test]
    fn test_right_inner_dense_direct() {
        let config = test_config();
        let a = DenseMultiVector::from_coefficients(
            config.clone(),
            vec![0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        );
        let b = DenseMultiVector::from_coefficients(
            config.clone(),
            vec![0.0, 3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        );

        let result = right_inner_dense(&a, &b);
        assert!((result.get_coefficient(0) - 6.0).abs() < 1e-10);
    }

    #[test]
    fn test_left_inner_dense_all_nonzero() {
        let config = test_config();
        let a = DenseMultiVector::from_coefficients(config.clone(), vec![1.0; 8]);
        let b = DenseMultiVector::from_coefficients(config.clone(), vec![1.0; 8]);
        let result = left_inner_dense(&a, &b);
        assert!(!result.coefficients.is_empty());
    }

    #[test]
    fn test_right_inner_dense_all_nonzero() {
        let config = test_config();
        let a = DenseMultiVector::from_coefficients(config.clone(), vec![1.0; 8]);
        let b = DenseMultiVector::from_coefficients(config.clone(), vec![1.0; 8]);
        let result = right_inner_dense(&a, &b);
        assert!(!result.coefficients.is_empty());
    }
}
