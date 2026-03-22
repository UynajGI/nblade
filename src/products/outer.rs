//! 外积模块
//!
//! 实现多向量的外积（Wedge Product / Outer Product）
//!
//! 公式：A∧B = Σ⟨ArBs⟩_{r+s}

use crate::basis::index::basis_outer_product;
use crate::multivector::{DenseMultiVector, MultiVector, SparseMultiVector};
use crate::products::should_use_parallel;
use rayon::prelude::*;

/// 密集多向量的外积（自适应并行调度）
pub fn outer_product_dense(a: &DenseMultiVector, b: &DenseMultiVector) -> DenseMultiVector {
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
        outer_product_dense_parallel(&non_zero_a, &non_zero_b, config, size)
    } else {
        outer_product_dense_sequential(&non_zero_a, &non_zero_b, config, size)
    }
}

fn outer_product_dense_sequential(
    non_zero_a: &[(usize, &f64)],
    non_zero_b: &[(usize, &f64)],
    config: crate::algebra_config::AlgebraConfigRef,
    size: usize,
) -> DenseMultiVector {
    let mut result = vec![0.0f64; size];

    for (i, &coef_a) in non_zero_a {
        for (j, &coef_b) in non_zero_b {
            if let Some((result_idx, sign)) = basis_outer_product(*i as u64, *j as u64) {
                result[result_idx as usize] += sign as f64 * coef_a * coef_b;
            }
        }
    }

    DenseMultiVector::from_coefficients(config, result)
}

fn outer_product_dense_parallel(
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
                    if let Some((result_idx, sign)) = basis_outer_product(*i as u64, *j as u64) {
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

/// 稀疏多向量的外积
pub fn outer_product_sparse(a: &SparseMultiVector, b: &SparseMultiVector) -> SparseMultiVector {
    let mut result_coeffs = std::collections::HashMap::new();

    for (&i, &coef_a) in &a.coefficients {
        for (&j, &coef_b) in &b.coefficients {
            if let Some((result, sign)) = basis_outer_product(i, j) {
                let product = sign as f64 * coef_a * coef_b;
                *result_coeffs.entry(result).or_insert(0.0) += product;
            }
        }
    }

    result_coeffs.retain(|_, &mut v| v.abs() > 1e-15);

    SparseMultiVector::from_map(a.config.clone(), result_coeffs)
}

impl MultiVector {
    /// 外积 A∧B
    pub fn outer_product(&self, other: &Self) -> Self {
        match (self, other) {
            (MultiVector::Dense(a), MultiVector::Dense(b)) => {
                MultiVector::Dense(outer_product_dense(a, b))
            }
            (MultiVector::Sparse(a), MultiVector::Sparse(b)) => {
                MultiVector::Sparse(outer_product_sparse(a, b))
            }
            (MultiVector::Dense(a), MultiVector::Sparse(b)) => {
                MultiVector::Dense(outer_product_dense(a, &b.to_dense()))
            }
            (MultiVector::Sparse(a), MultiVector::Dense(b)) => {
                MultiVector::Dense(outer_product_dense(&a.to_dense(), b))
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{AlgebraConfig, MultiVector};
    use std::sync::Arc;

    fn test_config() -> crate::AlgebraConfigRef {
        Arc::new(AlgebraConfig::euclidean(3))
    }

    #[test]
    fn test_outer_product_vectors() {
        let config = test_config();
        let e1 = MultiVector::basis_vector(config.clone(), 0);
        let e2 = MultiVector::basis_vector(config.clone(), 1);

        // e1 ∧ e2 = e12
        let result = e1.outer_product(&e2);
        assert!((result.get_coefficient(3) - 1.0).abs() < 1e-10);

        // e2 ∧ e1 = -e12
        let result = e2.outer_product(&e1);
        assert!((result.get_coefficient(3) + 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_outer_product_same_vector() {
        let config = test_config();
        let e1 = MultiVector::basis_vector(config.clone(), 0);

        // e1 ∧ e1 = 0
        let result = e1.outer_product(&e1);
        assert!(result.is_zero());
    }

    #[test]
    fn test_outer_product_trivector() {
        let config = test_config();
        let e1 = MultiVector::basis_vector(config.clone(), 0);
        let e2 = MultiVector::basis_vector(config.clone(), 1);
        let e3 = MultiVector::basis_vector(config.clone(), 2);

        // e1 ∧ e2 ∧ e3 = e123
        let e12 = e1.outer_product(&e2);
        let e123 = e12.outer_product(&e3);
        assert!((e123.get_coefficient(7) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_outer_product_scalar() {
        let config = test_config();
        let scalar = MultiVector::from_scalar(config.clone(), 2.0);
        let e1 = MultiVector::basis_vector(config.clone(), 0);

        // 2 ∧ e1 = 2e1
        let result = scalar.outer_product(&e1);
        assert!((result.get_coefficient(1) - 2.0).abs() < 1e-10);
    }

    // ============== Sparse outer product tests ==============

    #[test]
    fn test_outer_product_sparse_vectors() {
        use crate::multivector::SparseMultiVector;

        let config = test_config();
        let e1 = SparseMultiVector::basis_vector(config.clone(), 0);
        let e2 = SparseMultiVector::basis_vector(config.clone(), 1);

        // e1 ∧ e2 = e12 (index 0b011 = 3)
        let result = outer_product_sparse(&e1, &e2);
        assert!((result.get_coefficient(3) - 1.0).abs() < 1e-10);

        // e2 ∧ e1 = -e12
        let result = outer_product_sparse(&e2, &e1);
        assert!((result.get_coefficient(3) + 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_outer_product_sparse_same_vector() {
        use crate::multivector::SparseMultiVector;

        let config = test_config();
        let e1 = SparseMultiVector::basis_vector(config.clone(), 0);

        // e1 ∧ e1 = 0
        let result = outer_product_sparse(&e1, &e1);
        assert!(result.coefficients.is_empty());
    }

    #[test]
    fn test_outer_product_sparse_trivector() {
        use crate::multivector::SparseMultiVector;

        let config = test_config();
        let e1 = SparseMultiVector::basis_vector(config.clone(), 0);
        let e2 = SparseMultiVector::basis_vector(config.clone(), 1);
        let e3 = SparseMultiVector::basis_vector(config.clone(), 2);

        // e1 ∧ e2 ∧ e3 = e123 (index 0b111 = 7)
        let e12 = outer_product_sparse(&e1, &e2);
        let e123 = outer_product_sparse(&e12, &e3);
        assert!((e123.get_coefficient(7) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_outer_product_sparse_scalar() {
        use crate::multivector::SparseMultiVector;

        let config = test_config();
        let scalar = SparseMultiVector::from_scalar(config.clone(), 3.0);
        let e1 = SparseMultiVector::basis_vector(config.clone(), 0);

        // 3 ∧ e1 = 3e1
        let result = outer_product_sparse(&scalar, &e1);
        assert!((result.get_coefficient(1) - 3.0).abs() < 1e-10);

        // e1 ∧ 3 = 3e1
        let result = outer_product_sparse(&e1, &scalar);
        assert!((result.get_coefficient(1) - 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_outer_product_sparse_empty() {
        use crate::multivector::SparseMultiVector;

        let config = test_config();
        let empty = SparseMultiVector::zeros(config.clone());
        let e1 = SparseMultiVector::basis_vector(config.clone(), 0);

        // 0 ∧ e1 = 0
        let result = outer_product_sparse(&empty, &e1);
        assert!(result.coefficients.is_empty());

        // e1 ∧ 0 = 0
        let result = outer_product_sparse(&e1, &empty);
        assert!(result.coefficients.is_empty());
    }

    #[test]
    fn test_outer_product_sparse_bivector_vector() {
        use crate::multivector::SparseMultiVector;

        let config = test_config();
        let e1 = SparseMultiVector::basis_vector(config.clone(), 0);
        let e2 = SparseMultiVector::basis_vector(config.clone(), 1);
        let e3 = SparseMultiVector::basis_vector(config.clone(), 2);

        // e12 ∧ e3 = e123
        let e12 = outer_product_sparse(&e1, &e2);
        let result = outer_product_sparse(&e12, &e3);
        assert!((result.get_coefficient(7) - 1.0).abs() < 1e-10);

        // e12 ∧ e1 = 0 (linearly dependent)
        let result = outer_product_sparse(&e12, &e1);
        assert!(result.coefficients.is_empty());
    }

    // ============== Mixed Dense/Sparse outer product tests ==============

    #[test]
    fn test_outer_product_dense_sparse() {
        let config = test_config();
        let e1 = MultiVector::basis_vector(config.clone(), 0);
        let e2 = MultiVector::basis_vector(config.clone(), 1);

        // Convert e2 to sparse
        let sparse_e2 = MultiVector::Sparse(crate::multivector::SparseMultiVector::basis_vector(
            config.clone(),
            1,
        ));

        // Dense ∧ Sparse
        let result = e1.outer_product(&sparse_e2);
        assert!((result.get_coefficient(3) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_outer_product_sparse_dense() {
        let config = test_config();
        let e1 = MultiVector::basis_vector(config.clone(), 0);

        // Convert to sparse e1
        let sparse_e1 = MultiVector::Sparse(crate::multivector::SparseMultiVector::basis_vector(
            config.clone(),
            0,
        ));
        let e2 = MultiVector::basis_vector(config.clone(), 1);

        // Sparse ∧ Dense
        let result = sparse_e1.outer_product(&e2);
        assert!((result.get_coefficient(3) - 1.0).abs() < 1e-10);
    }

    // ============== Parallel path test (dimension >= 6) ==============

    #[test]
    fn test_outer_product_parallel_path() {
        // Dimension 6 has 64 basis elements, triggering parallel path
        let config = Arc::new(AlgebraConfig::euclidean(6));
        let e1 = MultiVector::basis_vector(config.clone(), 0);
        let e2 = MultiVector::basis_vector(config.clone(), 1);
        let e3 = MultiVector::basis_vector(config.clone(), 2);
        let e4 = MultiVector::basis_vector(config.clone(), 3);
        let e5 = MultiVector::basis_vector(config.clone(), 4);
        let e6 = MultiVector::basis_vector(config.clone(), 5);

        // e1 ∧ e2 = e12 (index 0b000011 = 3)
        let result = e1.outer_product(&e2);
        assert!((result.get_coefficient(3) - 1.0).abs() < 1e-10);

        // e1 ∧ e6 = e16 (index 0b100001 = 33)
        let result = e1.outer_product(&e6);
        assert!((result.get_coefficient(33) - 1.0).abs() < 1e-10);

        // Build a multivector with many coefficients to trigger parallel
        let mut coeffs = vec![0.0; config.basis_count()];
        for i in 0..20 {
            coeffs[i] = (i + 1) as f64 * 0.1;
        }
        let dense_a = MultiVector::from_coefficients(config.clone(), coeffs);
        let dense_b = MultiVector::basis_vector(config.clone(), 5);
        let _result = dense_a.outer_product(&dense_b);
        // Just verify it doesn't crash and produces some result
    }

    #[test]
    fn test_outer_product_high_dimension() {
        // Test with 7 dimensions to ensure parallel path
        let config = Arc::new(AlgebraConfig::euclidean(7));
        let e1 = MultiVector::basis_vector(config.clone(), 0);
        let e7 = MultiVector::basis_vector(config.clone(), 6);

        // e1 ∧ e7 = e17 (index 0b1000001 = 65)
        let result = e1.outer_product(&e7);
        assert!((result.get_coefficient(65) - 1.0).abs() < 1e-10);

        // e1 ∧ e1 = 0
        let result = e1.outer_product(&e1);
        assert!(result.is_zero());
    }

    // ============== Edge cases ==============

    #[test]
    fn test_outer_product_zero_multivector() {
        let config = test_config();
        let zero = MultiVector::zeros(config.clone());
        let e1 = MultiVector::basis_vector(config.clone(), 0);

        // 0 ∧ e1 = 0
        let result = zero.outer_product(&e1);
        assert!(result.is_zero());

        // e1 ∧ 0 = 0
        let result = e1.outer_product(&zero);
        assert!(result.is_zero());
    }

    #[test]
    fn test_outer_product_bivector_zero_result() {
        let config = test_config();
        let e1 = MultiVector::basis_vector(config.clone(), 0);
        let e2 = MultiVector::basis_vector(config.clone(), 1);

        // e1 ∧ e2 = e12
        let e12 = e1.outer_product(&e2);

        // e12 ∧ e1 = 0 (linearly dependent)
        let result = e12.outer_product(&e1);
        assert!(result.is_zero());
    }

    #[test]
    fn test_outer_product_full_grade() {
        let config = test_config();
        let e1 = MultiVector::basis_vector(config.clone(), 0);
        let e2 = MultiVector::basis_vector(config.clone(), 1);
        let e3 = MultiVector::basis_vector(config.clone(), 2);

        // e123 (grade 3) is the highest grade in 3D
        let e123 = e1.outer_product(&e2).outer_product(&e3);

        // e123 ∧ anything = 0 in 3D (grade overflow)
        let result = e123.outer_product(&e1);
        assert!(result.is_zero());
    }

    #[test]
    fn test_outer_product_coefficient_accumulation() {
        let config = test_config();
        // Create multivector with multiple terms: 2*e1 + 3*e2
        let e1 = MultiVector::basis_vector(config.clone(), 0);
        let e2 = MultiVector::basis_vector(config.clone(), 1);
        let e3 = MultiVector::basis_vector(config.clone(), 2);

        let a = e1.scale(2.0).add(&e2.scale(3.0));
        let b = e3.scale(4.0);

        // (2e1 + 3e2) ∧ (4e3) = 8*e13 + 12*e23
        let result = a.outer_product(&b);
        assert!((result.get_coefficient(5) - 8.0).abs() < 1e-10); // e13 index = 0b101 = 5
        assert!((result.get_coefficient(6) - 12.0).abs() < 1e-10); // e23 index = 0b110 = 6
    }
}
