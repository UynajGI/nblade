//! 几何积模块
//!
//! 实现多向量的几何积（Geometric Product）
//! 公式：AB = Σ⟨ArBs⟩_{|r-s|+2j}

use crate::algebra_config::AlgebraConfigRef;
use crate::multiplication_table::GeometricProductTable;
use crate::multivector::{DenseMultiVector, MultiVector, SparseMultiVector};
use crate::products::should_use_parallel;
use rayon::prelude::*;

#[cfg(feature = "pool")]
use crate::multivector::pool::BufferPool;

/// 密集多向量的几何积（使用预计算乘法表）
/// 支持自适应并行调度和 SIMD 加速
pub fn geometric_product_dense(
    a: &DenseMultiVector,
    b: &DenseMultiVector,
    table: &GeometricProductTable,
) -> DenseMultiVector {
    let config = a.config.clone();
    let n = config.dimension() as usize;

    #[cfg(feature = "simd")]
    if n <= 4 {
        return crate::products::geometric_simd::geometric_product_dense_optimized(a, b, table);
    }

    let size = 1 << n;

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
        geometric_product_dense_parallel(&non_zero_a, &non_zero_b, table, config, size)
    } else {
        geometric_product_dense_sequential(&non_zero_a, &non_zero_b, table, config, size)
    }
}

fn geometric_product_dense_sequential(
    non_zero_a: &[(usize, &f64)],
    non_zero_b: &[(usize, &f64)],
    table: &GeometricProductTable,
    config: AlgebraConfigRef,
    size: usize,
) -> DenseMultiVector {
    #[cfg(feature = "pool")]
    let mut result = BufferPool::get_buffer(size);

    #[cfg(not(feature = "pool"))]
    let mut result = vec![0.0f64; size];

    for (i, &coef_a) in non_zero_a {
        for (j, &coef_b) in non_zero_b {
            let entry = table.get(*i, *j);
            if !entry.is_zero() {
                let result_idx = entry.result_index as usize;
                result[result_idx] += entry.sign as f64 * entry.metric_factor * coef_a * coef_b;
            }
        }
    }

    #[cfg(feature = "pool")]
    {
        DenseMultiVector::from_pooled_buffer(config, result)
    }

    #[cfg(not(feature = "pool"))]
    {
        DenseMultiVector::from_coefficients(config, result)
    }
}

fn geometric_product_dense_parallel(
    non_zero_a: &[(usize, &f64)],
    non_zero_b: &[(usize, &f64)],
    table: &GeometricProductTable,
    config: AlgebraConfigRef,
    size: usize,
) -> DenseMultiVector {
    #[cfg(feature = "pool")]
    {
        let result: Vec<f64> = non_zero_a
            .par_iter()
            .fold(
                || BufferPool::get_buffer(size),
                |mut local_result, (i, &coef_a)| {
                    for (j, &coef_b) in non_zero_b {
                        let entry = table.get(*i, *j);
                        if !entry.is_zero() {
                            let result_idx = entry.result_index as usize;
                            local_result[result_idx] +=
                                entry.sign as f64 * entry.metric_factor * coef_a * coef_b;
                        }
                    }
                    local_result
                },
            )
            .reduce(
                || BufferPool::get_buffer(size),
                |mut a, b| {
                    for i in 0..size {
                        a[i] += b[i];
                    }
                    BufferPool::return_buffer(b);
                    a
                },
            );

        DenseMultiVector::from_pooled_buffer(config, result)
    }

    #[cfg(not(feature = "pool"))]
    {
        let result: Vec<f64> = non_zero_a
            .par_iter()
            .fold(
                || vec![0.0f64; size],
                |mut local_result, (i, &coef_a)| {
                    for (j, &coef_b) in non_zero_b {
                        let entry = table.get(*i, *j);
                        if !entry.is_zero() {
                            let result_idx = entry.result_index as usize;
                            local_result[result_idx] +=
                                entry.sign as f64 * entry.metric_factor * coef_a * coef_b;
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
}

/// 稀疏多向量的几何积
pub fn geometric_product_sparse(
    a: &SparseMultiVector,
    b: &SparseMultiVector,
    config: &AlgebraConfigRef,
) -> SparseMultiVector {
    use crate::basis::index::basis_geometric_product;

    let mut result_coeffs = std::collections::HashMap::new();

    for (&i, &coef_a) in &a.coefficients {
        for (&j, &coef_b) in &b.coefficients {
            let (result, sign, metric) = basis_geometric_product(i, j, &config.signature);

            if metric != 0.0 {
                let product = sign as f64 * metric * coef_a * coef_b;
                *result_coeffs.entry(result).or_insert(0.0) += product;
            }
        }
    }

    // 移除接近零的系数
    result_coeffs.retain(|_, &mut v| v.abs() > 1e-15);

    SparseMultiVector::from_map(config.clone(), result_coeffs)
}

impl MultiVector {
    /// 几何积
    pub fn geometric_product(&self, other: &Self) -> Self {
        match (self, other) {
            (MultiVector::Dense(a), MultiVector::Dense(b)) => {
                let table = GeometricProductTable::new(a.config());
                let result = geometric_product_dense(a, b, &table);
                MultiVector::Dense(result)
            }
            (MultiVector::Sparse(a), MultiVector::Sparse(b)) => {
                let result = geometric_product_sparse(a, b, self.config());
                MultiVector::Sparse(result)
            }
            (MultiVector::Dense(a), MultiVector::Sparse(b)) => {
                let b_dense = b.to_dense();
                let table = GeometricProductTable::new(a.config());
                let result = geometric_product_dense(a, &b_dense, &table);
                MultiVector::Dense(result)
            }
            (MultiVector::Sparse(a), MultiVector::Dense(b)) => {
                let a_dense = a.to_dense();
                let table = GeometricProductTable::new(a.config());
                let result = geometric_product_dense(&a_dense, b, &table);
                MultiVector::Dense(result)
            }
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

    #[test]
    fn test_geometric_product_vectors() {
        let config = test_config();
        let e1 = MultiVector::basis_vector(config.clone(), 0);
        let e2 = MultiVector::basis_vector(config.clone(), 1);

        // e1 * e1 = 1
        let result = e1.geometric_product(&e1);
        assert!((result.get_coefficient(0) - 1.0).abs() < 1e-10);

        // e1 * e2 = e12
        let result = e1.geometric_product(&e2);
        assert!((result.get_coefficient(3) - 1.0).abs() < 1e-10);

        // e2 * e1 = -e12
        let result = e2.geometric_product(&e1);
        assert!((result.get_coefficient(3) + 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_geometric_product_scalar() {
        let config = test_config();
        let scalar = MultiVector::from_scalar(config.clone(), 2.0);
        let e1 = MultiVector::basis_vector(config.clone(), 0);

        let result = scalar.geometric_product(&e1);
        assert!((result.get_coefficient(1) - 2.0).abs() < 1e-10);
    }
}
