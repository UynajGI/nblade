//! SIMD 加速的几何积模块 / SIMD-accelerated geometric product module
//!
//! 针对 2D、3D、4D 几何代数的 SIMD 优化实现
//! SIMD optimized implementation for 2D, 3D, 4D geometric algebra

#[cfg(feature = "simd")]
use wide::f64x4;

use crate::multiplication_table::GeometricProductTable;
use crate::multivector::DenseMultiVector;

/// 使用 SIMD 计算几何积的核心循环（2D 情况）
/// Core SIMD loop for geometric product (2D case)
#[cfg(feature = "simd")]
#[inline(always)]
fn geometric_product_simd_core_2d(
    a: &[f64; 4],
    b: &[f64; 4],
    table: &GeometricProductTable,
) -> [f64; 4] {
    let mut result = [0.0; 4];

    for i in 0..4 {
        let a_i = a[i];
        if a_i.abs() > 1e-15 {
            let a_splat = f64x4::splat(a_i);
            let b_vec = f64x4::from(*b);
            let scaled = a_splat * b_vec;
            let scaled_arr = scaled.as_array_ref();

            for j in 0..4 {
                let entry = table.get(i, j);
                if !entry.is_zero() {
                    result[entry.result_index as usize] +=
                        entry.sign as f64 * entry.metric_factor * scaled_arr[j];
                }
            }
        }
    }

    result
}

/// 使用 SIMD 计算几何积的核心循环（3D 情况）
/// Core SIMD loop for geometric product (3D case)
#[cfg(feature = "simd")]
#[inline(always)]
fn geometric_product_simd_core_3d(
    a: &[f64; 8],
    b: &[f64; 8],
    table: &GeometricProductTable,
) -> [f64; 8] {
    let mut result = [0.0; 8];

    for i in 0..8 {
        let a_i = a[i];
        if a_i.abs() > 1e-15 {
            let a_splat = f64x4::splat(a_i);

            // Process b in two chunks of 4
            let b_low: [f64; 4] = [b[0], b[1], b[2], b[3]];
            let b_high: [f64; 4] = [b[4], b[5], b[6], b[7]];

            let b_low_vec = f64x4::from(b_low);
            let b_high_vec = f64x4::from(b_high);

            let scaled_low = a_splat * b_low_vec;
            let scaled_high = a_splat * b_high_vec;

            let scaled_low_arr = scaled_low.as_array_ref();
            let scaled_high_arr = scaled_high.as_array_ref();

            for j in 0..4 {
                let entry = table.get(i, j);
                if !entry.is_zero() {
                    result[entry.result_index as usize] +=
                        entry.sign as f64 * entry.metric_factor * scaled_low_arr[j];
                }
            }
            for j in 4..8 {
                let entry = table.get(i, j);
                if !entry.is_zero() {
                    result[entry.result_index as usize] +=
                        entry.sign as f64 * entry.metric_factor * scaled_high_arr[j - 4];
                }
            }
        }
    }

    result
}

/// 使用 SIMD 计算几何积的核心循环（4D 情况）
/// Core SIMD loop for geometric product (4D case)
#[cfg(feature = "simd")]
#[inline(always)]
fn geometric_product_simd_core_4d(
    a: &[f64; 16],
    b: &[f64; 16],
    table: &GeometricProductTable,
) -> [f64; 16] {
    let mut result = [0.0; 16];

    // Pre-split b into chunks
    let b_chunks: [[f64; 4]; 4] = [
        [b[0], b[1], b[2], b[3]],
        [b[4], b[5], b[6], b[7]],
        [b[8], b[9], b[10], b[11]],
        [b[12], b[13], b[14], b[15]],
    ];

    for i in 0..16 {
        let a_i = a[i];
        if a_i.abs() > 1e-15 {
            let a_splat = f64x4::splat(a_i);

            for chunk in 0..4 {
                let j_start = chunk * 4;
                let b_chunk = b_chunks[chunk];
                let b_vec = f64x4::from(b_chunk);
                let scaled = a_splat * b_vec;
                let scaled_arr = scaled.as_array_ref();

                for j in 0..4 {
                    let entry = table.get(i, j_start + j);
                    if !entry.is_zero() {
                        result[entry.result_index as usize] +=
                            entry.sign as f64 * entry.metric_factor * scaled_arr[j];
                    }
                }
            }
        }
    }

    result
}

/// 标量计算几何积（2D 情况）
/// Scalar geometric product (2D case)
#[inline(always)]
fn geometric_product_scalar_core_2d(
    a: &[f64; 4],
    b: &[f64; 4],
    table: &GeometricProductTable,
) -> [f64; 4] {
    let mut result = [0.0; 4];

    for i in 0..4 {
        if a[i].abs() > 1e-15 {
            for j in 0..4 {
                if b[j].abs() > 1e-15 {
                    let entry = table.get(i, j);
                    if !entry.is_zero() {
                        result[entry.result_index as usize] +=
                            entry.sign as f64 * entry.metric_factor * a[i] * b[j];
                    }
                }
            }
        }
    }

    result
}

/// 标量计算几何积（3D 情况）
/// Scalar geometric product (3D case)
#[inline(always)]
fn geometric_product_scalar_core_3d(
    a: &[f64; 8],
    b: &[f64; 8],
    table: &GeometricProductTable,
) -> [f64; 8] {
    let mut result = [0.0; 8];

    for i in 0..8 {
        if a[i].abs() > 1e-15 {
            for j in 0..8 {
                if b[j].abs() > 1e-15 {
                    let entry = table.get(i, j);
                    if !entry.is_zero() {
                        result[entry.result_index as usize] +=
                            entry.sign as f64 * entry.metric_factor * a[i] * b[j];
                    }
                }
            }
        }
    }

    result
}

/// 标量计算几何积（4D 情况）
/// Scalar geometric product (4D case)
#[inline(always)]
fn geometric_product_scalar_core_4d(
    a: &[f64; 16],
    b: &[f64; 16],
    table: &GeometricProductTable,
) -> [f64; 16] {
    let mut result = [0.0; 16];

    for i in 0..16 {
        if a[i].abs() > 1e-15 {
            for j in 0..16 {
                if b[j].abs() > 1e-15 {
                    let entry = table.get(i, j);
                    if !entry.is_zero() {
                        result[entry.result_index as usize] +=
                            entry.sign as f64 * entry.metric_factor * a[i] * b[j];
                    }
                }
            }
        }
    }

    result
}

/// 为 DenseMultiVector 提供的 SIMD 加速几何积（2D）
/// SIMD-accelerated geometric product for DenseMultiVector (2D)
#[cfg(feature = "simd")]
pub fn geometric_product_dense_2d_simd(
    a: &DenseMultiVector,
    b: &DenseMultiVector,
    table: &GeometricProductTable,
) -> DenseMultiVector {
    let config = a.config.clone();
    let a_arr: [f64; 4] = [
        a.coefficients[0],
        a.coefficients[1],
        a.coefficients[2],
        a.coefficients[3],
    ];
    let b_arr: [f64; 4] = [
        b.coefficients[0],
        b.coefficients[1],
        b.coefficients[2],
        b.coefficients[3],
    ];

    let result = geometric_product_simd_core_2d(&a_arr, &b_arr, table);
    DenseMultiVector::from_coefficients(config, result.to_vec())
}

/// 为 DenseMultiVector 提供的 SIMD 加速几何积（3D）
/// SIMD-accelerated geometric product for DenseMultiVector (3D)
#[cfg(feature = "simd")]
pub fn geometric_product_dense_3d_simd(
    a: &DenseMultiVector,
    b: &DenseMultiVector,
    table: &GeometricProductTable,
) -> DenseMultiVector {
    let config = a.config.clone();
    let mut a_arr = [0.0; 8];
    let mut b_arr = [0.0; 8];
    a_arr.copy_from_slice(&a.coefficients.as_slice().unwrap()[..8]);
    b_arr.copy_from_slice(&b.coefficients.as_slice().unwrap()[..8]);

    let result = geometric_product_simd_core_3d(&a_arr, &b_arr, table);
    DenseMultiVector::from_coefficients(config, result.to_vec())
}

/// 为 DenseMultiVector 提供的 SIMD 加速几何积（4D）
/// SIMD-accelerated geometric product for DenseMultiVector (4D)
#[cfg(feature = "simd")]
pub fn geometric_product_dense_4d_simd(
    a: &DenseMultiVector,
    b: &DenseMultiVector,
    table: &GeometricProductTable,
) -> DenseMultiVector {
    let config = a.config.clone();
    let mut a_arr = [0.0; 16];
    let mut b_arr = [0.0; 16];
    a_arr.copy_from_slice(&a.coefficients.as_slice().unwrap()[..16]);
    b_arr.copy_from_slice(&b.coefficients.as_slice().unwrap()[..16]);

    let result = geometric_product_simd_core_4d(&a_arr, &b_arr, table);
    DenseMultiVector::from_coefficients(config, result.to_vec())
}

/// 为 DenseMultiVector 提供的标量几何积（2D 优化版本）
/// Scalar geometric product for DenseMultiVector (2D optimized)
pub fn geometric_product_dense_2d_scalar(
    a: &DenseMultiVector,
    b: &DenseMultiVector,
    table: &GeometricProductTable,
) -> DenseMultiVector {
    let config = a.config.clone();
    let a_arr: [f64; 4] = [
        a.coefficients[0],
        a.coefficients[1],
        a.coefficients[2],
        a.coefficients[3],
    ];
    let b_arr: [f64; 4] = [
        b.coefficients[0],
        b.coefficients[1],
        b.coefficients[2],
        b.coefficients[3],
    ];

    let result = geometric_product_scalar_core_2d(&a_arr, &b_arr, table);
    DenseMultiVector::from_coefficients(config, result.to_vec())
}

/// 为 DenseMultiVector 提供的标量几何积（3D 优化版本）
/// Scalar geometric product for DenseMultiVector (3D optimized)
pub fn geometric_product_dense_3d_scalar(
    a: &DenseMultiVector,
    b: &DenseMultiVector,
    table: &GeometricProductTable,
) -> DenseMultiVector {
    let config = a.config.clone();
    let mut a_arr = [0.0; 8];
    let mut b_arr = [0.0; 8];
    a_arr.copy_from_slice(&a.coefficients.as_slice().unwrap()[..8]);
    b_arr.copy_from_slice(&b.coefficients.as_slice().unwrap()[..8]);

    let result = geometric_product_scalar_core_3d(&a_arr, &b_arr, table);
    DenseMultiVector::from_coefficients(config, result.to_vec())
}

/// 为 DenseMultiVector 提供的标量几何积（4D 优化版本）
/// Scalar geometric product for DenseMultiVector (4D optimized)
pub fn geometric_product_dense_4d_scalar(
    a: &DenseMultiVector,
    b: &DenseMultiVector,
    table: &GeometricProductTable,
) -> DenseMultiVector {
    let config = a.config.clone();
    let mut a_arr = [0.0; 16];
    let mut b_arr = [0.0; 16];
    a_arr.copy_from_slice(&a.coefficients.as_slice().unwrap()[..16]);
    b_arr.copy_from_slice(&b.coefficients.as_slice().unwrap()[..16]);

    let result = geometric_product_scalar_core_4d(&a_arr, &b_arr, table);
    DenseMultiVector::from_coefficients(config, result.to_vec())
}

/// 分发函数：根据维度选择 SIMD 或标量实现
/// Dispatcher: select SIMD or scalar implementation based on dimension
pub fn geometric_product_dense_optimized(
    a: &DenseMultiVector,
    b: &DenseMultiVector,
    table: &GeometricProductTable,
) -> DenseMultiVector {
    let dim = a.config.dimension();

    match dim {
        #[cfg(feature = "simd")]
        2 => geometric_product_dense_2d_simd(a, b, table),
        #[cfg(feature = "simd")]
        3 => geometric_product_dense_3d_simd(a, b, table),
        #[cfg(feature = "simd")]
        4 => geometric_product_dense_4d_simd(a, b, table),
        #[cfg(not(feature = "simd"))]
        2 => geometric_product_dense_2d_scalar(a, b, table),
        #[cfg(not(feature = "simd"))]
        3 => geometric_product_dense_3d_scalar(a, b, table),
        #[cfg(not(feature = "simd"))]
        4 => geometric_product_dense_4d_scalar(a, b, table),
        _ => super::geometric::geometric_product_dense(a, b, table),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::algebra_config::{AlgebraConfig, AlgebraConfigRef};
    use std::sync::Arc;

    fn test_config_2d() -> AlgebraConfigRef {
        Arc::new(AlgebraConfig::euclidean(2))
    }

    fn test_config_3d() -> AlgebraConfigRef {
        Arc::new(AlgebraConfig::euclidean(3))
    }

    fn test_config_4d() -> AlgebraConfigRef {
        Arc::new(AlgebraConfig::euclidean(4))
    }

    #[test]
    fn test_2d_geometric_product() {
        let config = test_config_2d();
        let e1 = DenseMultiVector::basis_vector(config.clone(), 0);
        let e2 = DenseMultiVector::basis_vector(config.clone(), 1);

        let table = GeometricProductTable::new(&config);

        let e1_e1 = geometric_product_dense_optimized(&e1, &e1, &table);
        assert!(
            (e1_e1.coefficients[0] - 1.0).abs() < 1e-10,
            "e1*e1 should be 1"
        );

        let e1_e2 = geometric_product_dense_optimized(&e1, &e2, &table);
        assert!(
            (e1_e2.coefficients[3] - 1.0).abs() < 1e-10,
            "e1*e2 should be e12"
        );

        let e2_e1 = geometric_product_dense_optimized(&e2, &e1, &table);
        assert!(
            (e2_e1.coefficients[3] + 1.0).abs() < 1e-10,
            "e2*e1 should be -e12"
        );
    }

    #[test]
    fn test_3d_geometric_product() {
        let config = test_config_3d();
        let e1 = DenseMultiVector::basis_vector(config.clone(), 0);
        let e2 = DenseMultiVector::basis_vector(config.clone(), 1);
        let e3 = DenseMultiVector::basis_vector(config.clone(), 2);

        let table = GeometricProductTable::new(&config);

        let e1_e1 = geometric_product_dense_optimized(&e1, &e1, &table);
        assert!(
            (e1_e1.coefficients[0] - 1.0).abs() < 1e-10,
            "e1*e1 should be 1"
        );

        let e1_e2 = geometric_product_dense_optimized(&e1, &e2, &table);
        assert!(
            (e1_e2.coefficients[3] - 1.0).abs() < 1e-10,
            "e1*e2 should be e12"
        );

        let e2_e1 = geometric_product_dense_optimized(&e2, &e1, &table);
        assert!(
            (e2_e1.coefficients[3] + 1.0).abs() < 1e-10,
            "e2*e1 should be -e12"
        );

        let e1_e3 = geometric_product_dense_optimized(&e1, &e3, &table);
        assert!(
            (e1_e3.coefficients[5] - 1.0).abs() < 1e-10,
            "e1*e3 should be e13"
        );
    }

    #[test]
    fn test_4d_geometric_product() {
        let config = test_config_4d();
        let e1 = DenseMultiVector::basis_vector(config.clone(), 0);
        let e2 = DenseMultiVector::basis_vector(config.clone(), 1);

        let table = GeometricProductTable::new(&config);

        let e1_e1 = geometric_product_dense_optimized(&e1, &e1, &table);
        assert!(
            (e1_e1.coefficients[0] - 1.0).abs() < 1e-10,
            "e1*e1 should be 1"
        );

        let e1_e2 = geometric_product_dense_optimized(&e1, &e2, &table);
        assert!(
            (e1_e2.coefficients[3] - 1.0).abs() < 1e-10,
            "e1*e2 should be e12"
        );
    }

    #[test]
    fn test_3d_scalar_vs_optimized_consistency() {
        let config = test_config_3d();

        let mut a = DenseMultiVector::zeros(config.clone());
        a.coefficients[0] = 1.0;
        a.coefficients[1] = 2.0;
        a.coefficients[2] = 3.0;
        a.coefficients[3] = 4.0;
        a.coefficients[4] = 5.0;
        a.coefficients[5] = 6.0;
        a.coefficients[6] = 7.0;
        a.coefficients[7] = 8.0;

        let mut b = DenseMultiVector::zeros(config.clone());
        b.coefficients[0] = 0.5;
        b.coefficients[1] = 1.5;
        b.coefficients[2] = 2.5;
        b.coefficients[3] = 3.5;
        b.coefficients[4] = 4.5;
        b.coefficients[5] = 5.5;
        b.coefficients[6] = 6.5;
        b.coefficients[7] = 7.5;

        let table = GeometricProductTable::new(&config);

        let result_scalar = geometric_product_dense_3d_scalar(&a, &b, &table);
        let result_optimized = geometric_product_dense_optimized(&a, &b, &table);

        for i in 0..8 {
            let diff = (result_scalar.coefficients[i] - result_optimized.coefficients[i]).abs();
            assert!(
                diff < 1e-10,
                "Coefficient {} differs: {} vs {}",
                i,
                result_scalar.coefficients[i],
                result_optimized.coefficients[i]
            );
        }
    }

    #[test]
    fn test_3d_vs_original_consistency() {
        let config = test_config_3d();
        let e1 = DenseMultiVector::basis_vector(config.clone(), 0);
        let e2 = DenseMultiVector::basis_vector(config.clone(), 1);

        let table = GeometricProductTable::new(&config);

        let result_optimized = geometric_product_dense_optimized(&e1, &e2, &table);
        let result_original = super::super::geometric::geometric_product_dense(&e1, &e2, &table);

        for i in 0..8 {
            let diff = (result_optimized.coefficients[i] - result_original.coefficients[i]).abs();
            assert!(
                diff < 1e-10,
                "Coefficient {} differs: {} vs {}",
                i,
                result_optimized.coefficients[i],
                result_original.coefficients[i]
            );
        }
    }
}
