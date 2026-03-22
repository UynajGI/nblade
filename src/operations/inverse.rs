//! 逆运算模块
//! 
//! 实现多向量的逆：A⁻¹ = A† / |A|²
//! 
//! 注意：并非所有多向量都有逆，只有非零范数的多向量才可逆

use crate::multivector::{DenseMultiVector, SparseMultiVector, MultiVector};
use crate::operations::{involution::reversion_dense, norm::{norm_squared_dense, norm_squared_sparse}};

/// 密集多向量的逆
pub fn inverse_dense(mv: &DenseMultiVector) -> Result<DenseMultiVector, String> {
    let norm_sq = norm_squared_dense(mv);
    
    if norm_sq.abs() < 1e-15 {
        return Err(format!("MultiVector has zero norm ({}), cannot invert", norm_sq));
    }
    
    let mv_rev = reversion_dense(mv);
    let result = mv_rev.scale(1.0 / norm_sq);
    
    Ok(result)
}

/// 稀疏多向量的逆
pub fn inverse_sparse(mv: &SparseMultiVector) -> Result<SparseMultiVector, String> {
    let norm_sq = norm_squared_sparse(mv);
    
    if norm_sq.abs() < 1e-15 {
        return Err(format!("MultiVector has zero norm ({}), cannot invert", norm_sq));
    }
    
    let mv_rev = crate::operations::involution::reversion_sparse(mv);
    let result = mv_rev.scale(1.0 / norm_sq);
    
    Ok(result)
}

impl MultiVector {
    /// 逆 A⁻¹
    pub fn inverse(&self) -> Result<Self, String> {
        match self {
            MultiVector::Dense(d) => inverse_dense(d).map(MultiVector::Dense),
            MultiVector::Sparse(s) => inverse_sparse(s).map(MultiVector::Sparse),
        }
    }

    /// 检查是否可逆
    pub fn is_invertible(&self) -> bool {
        self.norm_squared().abs() > 1e-15
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
    fn test_inverse_scalar() {
        let config = test_config();
        let scalar = MultiVector::from_scalar(config.clone(), 2.0);
        
        let inv = scalar.inverse().unwrap();
        assert!((inv.get_coefficient(0) - 0.5).abs() < 1e-10);
        
        // scalar * inverse = 1
        let product = scalar.geometric_product(&inv);
        assert!((product.get_coefficient(0) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_inverse_vector() {
        let config = test_config();
        let e1 = MultiVector::basis_vector(config.clone(), 0);
        
        let inv = e1.inverse().unwrap();
        // e1⁻¹ = e1 (因为 e1² = 1)
        assert!((inv.get_coefficient(1) - 1.0).abs() < 1e-10);
        
        // e1 * e1⁻¹ = 1
        let product = e1.geometric_product(&inv);
        assert!((product.get_coefficient(0) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_inverse_bivector() {
        let config = test_config();
        let e1 = MultiVector::basis_vector(config.clone(), 0);
        let e2 = MultiVector::basis_vector(config.clone(), 1);
        let e12 = e1.outer_product(&e2);
        
        let inv = e12.inverse().unwrap();
        // (e1∧e2)⁻¹ = -(e1∧e2) (因为 (e1∧e2)² = -1)
        assert!((inv.get_coefficient(3) + 1.0).abs() < 1e-10);
        
        // e12 * inv = 1
        let product = e12.geometric_product(&inv);
        assert!((product.get_coefficient(0) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_inverse_zero() {
        let config = test_config();
        let zero = MultiVector::zeros(config.clone());
        
        assert!(zero.inverse().is_err());
    }

    #[test]
    fn test_is_invertible() {
        let config = test_config();
        let e1 = MultiVector::basis_vector(config.clone(), 0);
        let zero = MultiVector::zeros(config.clone());
        
        assert!(e1.is_invertible());
        assert!(!zero.is_invertible());
    }
}
