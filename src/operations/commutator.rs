//! 交换子模块
//! 
//! 实现交换子运算：A × B = ½(AB - BA)
//! 
//! 交换子满足：
//! - 反对称性：A × B = -B × A
//! - Jacobi 恒等式：A × (B × C) + B × (C × A) + C × (A × B) = 0
//! - 与二重向量的交换子是保持阶次的运算

use crate::multivector::MultiVector;

/// 计算交换子
pub fn commutator(a: &MultiVector, b: &MultiVector) -> MultiVector {
    // A × B = ½(AB - BA)
    let ab = a.geometric_product(b);
    let ba = b.geometric_product(a);
    let diff = ab.sub(&ba);
    diff.scale(0.5)
}

impl MultiVector {
    /// 交换子 A × B
    pub fn commutator(&self, other: &Self) -> Self {
        commutator(self, other)
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
    fn test_commutator_scalars() {
        let config = test_config();
        let a = MultiVector::from_scalar(config.clone(), 2.0);
        let b = MultiVector::from_scalar(config.clone(), 3.0);
        
        // 标量交换子为 0
        let result = a.commutator(&b);
        assert!(result.is_zero());
    }

    #[test]
    fn test_commutator_vectors() {
        let config = test_config();
        let e1 = MultiVector::basis_vector(config.clone(), 0);
        let e2 = MultiVector::basis_vector(config.clone(), 1);
        
        // e1 × e2 = ½(e1e2 - e2e1) = ½(e12 - (-e12)) = e12
        let result = e1.commutator(&e2);
        assert!((result.get_coefficient(3) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_commutator_antisymmetry() {
        let config = test_config();
        let e1 = MultiVector::basis_vector(config.clone(), 0);
        let e2 = MultiVector::basis_vector(config.clone(), 1);
        
        // e1 × e2 = -(e2 × e1)
        let result1 = e1.commutator(&e2);
        let result2 = e2.commutator(&e1);
        
        let sum = result1.add(&result2);
        assert!(sum.is_zero());
    }

    #[test]
    fn test_commutator_bivector_vector() {
        let config = test_config();
        let e1 = MultiVector::basis_vector(config.clone(), 0);
        let e2 = MultiVector::basis_vector(config.clone(), 1);
        let e12 = e1.outer_product(&e2);
        
        // e12 × e1 = ½(e12e1 - e1e12)
        // e12e1 = e1e2e1 = -e1e1e2 = -e2
        // e1e12 = e1e1e2 = e2
        // 所以 e12 × e1 = ½(-e2 - e2) = -e2
        let result = e12.commutator(&e1);
        assert!((result.get_coefficient(2) + 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_commutator_same_element() {
        let config = test_config();
        let e1 = MultiVector::basis_vector(config.clone(), 0);
        
        // A × A = 0
        let result = e1.commutator(&e1);
        assert!(result.is_zero());
    }
}
