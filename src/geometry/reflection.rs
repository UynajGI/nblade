//! 反射模块
//!
//! 实现多重向量在子空间中的反射
//!
//! 公式：B' = Ar B*ʳ Ar⁻¹
//! 其中 r 是 Ar 的阶次，*ʳ 表示 r 次阶次对合

use crate::multivector::MultiVector;

/// 反射多重向量
///
/// 公式：B' = Ar B*ʳ Ar⁻¹
pub fn reflect_in(mv: &MultiVector, blade: &MultiVector) -> Result<MultiVector, String> {
    // 获取 blade 的阶次
    let blade_grade = get_blade_grade(blade);

    // 计算 blade 的逆
    let blade_inv = blade.inverse()?;

    // 计算 B*ʳ (r 次阶次对合)
    let mv_transformed = grade_involution_power(mv, blade_grade);

    // Ar B*ʳ
    let left_product = blade.geometric_product(&mv_transformed);

    // (Ar B*ʳ) Ar⁻¹
    let result = left_product.geometric_product(&blade_inv);

    Ok(result)
}

/// 获取 blade 的阶次
fn get_blade_grade(blade: &MultiVector) -> u32 {
    // 找到第一个非零系数的阶次
    for i in 0..blade.config().basis_count() {
        if blade.get_coefficient(i).abs() > 1e-15 {
            return blade.config().grade_of_index(i) as u32;
        }
    }
    0
}

/// 计算 r 次阶次对合
/// 反射公式：对于 grade(A) = r
/// - r 为奇数：不应用对合 (X' = A X A⁻¹)
/// - r 为偶数：应用一次对合 (X' = A X* A⁻¹)
fn grade_involution_power(mv: &MultiVector, r: u32) -> MultiVector {
    if r % 2 == 0 {
        // 偶数阶 blade：应用一次阶次对合
        mv.grade_involution()
    } else {
        // 奇数阶 blade：不变
        mv.clone()
    }
}

impl MultiVector {
    /// 在由 blade 表示的子空间中反射
    pub fn reflect_in(&self, blade: &Self) -> Result<Self, String> {
        reflect_in(self, blade)
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
    fn test_reflect_vector_in_vector() {
        let config = test_config();
        let e1 = MultiVector::basis_vector(config.clone(), 0);
        let e2 = MultiVector::basis_vector(config.clone(), 1);

        // e2 在 e1 中反射 = -e2
        let reflected = e2.reflect_in(&e1).unwrap();
        assert!((reflected.get_coefficient(2) + 1.0).abs() < 1e-10);

        // e1 在 e1 中反射 = e1
        let reflected = e1.reflect_in(&e1).unwrap();
        assert!((reflected.get_coefficient(1) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_reflect_vector_in_plane() {
        let config = test_config();
        let e1 = MultiVector::basis_vector(config.clone(), 0);
        let e2 = MultiVector::basis_vector(config.clone(), 1);
        let e3 = MultiVector::basis_vector(config.clone(), 2);
        let e12 = e1.outer_product(&e2);

        // e3 在 e12 平面中反射 = -e3 (垂直于平面)
        // e3 has index 4 (0b100)
        let reflected = e3.reflect_in(&e12).unwrap();
        assert!((reflected.get_coefficient(4) + 1.0).abs() < 1e-10);

        // e1 在 e12 平面中反射 = e1 (在平面内)
        let reflected = e1.reflect_in(&e12).unwrap();
        assert!((reflected.get_coefficient(1) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_reflect_twice() {
        let config = test_config();
        let e1 = MultiVector::basis_vector(config.clone(), 0);
        let e2 = MultiVector::basis_vector(config.clone(), 1);

        // 反射两次应该回到原向量
        let reflected1 = e2.reflect_in(&e1).unwrap();
        let reflected2 = reflected1.reflect_in(&e1).unwrap();

        assert!((reflected2.get_coefficient(2) - 1.0).abs() < 1e-10);
    }
}
