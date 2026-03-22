//! 旋转模块
//!
//! 实现使用转子（Rotor）的旋转
//!
//! 公式：B' = R B R⁻¹ = R B R†
//! 其中 R 是转子，满足 R†R = 1
//!
//! 转子形式：R = exp(-Bθ/2) = cos(θ/2) - B sin(θ/2)
//! 其中 B 是单位二重向量（旋转平面）

use crate::multivector::MultiVector;

/// 创建转子
///
/// # 参数
/// * `plane` - 旋转平面（二重向量）
/// * `angle` - 旋转角度（弧度）
///
/// # 返回
/// 转子 R = cos(θ/2) - B sin(θ/2)
pub fn create_rotor(plane: &MultiVector, angle: f64) -> MultiVector {
    let config = plane.config().clone();

    // 归一化平面
    let plane_norm = plane.norm();
    let unit_plane = if plane_norm.abs() > 1e-15 {
        plane.scale(1.0 / plane_norm)
    } else {
        plane.clone()
    };

    // R = cos(θ/2) - B sin(θ/2)
    let half_angle = angle / 2.0;
    let cos_half = half_angle.cos();
    let sin_half = half_angle.sin();

    let scalar_part = MultiVector::from_scalar(config.clone(), cos_half);
    let bivector_part = unit_plane.scale(-sin_half);

    scalar_part.add(&bivector_part)
}

/// 使用转子旋转多重向量
///
/// 公式：B' = R B R⁻¹
pub fn rotate_by(mv: &MultiVector, rotor: &MultiVector) -> Result<MultiVector, String> {
    // 计算转子逆
    let rotor_inv = rotor.inverse()?;

    // R B R⁻¹
    let temp = rotor.geometric_product(mv);
    let result = temp.geometric_product(&rotor_inv);

    Ok(result)
}

/// 在平面中旋转向量
///
/// # 参数
/// * `vector` - 要旋转的向量
/// * `plane` - 旋转平面（二重向量）
/// * `angle` - 旋转角度（弧度）
pub fn rotate_vector_in_plane(
    vector: &MultiVector,
    plane: &MultiVector,
    angle: f64,
) -> Result<MultiVector, String> {
    let rotor = create_rotor(plane, angle);
    rotate_by(vector, &rotor)
}

impl MultiVector {
    /// 使用转子旋转
    pub fn rotate_by(&self, rotor: &Self) -> Result<Self, String> {
        rotate_by(self, rotor)
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
    fn test_rotor_creation() {
        let config = test_config();
        let e1 = MultiVector::basis_vector(config.clone(), 0);
        let e2 = MultiVector::basis_vector(config.clone(), 1);
        let e12 = e1.outer_product(&e2);

        // 创建 90 度转子
        let rotor = create_rotor(&e12, std::f64::consts::PI / 2.0);

        // 应该有标量和二重向量部分
        assert!(rotor.get_coefficient(0).abs() > 1e-10);
        assert!(rotor.get_coefficient(3).abs() > 1e-10);
    }

    #[test]
    fn test_rotate_vector_90_degrees() {
        let config = test_config();
        let e1 = MultiVector::basis_vector(config.clone(), 0);
        let e2 = MultiVector::basis_vector(config.clone(), 1);
        let e12 = e1.outer_product(&e2);

        // 旋转 90 度
        let rotor = create_rotor(&e12, std::f64::consts::PI / 2.0);
        let rotated = e1.rotate_by(&rotor).unwrap();

        // e1 旋转 90 度应该变成 e2
        assert!((rotated.get_coefficient(2) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_rotate_vector_180_degrees() {
        let config = test_config();
        let e1 = MultiVector::basis_vector(config.clone(), 0);
        let e2 = MultiVector::basis_vector(config.clone(), 1);
        let e12 = e1.outer_product(&e2);

        // 旋转 180 度
        let rotor = create_rotor(&e12, std::f64::consts::PI);
        let rotated = e1.rotate_by(&rotor).unwrap();

        // e1 旋转 180 度应该变成 -e1
        assert!((rotated.get_coefficient(1) + 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_rotate_preserves_norm() {
        let config = test_config();
        let e1 = MultiVector::basis_vector(config.clone(), 0);
        let e2 = MultiVector::basis_vector(config.clone(), 1);
        let e12 = e1.outer_product(&e2);

        let rotor = create_rotor(&e12, std::f64::consts::PI / 4.0);
        let rotated = e1.rotate_by(&rotor).unwrap();

        // 旋转应该保持范数
        assert!((rotated.norm() - e1.norm()).abs() < 1e-10);
    }

    #[test]
    fn test_rotate_vector_out_of_plane() {
        let config = test_config();
        let e1 = MultiVector::basis_vector(config.clone(), 0);
        let e2 = MultiVector::basis_vector(config.clone(), 1);
        let e3 = MultiVector::basis_vector(config.clone(), 2);
        let e12 = e1.outer_product(&e2);

        // 旋转 e3（垂直于旋转平面）应该不变
        let rotor = create_rotor(&e12, std::f64::consts::PI / 2.0);
        let rotated = e3.rotate_by(&rotor).unwrap();

        // e3 has index 4 (0b100)
        assert!((rotated.get_coefficient(4) - 1.0).abs() < 1e-10);
    }
}
