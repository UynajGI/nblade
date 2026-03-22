//! 标架和互逆标架模块 / Frame and Reciprocal Frame Module
//!
//! 实现几何代数中互逆标架的计算
//!
//! # 公式 / Formulas
//!
//! 给定标架 {a₁, a₂, ..., aₙ}，互逆标架 {a¹, a², ..., aⁿ} 定义为：
//!
//! aⁱ = (-1)^(i-1) (a₁ ∧ a₂ ∧ ... ∧ âᵢ ∧ ... ∧ aₙ) a_N⁻¹
//!
//! 其中：
//! - âᵢ 表示 aᵢ 被省略
//! - a_N = a₁ ∧ a₂ ∧ ... ∧ aₙ 是体积元素
//! - a_N⁻¹ 是体积元素的逆
//!
//! 互逆标架满足：
//! aⁱ⌋aⱼ = δⁱⱼ (Kronecker delta)

use crate::multivector::MultiVector;

/// 计算互逆标架 / Compute reciprocal frame
///
/// 给定一个标架 {a₁, a₂, ..., aₙ}，计算互逆标架 {a¹, a², ..., aⁿ}
/// 使得 aⁱ⌋aⱼ = δⁱⱼ
///
/// # 参数 / Arguments
///
/// * `vectors` - 标架向量切片，必须是非空的 1-向量
///
/// # 返回 / Returns
///
/// * `Ok(Vec<MultiVector>)` - 互逆标架向量
/// * `Err(String)` - 错误信息（标架为空、向量数量不匹配维度、体积元素不可逆等）
///
/// # 公式 / Formula
///
/// aⁱ = (-1)^(i-1) (a₁ ∧ a₂ ∧ ... ∧ âᵢ ∧ ... ∧ aₙ) a_N⁻¹
///
/// # 示例 / Example
///
/// ```rust
/// use nblade::{AlgebraConfig, MultiVector};
/// use nblade::basis::reciprocal_frame;
/// use std::sync::Arc;
///
/// let config = Arc::new(AlgebraConfig::euclidean(3));
///
/// // 正交基的互逆标架就是它自己
/// let e1 = MultiVector::basis_vector(config.clone(), 0);
/// let e2 = MultiVector::basis_vector(config.clone(), 1);
/// let e3 = MultiVector::basis_vector(config.clone(), 2);
///
/// let reciprocal = reciprocal_frame(&[e1.clone(), e2.clone(), e3.clone()]).unwrap();
///
/// // 验证 aⁱ⌋aⱼ = δⁱⱼ
/// for i in 0..3 {
///     for j in 0..3 {
///         let contraction = reciprocal[i].left_inner(&[e1.clone(), e2.clone(), e3.clone()][j]);
///         let expected = if i == j { 1.0 } else { 0.0 };
///         assert!((contraction.scalar_part() - expected).abs() < 1e-10);
///     }
/// }
/// ```
pub fn reciprocal_frame(vectors: &[MultiVector]) -> Result<Vec<MultiVector>, String> {
    // 验证输入
    if vectors.is_empty() {
        return Err("Frame cannot be empty".to_string());
    }

    // 获取代数配置
    let config = vectors[0].config().clone();
    let n = vectors.len();

    // 验证向量数量与维度匹配
    let dimension = config.dimension() as usize;
    if n != dimension {
        return Err(format!(
            "Frame size {} must match algebra dimension {}",
            n, dimension
        ));
    }

    // 验证所有向量使用相同的配置
    for (i, v) in vectors.iter().enumerate() {
        if v.config() != &config {
            return Err(format!("Vector {} has different algebra configuration", i));
        }
    }

    // 计算体积元素 a_N = a₁ ∧ a₂ ∧ ... ∧ aₙ
    let volume_element = compute_volume_element(vectors)?;

    // 计算体积元素的逆 a_N⁻¹
    let volume_inverse = volume_element
        .inverse()
        .map_err(|e| format!("Cannot invert volume element: {}", e))?;

    // 计算每个互逆向量
    let mut reciprocal = Vec::with_capacity(n);

    for i in 0..n {
        // 计算除了第 i 个向量外所有向量的外积
        let wedge_except_i = compute_wedge_except(vectors, i)?;

        // 应用符号因子 (-1)^(i-1)，注意 i 是 0-indexed
        let sign = if i % 2 == 0 { 1.0 } else { -1.0 }; // (-1)^i 等价于 (-1)^(i+1-1) = (-1)^i

        // aⁱ = (-1)^i × (wedge_except_i) × (volume_inverse)
        // 注意：公式中是 (-1)^(i-1)，但 i 是 1-indexed
        // 我们的 i 是 0-indexed，所以是 (-1)^i
        let scaled_wedge = wedge_except_i.scale(sign);
        let reciprocal_vector = scaled_wedge.geometric_product(&volume_inverse);

        reciprocal.push(reciprocal_vector);
    }

    Ok(reciprocal)
}

/// 计算体积元素（所有向量的外积）
///
/// a_N = a₁ ∧ a₂ ∧ ... ∧ aₙ
fn compute_volume_element(vectors: &[MultiVector]) -> Result<MultiVector, String> {
    if vectors.is_empty() {
        return Err("Cannot compute volume element of empty frame".to_string());
    }

    let mut result = vectors[0].clone();
    for v in vectors.iter().skip(1) {
        result = result.outer_product(v);
    }

    // 检查体积元素是否为零（线性相关向量）
    if result.is_zero() {
        return Err("Frame vectors are linearly dependent (volume element is zero)".to_string());
    }

    Ok(result)
}

/// 计算除了第 i 个向量外所有向量的外积
///
/// 返回 a₁ ∧ ... ∧ âᵢ ∧ ... ∧ aₙ（省略第 i 个向量）
fn compute_wedge_except(
    vectors: &[MultiVector],
    except_index: usize,
) -> Result<MultiVector, String> {
    if vectors.len() <= 1 {
        // 如果只有一个向量，省略后得到单位标量
        let config = vectors[0].config().clone();
        return Ok(MultiVector::one(config));
    }

    // 找到第一个不是 except_index 的向量作为起始
    let first_idx = if except_index == 0 { 1 } else { 0 };
    let mut result = vectors[first_idx].clone();

    for (i, v) in vectors.iter().enumerate() {
        if i == except_index || i == first_idx {
            continue;
        }
        result = result.outer_product(v);
    }

    // 如果所有向量都被省略了（不可能，因为已检查 len > 1）
    // 或者结果为零
    if result.is_zero() && vectors.len() > 2 {
        // 这通常表示向量线性相关
        return Err(format!(
            "Wedge product excluding vector {} resulted in zero (possibly linearly dependent)",
            except_index
        ));
    }

    Ok(result)
}

/// 验证互逆标架条件
///
/// 检查 aⁱ⌋aⱼ = δⁱⱼ 是否成立
///
/// # 参数 / Arguments
///
/// * `frame` - 原始标架
/// * `reciprocal` - 互逆标架
/// * `tolerance` - 容差（默认 1e-10）
///
/// # 返回 / Returns
///
/// * `Ok(())` - 验证通过
/// * `Err(String)` - 验证失败，包含错误详情
pub fn verify_reciprocal_frame(
    frame: &[MultiVector],
    reciprocal: &[MultiVector],
    tolerance: f64,
) -> Result<(), String> {
    if frame.len() != reciprocal.len() {
        return Err(format!(
            "Frame size {} does not match reciprocal size {}",
            frame.len(),
            reciprocal.len()
        ));
    }

    let n = frame.len();
    let mut errors = Vec::new();

    for i in 0..n {
        for j in 0..n {
            let contraction = reciprocal[i].left_inner(&frame[j]);
            let expected = if i == j { 1.0 } else { 0.0 };
            let actual = contraction.scalar_part();

            if (actual - expected).abs() > tolerance {
                errors.push(format!(
                    "a^{}⌋a_{} = {} (expected {})",
                    i + 1,
                    j + 1,
                    actual,
                    expected
                ));
            }
        }
    }

    if errors.is_empty() {
        Ok(())
    } else {
        Err(format!(
            "Reciprocal frame verification failed:\n{}",
            errors.join("\n")
        ))
    }
}

/// 计算标架的度量张量
///
/// 度量张量 g_ij = a_i · a_j（标量积）
///
/// # 参数 / Arguments
///
/// * `vectors` - 标架向量
///
/// # 返回 / Returns
///
/// 度量张量矩阵（n × n）
pub fn metric_tensor(vectors: &[MultiVector]) -> Result<Vec<Vec<f64>>, String> {
    if vectors.is_empty() {
        return Err("Cannot compute metric tensor of empty frame".to_string());
    }

    let n = vectors.len();
    let mut g = vec![vec![0.0; n]; n];

    for i in 0..n {
        for j in 0..n {
            // 对于向量，a·b = a⌋b
            let inner = vectors[i].left_inner(&vectors[j]);
            g[i][j] = inner.scalar_part();
        }
    }

    Ok(g)
}

// ============================================================================
// 基展开 / Basis Expansion
// ============================================================================

use crate::basis::index::BasisIndex;

/// 多向量基展开 / Multivector Basis Expansion
///
/// 将多向量展开为基刃的线性组合：A = Σᴵ Aᴵ aᴵ
/// 返回各基刃的系数 Aᴵ = A * aᴵ（标量积）
///
/// # 公式 / Formula
///
/// A = Σᴵ Aᴵ aᴵ 其中 Aᴵ = A * aᴵ (文档公式 123)
///
/// 对于正交归一基，aᴵ = a_I，因此：
/// Aᴵ = A * a_I（标量积）
///
/// # 参数 / Arguments
///
/// * `mv` - 待展开的多向量
///
/// # 返回 / Returns
///
/// Vec<(BasisIndex, f64)> - 基索引到系数的映射（仅包含非零系数）
///
/// # 示例 / Example
///
/// ```rust
/// use nblade::{AlgebraConfig, MultiVector};
/// use nblade::basis::basis_expansion;
/// use std::sync::Arc;
///
/// let config = Arc::new(AlgebraConfig::euclidean(3));
///
/// // 标量展开
/// let scalar = MultiVector::from_scalar(config.clone(), 5.0);
/// let coeffs = basis_expansion(&scalar);
/// assert_eq!(coeffs.len(), 1);
/// assert_eq!(coeffs[0].0, 0); // 标量索引
/// assert!((coeffs[0].1 - 5.0).abs() < 1e-10);
///
/// // 向量展开
/// let e1 = MultiVector::basis_vector(config.clone(), 0);
/// let coeffs = basis_expansion(&e1);
/// assert_eq!(coeffs.len(), 1);
/// assert_eq!(coeffs[0].0, 1); // e₁ 的索引
/// assert!((coeffs[0].1 - 1.0).abs() < 1e-10);
/// ```
pub fn basis_expansion(mv: &MultiVector) -> Vec<(BasisIndex, f64)> {
    let config = mv.config();
    let basis_count = config.basis_count();

    let mut coefficients = Vec::new();

    // 遍历所有可能的基刃索引
    for index in 0..basis_count as BasisIndex {
        // 创建该基刃
        let basis_blade = create_basis_blade(config.clone(), index);

        // 计算系数 Aᴵ = A * a_I（标量积）
        let coef = mv.scalar_product(&basis_blade);

        // 只保留非零系数
        if coef.abs() > 1e-15 {
            coefficients.push((index, coef));
        }
    }

    coefficients
}

/// 使用指定标架的基展开 / Basis Expansion with Specified Frame
///
/// 给定标架 {a_I} 和互逆标架 {aᴵ}，计算多向量在该标架下的展开系数
///
/// # 公式 / Formula
///
/// A = Σᴵ Aᴵ aᴵ 其中 Aᴵ = A * aᴵ (文档公式 123)
///
/// # 参数 / Arguments
///
/// * `mv` - 待展开的多向量
/// * `reciprocal_blades` - 互逆标架的基刃切片，长度必须等于代数基数量
///
/// # 返回 / Returns
///
/// Vec<(BasisIndex, f64)> - 基索引到系数的映射
pub fn basis_expansion_with_reciprocal(
    mv: &MultiVector,
    reciprocal_blades: &[MultiVector],
) -> Vec<(BasisIndex, f64)> {
    let config = mv.config();
    let basis_count = config.basis_count();

    assert_eq!(
        reciprocal_blades.len(),
        basis_count,
        "Reciprocal blades count {} must match basis count {}",
        reciprocal_blades.len(),
        basis_count
    );

    let mut coefficients = Vec::new();

    for (index, reciprocal_blade) in reciprocal_blades.iter().enumerate() {
        // Aᴵ = A * aᴵ（标量积）
        let coef = mv.scalar_product(reciprocal_blade);

        if coef.abs() > 1e-15 {
            coefficients.push((index as BasisIndex, coef));
        }
    }

    coefficients
}

/// 从展开系数重构多向量 / Reconstruct Multivector from Expansion
///
/// 给定基展开系数，重构原始多向量
///
/// # 公式 / Formula
///
/// A = Σᴵ Aᴵ aᴵ
///
/// # 参数 / Arguments
///
/// * `config` - 代数配置
/// * `coefficients` - 基索引到系数的映射
///
/// # 返回 / Returns
///
/// 重构的多向量
pub fn basis_reconstruction(
    config: crate::algebra_config::AlgebraConfigRef,
    coefficients: &[(BasisIndex, f64)],
) -> MultiVector {
    let mut result = MultiVector::zeros(config.clone());

    for &(index, coef) in coefficients {
        let basis_blade = create_basis_blade(config.clone(), index);
        result = result.add(&basis_blade.scale(coef));
    }

    result
}

/// 创建指定索引的基刃 / Create Basis Blade at Specified Index
///
/// 基索引编码：位 i 置位表示包含 e_{i+1}
///
/// # 示例 / Examples
///
/// - 索引 0 = 标量 1
/// - 索引 1 = e₁ (0b001)
/// - 索引 2 = e₂ (0b010)
/// - 索引 3 = e₁∧e₂ (0b011)
/// - 索引 7 = e₁∧e₂∧e₃ (0b111)
fn create_basis_blade(
    config: crate::algebra_config::AlgebraConfigRef,
    index: BasisIndex,
) -> MultiVector {
    if index == 0 {
        return MultiVector::one(config);
    }

    // 找到所有置位的维度索引
    let mut blade = MultiVector::one(config.clone());
    let mut first = true;

    for bit in 0..config.dimension() {
        if index & (1u64 << bit) != 0 {
            let basis_vec = MultiVector::basis_vector(config.clone(), bit);
            if first {
                blade = basis_vec;
                first = false;
            } else {
                blade = blade.outer_product(&basis_vec);
            }
        }
    }

    blade
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{AlgebraConfig, AlgebraConfigRef, MultiVector};
    use std::sync::Arc;

    fn test_config() -> AlgebraConfigRef {
        Arc::new(AlgebraConfig::euclidean(3))
    }

    #[test]
    fn test_reciprocal_frame_orthonormal() {
        // 正交基的互逆标架就是它自己
        let config = test_config();
        let e1 = MultiVector::basis_vector(config.clone(), 0);
        let e2 = MultiVector::basis_vector(config.clone(), 1);
        let e3 = MultiVector::basis_vector(config.clone(), 2);

        let frame = vec![e1.clone(), e2.clone(), e3.clone()];
        let reciprocal = reciprocal_frame(&frame).unwrap();

        // 验证 aⁱ⌋aⱼ = δⁱⱼ
        for i in 0..3 {
            for j in 0..3 {
                let contraction = reciprocal[i].left_inner(&frame[j]);
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(
                    (contraction.scalar_part() - expected).abs() < 1e-10,
                    "a^{}⌋a_{} = {} (expected {})",
                    i + 1,
                    j + 1,
                    contraction.scalar_part(),
                    expected
                );
            }
        }
    }

    #[test]
    fn test_reciprocal_frame_non_orthogonal_2d() {
        // 非正交基的 2D 情况
        let config = Arc::new(AlgebraConfig::euclidean(2));

        // 构造两个非正交向量
        // a₁ = e₁
        // a₂ = e₁ + e₂
        let e1 = MultiVector::basis_vector(config.clone(), 0);
        let e2 = MultiVector::basis_vector(config.clone(), 1);

        let a1 = e1.clone();
        let a2 = e1.add(&e2);

        let frame = vec![a1.clone(), a2.clone()];
        let reciprocal = reciprocal_frame(&frame).unwrap();

        // 验证 aⁱ⌋aⱼ = δⁱⱼ
        for i in 0..2 {
            for j in 0..2 {
                let contraction = reciprocal[i].left_inner(&frame[j]);
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(
                    (contraction.scalar_part() - expected).abs() < 1e-10,
                    "a^{}⌋a_{} = {} (expected {})",
                    i + 1,
                    j + 1,
                    contraction.scalar_part(),
                    expected
                );
            }
        }
    }

    #[test]
    fn test_reciprocal_frame_scaled_vectors() {
        // 缩放向量的互逆标架
        let config = test_config();

        // a₁ = 2e₁, a₂ = 3e₂, a₃ = e₃
        let e1 = MultiVector::basis_vector(config.clone(), 0);
        let e2 = MultiVector::basis_vector(config.clone(), 1);
        let e3 = MultiVector::basis_vector(config.clone(), 2);

        let a1 = e1.scale(2.0);
        let a2 = e2.scale(3.0);
        let a3 = e3;

        let frame = vec![a1.clone(), a2.clone(), a3.clone()];
        let reciprocal = reciprocal_frame(&frame).unwrap();

        // 验证 aⁱ⌋aⱼ = δⁱⱼ
        for i in 0..3 {
            for j in 0..3 {
                let contraction = reciprocal[i].left_inner(&frame[j]);
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(
                    (contraction.scalar_part() - expected).abs() < 1e-10,
                    "a^{}⌋a_{} = {} (expected {})",
                    i + 1,
                    j + 1,
                    contraction.scalar_part(),
                    expected
                );
            }
        }

        // 互逆向量应该分别是 (1/2)e₁, (1/3)e₂, e₃
        assert!((reciprocal[0].get_coefficient(1) - 0.5).abs() < 1e-10); // (1/2)e₁
        assert!((reciprocal[1].get_coefficient(2) - 1.0 / 3.0).abs() < 1e-10); // (1/3)e₂
        assert!((reciprocal[2].get_coefficient(4) - 1.0).abs() < 1e-10); // e₃
    }

    #[test]
    fn test_reciprocal_frame_4d() {
        // 4D 正交基
        let config = Arc::new(AlgebraConfig::euclidean(4));

        let e1 = MultiVector::basis_vector(config.clone(), 0);
        let e2 = MultiVector::basis_vector(config.clone(), 1);
        let e3 = MultiVector::basis_vector(config.clone(), 2);
        let e4 = MultiVector::basis_vector(config.clone(), 3);

        let frame = vec![e1.clone(), e2.clone(), e3.clone(), e4.clone()];
        let reciprocal = reciprocal_frame(&frame).unwrap();

        // 验证 aⁱ⌋aⱼ = δⁱⱼ
        verify_reciprocal_frame(&frame, &reciprocal, 1e-10).unwrap();
    }

    #[test]
    fn test_reciprocal_frame_empty_error() {
        let result = reciprocal_frame(&[]);
        assert!(result.is_err());
        assert!(result.unwrap_err().contains("empty"));
    }

    #[test]
    fn test_reciprocal_frame_dimension_mismatch() {
        let config = test_config();
        let e1 = MultiVector::basis_vector(config.clone(), 0);
        let e2 = MultiVector::basis_vector(config.clone(), 1);

        // 只提供 2 个向量，但维度是 3
        let result = reciprocal_frame(&[e1, e2]);
        assert!(result.is_err());
        assert!(result.unwrap_err().contains("must match"));
    }

    #[test]
    fn test_verify_reciprocal_frame() {
        let config = test_config();
        let e1 = MultiVector::basis_vector(config.clone(), 0);
        let e2 = MultiVector::basis_vector(config.clone(), 1);
        let e3 = MultiVector::basis_vector(config.clone(), 2);

        let frame = vec![e1.clone(), e2.clone(), e3.clone()];
        let reciprocal = reciprocal_frame(&frame).unwrap();

        // 正确的互逆标架应该验证通过
        assert!(verify_reciprocal_frame(&frame, &reciprocal, 1e-10).is_ok());

        // 错误的互逆标架应该验证失败
        let wrong_reciprocal = vec![e1.clone(), e1.clone(), e3.clone()];
        assert!(verify_reciprocal_frame(&frame, &wrong_reciprocal, 1e-10).is_err());
    }

    #[test]
    fn test_metric_tensor_orthonormal() {
        let config = test_config();
        let e1 = MultiVector::basis_vector(config.clone(), 0);
        let e2 = MultiVector::basis_vector(config.clone(), 1);
        let e3 = MultiVector::basis_vector(config.clone(), 2);

        let frame = vec![e1, e2, e3];
        let g = metric_tensor(&frame).unwrap();

        // 正交基的度量张量是单位矩阵
        for i in 0..3 {
            for j in 0..3 {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!((g[i][j] - expected).abs() < 1e-10);
            }
        }
    }

    #[test]
    fn test_metric_tensor_non_orthogonal() {
        let config = Arc::new(AlgebraConfig::euclidean(2));
        let e1 = MultiVector::basis_vector(config.clone(), 0);
        let e2 = MultiVector::basis_vector(config.clone(), 1);

        // a₁ = e₁, a₂ = e₁ + e₂
        let a1 = e1.clone();
        let a2 = e1.add(&e2);

        let frame = vec![a1, a2];
        let g = metric_tensor(&frame).unwrap();

        // g₁₁ = a₁·a₁ = 1
        // g₁₂ = a₁·a₂ = 1
        // g₂₁ = a₂·a₁ = 1
        // g₂₂ = a₂·a₂ = 2
        assert!((g[0][0] - 1.0).abs() < 1e-10);
        assert!((g[0][1] - 1.0).abs() < 1e-10);
        assert!((g[1][0] - 1.0).abs() < 1e-10);
        assert!((g[1][1] - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_volume_element_computation() {
        let config = test_config();
        let e1 = MultiVector::basis_vector(config.clone(), 0);
        let e2 = MultiVector::basis_vector(config.clone(), 1);
        let e3 = MultiVector::basis_vector(config.clone(), 2);

        let frame = vec![e1, e2, e3];
        let volume = compute_volume_element(&frame).unwrap();

        // 体积元素应该是 e₁∧e₂∧e₃，索引为 7 (0b111)
        assert!((volume.get_coefficient(7) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_wedge_except() {
        let config = test_config();
        let e1 = MultiVector::basis_vector(config.clone(), 0);
        let e2 = MultiVector::basis_vector(config.clone(), 1);
        let e3 = MultiVector::basis_vector(config.clone(), 2);

        let frame = vec![e1.clone(), e2.clone(), e3.clone()];

        // 省略第 0 个：e₂∧e₃（索引 6 = 0b110）
        let wedge0 = compute_wedge_except(&frame, 0).unwrap();
        assert!((wedge0.get_coefficient(6) - 1.0).abs() < 1e-10);

        // 省略第 1 个：e₁∧e₃（索引 5 = 0b101）
        let wedge1 = compute_wedge_except(&frame, 1).unwrap();
        // 注意：e₁∧e₃ 的索引是 5，但需要检查符号
        // e₁ ∧ e₃: 由于 e₁ 的索引是 0，e₃ 的索引是 2
        // 外积 e₁∧e₃ 没有交换，所以系数应该是 +1
        assert!((wedge1.get_coefficient(5) - 1.0).abs() < 1e-10);

        // 省略第 2 个：e₁∧e₂（索引 3 = 0b011）
        let wedge2 = compute_wedge_except(&frame, 2).unwrap();
        assert!((wedge2.get_coefficient(3) - 1.0).abs() < 1e-10);
    }

    // ========================================================================
    // 基展开测试 / Basis Expansion Tests
    // ========================================================================

    #[test]
    fn test_basis_expansion_scalar() {
        let config = test_config();
        let scalar = MultiVector::from_scalar(config.clone(), 5.0);
        let coeffs = basis_expansion(&scalar);

        assert_eq!(coeffs.len(), 1);
        assert_eq!(coeffs[0].0, 0); // 标量索引
        assert!((coeffs[0].1 - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_basis_expansion_vector() {
        let config = test_config();
        let e1 = MultiVector::basis_vector(config.clone(), 0);
        let coeffs = basis_expansion(&e1);

        assert_eq!(coeffs.len(), 1);
        assert_eq!(coeffs[0].0, 1); // e₁ 索引
        assert!((coeffs[0].1 - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_basis_expansion_bivector() {
        let config = test_config();
        let e1 = MultiVector::basis_vector(config.clone(), 0);
        let e2 = MultiVector::basis_vector(config.clone(), 1);
        let e12 = e1.outer_product(&e2);

        let coeffs = basis_expansion(&e12);

        assert_eq!(coeffs.len(), 1);
        assert_eq!(coeffs[0].0, 3); // e₁∧e₂ 索引 (0b011 = 3)
        assert!((coeffs[0].1 - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_basis_expansion_trivector() {
        let config = test_config();
        let e1 = MultiVector::basis_vector(config.clone(), 0);
        let e2 = MultiVector::basis_vector(config.clone(), 1);
        let e3 = MultiVector::basis_vector(config.clone(), 2);
        let e123 = e1.outer_product(&e2).outer_product(&e3);

        let coeffs = basis_expansion(&e123);

        assert_eq!(coeffs.len(), 1);
        assert_eq!(coeffs[0].0, 7); // e₁∧e₂∧e₃ 索引 (0b111 = 7)
        assert!((coeffs[0].1 - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_basis_expansion_mixed() {
        let config = test_config();

        // 创建 A = 2 + 3e₁ + 4e₁₂
        let scalar = MultiVector::from_scalar(config.clone(), 2.0);
        let e1 = MultiVector::basis_vector(config.clone(), 0);
        let e2 = MultiVector::basis_vector(config.clone(), 1);
        let e12 = e1.outer_product(&e2);

        let a = scalar.add(&e1.scale(3.0)).add(&e12.scale(4.0));
        let coeffs = basis_expansion(&a);

        assert_eq!(coeffs.len(), 3);

        // 验证各系数
        let get_coef = |idx: BasisIndex| -> f64 {
            coeffs
                .iter()
                .find(|(i, _)| *i == idx)
                .map(|(_, c)| *c)
                .unwrap_or(0.0)
        };

        assert!((get_coef(0) - 2.0).abs() < 1e-10); // 标量
        assert!((get_coef(1) - 3.0).abs() < 1e-10); // e₁
        assert!((get_coef(3) - 4.0).abs() < 1e-10); // e₁₂
    }

    #[test]
    fn test_basis_reconstruction() {
        let config = test_config();

        // 原始多向量
        let e1 = MultiVector::basis_vector(config.clone(), 0);
        let e2 = MultiVector::basis_vector(config.clone(), 1);
        let e3 = MultiVector::basis_vector(config.clone(), 2);
        let original = e1.scale(2.0).add(&e2.scale(3.0)).add(&e3.scale(4.0));

        // 展开并重构
        let coeffs = basis_expansion(&original);
        let reconstructed = basis_reconstruction(config, &coeffs);

        // 验证重构的多向量与原始相等
        let diff = original.sub(&reconstructed);
        assert!(diff.is_zero(), "Reconstruction should match original");
    }

    #[test]
    fn test_basis_expansion_roundtrip() {
        let config = test_config();

        // 创建复杂多向量
        let mut mv = MultiVector::from_scalar(config.clone(), 1.0);
        let e1 = MultiVector::basis_vector(config.clone(), 0);
        let e2 = MultiVector::basis_vector(config.clone(), 1);
        let e3 = MultiVector::basis_vector(config.clone(), 2);

        mv = mv.add(&e1.scale(2.0));
        mv = mv.add(&e2.scale(3.0));
        mv = mv.add(&e3.scale(4.0));
        mv = mv.add(&e1.outer_product(&e2).scale(5.0));
        mv = mv.add(&e2.outer_product(&e3).scale(6.0));
        mv = mv.add(&e1.outer_product(&e2).outer_product(&e3).scale(7.0));

        // 展开 -> 重构 -> 比较
        let coeffs = basis_expansion(&mv);
        let reconstructed = basis_reconstruction(config, &coeffs);
        let diff = mv.sub(&reconstructed);

        assert!(diff.is_zero(), "Roundtrip should preserve multivector");
    }

    #[test]
    fn test_create_basis_blade() {
        let config = test_config();

        // 标量
        let scalar = create_basis_blade(config.clone(), 0);
        assert!((scalar.scalar_part() - 1.0).abs() < 1e-10);

        // e₁
        let e1 = create_basis_blade(config.clone(), 1);
        assert!((e1.get_coefficient(1) - 1.0).abs() < 1e-10);

        // e₁∧e₂ (索引 3 = 0b011)
        let e12 = create_basis_blade(config.clone(), 3);
        assert!((e12.get_coefficient(3) - 1.0).abs() < 1e-10);

        // e₁∧e₂∧e₃ (索引 7 = 0b111)
        let e123 = create_basis_blade(config.clone(), 7);
        assert!((e123.get_coefficient(7) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_basis_expansion_with_reciprocal_orthonormal() {
        let config = test_config();

        // 对于正交归一基，互逆标架就是原标架
        let e1 = MultiVector::basis_vector(config.clone(), 0);
        let e2 = MultiVector::basis_vector(config.clone(), 1);
        let e3 = MultiVector::basis_vector(config.clone(), 2);

        // 创建所有基刃（简化：只测试向量部分）
        let basis_count = config.basis_count();
        let mut reciprocal_blades = Vec::with_capacity(basis_count);
        for i in 0..basis_count {
            reciprocal_blades.push(create_basis_blade(config.clone(), i as BasisIndex));
        }

        // 测试向量 e1 的展开
        let coeffs = basis_expansion_with_reciprocal(&e1, &reciprocal_blades);
        assert_eq!(coeffs.len(), 1);
        assert_eq!(coeffs[0].0, 1);
        assert!((coeffs[0].1 - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_basis_expansion_zero() {
        let config = test_config();
        let zero = MultiVector::zeros(config);
        let coeffs = basis_expansion(&zero);
        assert!(
            coeffs.is_empty(),
            "Zero multivector should have no coefficients"
        );
    }
}
