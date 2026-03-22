//! 任意维度测试套件
//!
//! 验证几何代数库对任意维度 n 的正确性
//!
//! 测试覆盖：
//! - 1D 到 10D 的欧几里得空间
//! - 不同签名 (p, q, r)
//! - 边界情况

use nblade::{AlgebraConfig, MultiVector, Signature};
use std::sync::Arc;

/// 测试 1D 到 10D 的欧几里得空间
#[test]
fn test_dimensions_1_to_10() {
    for dim in 1..=10 {
        let config = Arc::new(AlgebraConfig::euclidean(dim));

        // 验证基向量数量
        assert_eq!(
            config.basis_count(),
            1 << dim,
            "维度 {} 的基向量数量错误",
            dim
        );

        // 创建所有基向量
        let basis: Vec<MultiVector> = (0..dim)
            .map(|i| MultiVector::basis_vector(config.clone(), i))
            .collect();

        // 验证基向量的平方都为 1
        for (i, e_i) in basis.iter().enumerate() {
            let e_i_sq = e_i.geometric_product(e_i);
            assert!(
                (e_i_sq.get_coefficient(0) - 1.0).abs() < 1e-10,
                "维度 {} 的基向量 e{} 平方不为 1",
                dim,
                i
            );
        }

        // 验证不同基向量正交
        for i in 0..dim as usize {
            for j in (i + 1)..dim as usize {
                let inner = basis[i].left_inner(&basis[j]);
                assert!(
                    inner.is_zero(),
                    "维度 {} 的基向量 e{} 和 e{} 不正交",
                    dim,
                    i,
                    j
                );
            }
        }
    }
}

/// 测试不同签名
#[test]
fn test_different_signatures() {
    // G(2, 0, 0) - 2D 欧几里得
    let config = Arc::new(AlgebraConfig::new(2, Signature::new(2, 0, 0)));
    let e1 = MultiVector::basis_vector(config.clone(), 0);
    assert!((e1.geometric_product(&e1).get_coefficient(0) - 1.0).abs() < 1e-10);

    // G(1, 1, 0) - 2D 闵可夫斯基
    let config = Arc::new(AlgebraConfig::new(2, Signature::new(1, 1, 0)));
    let e0 = MultiVector::basis_vector(config.clone(), 0);
    let e1 = MultiVector::basis_vector(config.clone(), 1);
    assert!((e0.geometric_product(&e0).get_coefficient(0) - 1.0).abs() < 1e-10);
    assert!((e1.geometric_product(&e1).get_coefficient(0) + 1.0).abs() < 1e-10);

    // G(3, 0, 0) - 3D 欧几里得
    let config = Arc::new(AlgebraConfig::new(3, Signature::new(3, 0, 0)));
    for i in 0..3 {
        let e_i = MultiVector::basis_vector(config.clone(), i);
        assert!((e_i.geometric_product(&e_i).get_coefficient(0) - 1.0).abs() < 1e-10);
    }

    // G(1, 3, 0) - 4D 时空
    let config = Arc::new(AlgebraConfig::new(4, Signature::new(1, 3, 0)));
    let e0 = MultiVector::basis_vector(config.clone(), 0);
    assert!((e0.geometric_product(&e0).get_coefficient(0) - 1.0).abs() < 1e-10);
    for i in 1..4 {
        let e_i = MultiVector::basis_vector(config.clone(), i);
        assert!((e_i.geometric_product(&e_i).get_coefficient(0) + 1.0).abs() < 1e-10);
    }

    // G(4, 1, 0) - 5D 共形几何代数
    let config = Arc::new(AlgebraConfig::new(5, Signature::new(4, 1, 0)));
    assert_eq!(config.basis_count(), 32);
}

/// 测试外积的结合律
#[test]
fn test_outer_product_associativity() {
    for dim in 3..=6 {
        let config = Arc::new(AlgebraConfig::euclidean(dim));
        let e1 = MultiVector::basis_vector(config.clone(), 0);
        let e2 = MultiVector::basis_vector(config.clone(), 1);
        let e3 = MultiVector::basis_vector(config.clone(), 2);

        // (e1 ∧ e2) ∧ e3 = e1 ∧ (e2 ∧ e3)
        let left = e1.outer_product(&e2).outer_product(&e3);
        let right = e1.outer_product(&e2.outer_product(&e3));

        for i in 0..config.basis_count() {
            assert!(
                (left.get_coefficient(i) - right.get_coefficient(i)).abs() < 1e-10,
                "维度 {} 外积结合律不成立",
                dim
            );
        }
    }
}

/// 测试几何积的分配律
#[test]
fn test_geometric_product_distributivity() {
    for dim in 2..=5 {
        let config = Arc::new(AlgebraConfig::euclidean(dim));

        let a = MultiVector::basis_vector(config.clone(), 0);
        let b = MultiVector::basis_vector(config.clone(), 1);
        let c = if dim > 2 {
            MultiVector::basis_vector(config.clone(), 2)
        } else {
            MultiVector::from_scalar(config.clone(), 1.0)
        };

        let b_plus_c = b.add(&c);

        // a(b + c) = ab + ac
        let left = a.geometric_product(&b_plus_c);
        let right = a.geometric_product(&b).add(&a.geometric_product(&c));

        for i in 0..config.basis_count() {
            assert!(
                (left.get_coefficient(i) - right.get_coefficient(i)).abs() < 1e-10,
                "维度 {} 几何积分配律不成立",
                dim
            );
        }
    }
}

/// 测试阶次投影
#[test]
fn test_grade_projection() {
    for dim in 3..=8 {
        let config = Arc::new(AlgebraConfig::euclidean(dim));

        // 创建混合阶次多向量
        let mut mv = MultiVector::from_scalar(config.clone(), 1.0); // 标量
        mv = mv.add(&MultiVector::basis_vector(config.clone(), 0)); // 向量
        if dim >= 2 {
            let e1 = MultiVector::basis_vector(config.clone(), 0);
            let e2 = MultiVector::basis_vector(config.clone(), 1);
            mv = mv.add(&e1.outer_product(&e2)); // 二重向量
        }

        // 验证各阶次部分
        let grade_0 = mv.grade_projection(0);
        assert!((grade_0.get_coefficient(0) - 1.0).abs() < 1e-10);

        let grade_1 = mv.grade_projection(1);
        assert!((grade_1.get_coefficient(1) - 1.0).abs() < 1e-10);

        if dim >= 2 {
            let grade_2 = mv.grade_projection(2);
            assert!((grade_2.get_coefficient(3) - 1.0).abs() < 1e-10);
        }
    }
}

/// 测试偶部和奇部
#[test]
fn test_even_odd_parts() {
    for dim in 2..=6 {
        let config = Arc::new(AlgebraConfig::euclidean(dim));

        let mut mv = MultiVector::from_scalar(config.clone(), 1.0);
        mv = mv.add(&MultiVector::basis_vector(config.clone(), 0));
        if dim >= 2 {
            let e1 = MultiVector::basis_vector(config.clone(), 0);
            let e2 = MultiVector::basis_vector(config.clone(), 1);
            mv = mv.add(&e1.outer_product(&e2));
        }

        let even = mv.even_part();
        let odd = mv.odd_part();
        let sum = even.add(&odd);

        // 验证 mv = even + odd
        for i in 0..config.basis_count() {
            assert!(
                (sum.get_coefficient(i) - mv.get_coefficient(i)).abs() < 1e-10,
                "维度 {} 的偶部 + 奇部不等于原向量",
                dim
            );
        }
    }
}

/// 测试反转性质
#[test]
fn test_reversion_properties() {
    for dim in 2..=6 {
        let config = Arc::new(AlgebraConfig::euclidean(dim));

        let e1 = MultiVector::basis_vector(config.clone(), 0);
        let e2 = MultiVector::basis_vector(config.clone(), 1);
        let e12 = e1.outer_product(&e2);

        // (AB)† = B†A†
        let ab = e1.geometric_product(&e2);
        let ab_rev = ab.reversion();
        let b_rev_a_rev = e2.reversion().geometric_product(&e1.reversion());

        for i in 0..config.basis_count() {
            assert!(
                (ab_rev.get_coefficient(i) - b_rev_a_rev.get_coefficient(i)).abs() < 1e-10,
                "维度 {} 的反转性质 (AB)† = B†A†不成立",
                dim
            );
        }
    }
}

/// 测试对偶性质
#[test]
fn test_dual_properties() {
    for dim in 3..=6 {
        let config = Arc::new(AlgebraConfig::euclidean(dim));

        let e1 = MultiVector::basis_vector(config.clone(), 0);

        // 双重对偶：A⊥⊥ = ±A
        let dual = e1.dual();
        let dual_dual = dual.dual();

        // 在欧几里得空间中，双重对偶应该是 ±A
        // 具体符号取决于维度
        let norm_orig = e1.norm_squared();
        let norm_dual_dual = dual_dual.norm_squared();

        assert!(
            (norm_orig.abs() - norm_dual_dual.abs()).abs() < 1e-10,
            "维度 {} 的双重对偶范数不相等",
            dim
        );
    }
}

/// 测试逆的性质
#[test]
fn test_inverse_properties() {
    for dim in 2..=5 {
        let config = Arc::new(AlgebraConfig::euclidean(dim));

        // 测试向量的逆
        let e1 = MultiVector::basis_vector(config.clone(), 0);
        let e1_inv = e1.inverse().expect("向量应该可逆");

        // e1 * e1⁻¹ = 1
        let product = e1.geometric_product(&e1_inv);
        assert!(
            (product.get_coefficient(0) - 1.0).abs() < 1e-10,
            "维度 {} 的向量逆性质不成立",
            dim
        );

        // 测试二重向量的逆
        if dim >= 2 {
            let e2 = MultiVector::basis_vector(config.clone(), 1);
            let e12 = e1.outer_product(&e2);
            let e12_inv = e12.inverse().expect("二重向量应该可逆");

            // e12 * e12⁻¹ = 1
            let product = e12.geometric_product(&e12_inv);
            assert!(
                (product.get_coefficient(0) - 1.0).abs() < 1e-10,
                "维度 {} 的二重向量逆性质不成立",
                dim
            );
        }
    }
}

/// 测试旋转
#[test]
fn test_rotation() {
    use nblade::geometry::rotation::create_rotor;

    for dim in 3..=5 {
        let config = Arc::new(AlgebraConfig::euclidean(dim));

        let e1 = MultiVector::basis_vector(config.clone(), 0);
        let e2 = MultiVector::basis_vector(config.clone(), 1);
        let plane = e1.outer_product(&e2);

        // 创建 90 度转子
        let rotor = create_rotor(&plane, std::f64::consts::PI / 2.0);

        // 旋转 e1 应该得到 e2
        let rotated = e1.rotate_by(&rotor).expect("旋转应该成功");

        assert!(
            (rotated.get_coefficient(2) - 1.0).abs() < 1e-10,
            "维度 {} 的 90 度旋转不正确",
            dim
        );
    }
}

/// 测试反射
#[test]
fn test_reflection() {
    for dim in 2..=5 {
        let config = Arc::new(AlgebraConfig::euclidean(dim));

        let e1 = MultiVector::basis_vector(config.clone(), 0);
        let e2 = MultiVector::basis_vector(config.clone(), 1);

        // e2 在 e1 中反射应该得到 -e2
        let reflected = e2.reflect_in(&e1).expect("反射应该成功");

        assert!(
            (reflected.get_coefficient(2) + 1.0).abs() < 1e-10,
            "维度 {} 的反射不正确",
            dim
        );

        // 反射两次应该回到原向量
        let reflected_twice = reflected.reflect_in(&e1).expect("第二次反射应该成功");
        assert!(
            (reflected_twice.get_coefficient(2) - 1.0).abs() < 1e-10,
            "维度 {} 的双重反射不正确",
            dim
        );
    }
}

/// 测试投影
#[test]
fn test_projection() {
    for dim in 3..=5 {
        let config = Arc::new(AlgebraConfig::euclidean(dim));

        let e1 = MultiVector::basis_vector(config.clone(), 0);
        let e2 = MultiVector::basis_vector(config.clone(), 1);
        let plane = e1.outer_product(&e2);

        // e1 投影到 e1∧e2 平面应该还是 e1
        let proj = e1.project_to(&plane).expect("投影应该成功");
        assert!(
            (proj.get_coefficient(1) - 1.0).abs() < 1e-10,
            "维度 {} 的投影不正确",
            dim
        );

        // e3 (如果存在) 投影到 e1∧e2 平面应该为 0
        if dim >= 3 {
            let e3 = MultiVector::basis_vector(config.clone(), 2);
            let proj = e3.project_to(&plane).expect("投影应该成功");
            assert!(proj.is_zero(), "维度 {} 的正交向量投影不为零", dim);
        }
    }
}

/// 测试标量积的性质
#[test]
fn test_scalar_product_properties() {
    for dim in 2..=6 {
        let config = Arc::new(AlgebraConfig::euclidean(dim));

        let a = MultiVector::from_scalar(config.clone(), 1.0)
            .add(&MultiVector::basis_vector(config.clone(), 0));
        let b = MultiVector::from_scalar(config.clone(), 2.0)
            .add(&MultiVector::basis_vector(config.clone(), 1));

        // 对称性：A * B = B * A
        let ab = a.scalar_product(&b);
        let ba = b.scalar_product(&a);
        assert!((ab - ba).abs() < 1e-10, "维度 {} 的标量积不对称", dim);

        // 对合不变：A * B = A* * B*
        let a_star = a.grade_involution();
        let b_star = b.grade_involution();
        let ab_star = a_star.scalar_product(&b_star);
        assert!(
            (ab - ab_star).abs() < 1e-10,
            "维度 {} 的标量积对合不变性不成立",
            dim
        );
    }
}

/// 测试范数不变性
#[test]
fn test_norm_invariance() {
    for dim in 2..=6 {
        let config = Arc::new(AlgebraConfig::euclidean(dim));

        let mut mv = MultiVector::from_scalar(config.clone(), 1.0);
        for i in 0..dim.min(4) {
            mv = mv.add(&MultiVector::basis_vector(config.clone(), i));
        }
        if dim >= 2 {
            let e1 = MultiVector::basis_vector(config.clone(), 0);
            let e2 = MultiVector::basis_vector(config.clone(), 1);
            mv = mv.add(&e1.outer_product(&e2));
        }

        let norm_sq = mv.norm_squared();
        let norm_sq_star = mv.grade_involution().norm_squared();
        let norm_sq_dagger = mv.reversion().norm_squared();
        let norm_sq_conj = mv.clifford_conjugate().norm_squared();

        assert!(
            (norm_sq - norm_sq_star).abs() < 1e-10,
            "维度 {} 的阶次对合范数不变性不成立",
            dim
        );
        assert!(
            (norm_sq - norm_sq_dagger).abs() < 1e-10,
            "维度 {} 的反转范数不变性不成立",
            dim
        );
        assert!(
            (norm_sq - norm_sq_conj).abs() < 1e-10,
            "维度 {} 的 Clifford 共轭范数不变性不成立",
            dim
        );
    }
}

/// 测试交换子
#[test]
fn test_commutator() {
    for dim in 2..=5 {
        let config = Arc::new(AlgebraConfig::euclidean(dim));

        let e1 = MultiVector::basis_vector(config.clone(), 0);
        let e2 = MultiVector::basis_vector(config.clone(), 1);

        // 反对称性：A × B = -B × A
        let ab = e1.commutator(&e2);
        let ba = e2.commutator(&e1);

        for i in 0..config.basis_count() {
            assert!(
                (ab.get_coefficient(i) + ba.get_coefficient(i)).abs() < 1e-10,
                "维度 {} 的交换子不对称",
                dim
            );
        }

        // A × A = 0
        let aa = e1.commutator(&e1);
        assert!(aa.is_zero(), "维度 {} 的自交换子不为零", dim);
    }
}

/// 测试高维空间 (n > 10)
#[test]
fn test_high_dimensions() {
    // 测试 10D, 12D, 16D
    for dim in [10, 12, 16] {
        let config = Arc::new(AlgebraConfig::euclidean(dim));

        // 验证基向量数量
        assert_eq!(config.basis_count(), 1 << dim);

        // 创建基向量并验证平方
        let e1 = MultiVector::basis_vector(config.clone(), 0);
        let e2 = MultiVector::basis_vector(config.clone(), 1);

        assert!((e1.geometric_product(&e1).get_coefficient(0) - 1.0).abs() < 1e-10);
        assert!(e1.left_inner(&e2).is_zero());

        // 外积
        let e12 = e1.outer_product(&e2);
        assert!((e12.norm_squared() - 1.0).abs() < 1e-10);
    }
}

/// 测试基展开 / Test Basis Expansion
#[test]
fn test_basis_expansion_integration() {
    use nblade::basis::{basis_expansion, basis_reconstruction};

    // 测试 2D 到 5D
    for dim in 2..=5 {
        let config = Arc::new(AlgebraConfig::euclidean(dim));

        // 创建复杂多向量
        let mut mv = MultiVector::from_scalar(config.clone(), 1.0);
        for i in 0..dim.min(3) {
            let e = MultiVector::basis_vector(config.clone(), i);
            mv = mv.add(&e.scale((i + 1) as f64));
        }

        // 添加一个 bivector
        if dim >= 2 {
            let e1 = MultiVector::basis_vector(config.clone(), 0);
            let e2 = MultiVector::basis_vector(config.clone(), 1);
            let e12 = e1.outer_product(&e2);
            mv = mv.add(&e12.scale(5.0));
        }

        // 基展开
        let coeffs = basis_expansion(&mv);

        // 重构
        let reconstructed = basis_reconstruction(config.clone(), &coeffs);

        // 验证重构的多向量与原始相等
        let diff = mv.sub(&reconstructed);
        assert!(diff.is_zero(), "维度 {} 的基展开-重构循环失败", dim);
    }
}

/// 测试基展开与标量积的关系 / Test Basis Expansion Scalar Product Relationship
#[test]
fn test_basis_expansion_scalar_product() {
    let config = Arc::new(AlgebraConfig::euclidean(3));

    // 创建 A = 2e₁ + 3e₂
    let e1 = MultiVector::basis_vector(config.clone(), 0);
    let e2 = MultiVector::basis_vector(config.clone(), 1);
    let a = e1.scale(2.0).add(&e2.scale(3.0));

    // 验证 A * e₁ = 2（标量积给出 e₁ 的系数）
    let sp_e1 = a.scalar_product(&e1);
    assert!((sp_e1 - 2.0).abs() < 1e-10, "A * e₁ 应为 2");

    // 验证 A * e₂ = 3
    let sp_e2 = a.scalar_product(&e2);
    assert!((sp_e2 - 3.0).abs() < 1e-10, "A * e₂ 应为 3");

    // 验证 A * e₃ = 0
    let e3 = MultiVector::basis_vector(config.clone(), 2);
    let sp_e3 = a.scalar_product(&e3);
    assert!(sp_e3.abs() < 1e-10, "A * e₃ 应为 0");
}
