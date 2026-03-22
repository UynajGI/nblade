//! Property-Based Tests for Geometric Algebra Operations
//!
//! Verifies formulas from summary.md using proptest for random test generation.

use nblade::{AlgebraConfig, MultiVector, Signature};
use proptest::prelude::*;
use std::sync::Arc;

fn test_config(dim: u32) -> Arc<AlgebraConfig> {
    Arc::new(AlgebraConfig::euclidean(dim))
}

// ============================================================================
// GEOMETRIC PRODUCT PROPERTIES
// ============================================================================

mod geometric_product_tests {
    use super::*;

    #[test]
    fn test_geometric_product_associativity() {
        proptest!(|(dim in 2u32..=5, a_idx in 0u64..8, a_c in -10.0f64..10.0f64,
                    b_idx in 0u64..8, b_c in -10.0f64..10.0f64,
                    c_idx in 0u64..8, c_c in -10.0f64..10.0f64)| {
            let config = test_config(dim);
            let basis_count = config.basis_count();
            let a_idx = a_idx.min(basis_count as u64 - 1);
            let b_idx = b_idx.min(basis_count as u64 - 1);
            let c_idx = c_idx.min(basis_count as u64 - 1);

            let a = {
                let mut coeffs = vec![0.0; basis_count];
                coeffs[a_idx as usize] = a_c;
                MultiVector::from_coefficients(config.clone(), coeffs)
            };
            let b = {
                let mut coeffs = vec![0.0; basis_count];
                coeffs[b_idx as usize] = b_c;
                MultiVector::from_coefficients(config.clone(), coeffs)
            };
            let c = {
                let mut coeffs = vec![0.0; basis_count];
                coeffs[c_idx as usize] = c_c;
                MultiVector::from_coefficients(config.clone(), coeffs)
            };

            let ab_c = a.geometric_product(&b).geometric_product(&c);
            let a_bc = a.geometric_product(&b.geometric_product(&c));

            for i in 0..basis_count {
                let diff = (ab_c.get_coefficient(i) - a_bc.get_coefficient(i)).abs();
                prop_assert!(diff < 1e-8);
            }
        });
    }

    #[test]
    fn test_geometric_product_distributivity_left() {
        proptest!(|(dim in 2u32..=5, a_idx in 0u64..8, a_c in -5.0f64..5.0f64,
                    b_idx in 0u64..8, b_c in -5.0f64..5.0f64,
                    c_idx in 0u64..8, c_c in -5.0f64..5.0f64)| {
            let config = test_config(dim);
            let basis_count = config.basis_count();
            let a_idx = a_idx.min(basis_count as u64 - 1);
            let b_idx = b_idx.min(basis_count as u64 - 1);
            let c_idx = c_idx.min(basis_count as u64 - 1);

            let a = {
                let mut coeffs = vec![0.0; basis_count];
                coeffs[a_idx as usize] = a_c;
                MultiVector::from_coefficients(config.clone(), coeffs)
            };
            let b = {
                let mut coeffs = vec![0.0; basis_count];
                coeffs[b_idx as usize] = b_c;
                MultiVector::from_coefficients(config.clone(), coeffs)
            };
            let c = {
                let mut coeffs = vec![0.0; basis_count];
                coeffs[c_idx as usize] = c_c;
                MultiVector::from_coefficients(config.clone(), coeffs)
            };

            let left = a.geometric_product(&b.add(&c));
            let right = a.geometric_product(&b).add(&a.geometric_product(&c));

            for i in 0..basis_count {
                let diff = (left.get_coefficient(i) - right.get_coefficient(i)).abs();
                prop_assert!(diff < 1e-8);
            }
        });
    }

    #[test]
    fn test_geometric_product_distributivity_right() {
        proptest!(|(dim in 2u32..=5, a_idx in 0u64..8, a_c in -5.0f64..5.0f64,
                    b_idx in 0u64..8, b_c in -5.0f64..5.0f64,
                    c_idx in 0u64..8, c_c in -5.0f64..5.0f64)| {
            let config = test_config(dim);
            let basis_count = config.basis_count();
            let a_idx = a_idx.min(basis_count as u64 - 1);
            let b_idx = b_idx.min(basis_count as u64 - 1);
            let c_idx = c_idx.min(basis_count as u64 - 1);

            let a = {
                let mut coeffs = vec![0.0; basis_count];
                coeffs[a_idx as usize] = a_c;
                MultiVector::from_coefficients(config.clone(), coeffs)
            };
            let b = {
                let mut coeffs = vec![0.0; basis_count];
                coeffs[b_idx as usize] = b_c;
                MultiVector::from_coefficients(config.clone(), coeffs)
            };
            let c = {
                let mut coeffs = vec![0.0; basis_count];
                coeffs[c_idx as usize] = c_c;
                MultiVector::from_coefficients(config.clone(), coeffs)
            };

            let left = a.add(&b).geometric_product(&c);
            let right = a.geometric_product(&c).add(&b.geometric_product(&c));

            for i in 0..basis_count {
                let diff = (left.get_coefficient(i) - right.get_coefficient(i)).abs();
                prop_assert!(diff < 1e-8);
            }
        });
    }

    #[test]
    fn test_scalar_commutativity() {
        proptest!(|(dim in 2u32..=5, s in -10.0f64..10.0f64,
                    mv_idx in 0u64..8, mv_c in -5.0f64..5.0f64)| {
            let config = test_config(dim);
            let basis_count = config.basis_count();
            let mv_idx = mv_idx.min(basis_count as u64 - 1);

            let scalar = MultiVector::from_scalar(config.clone(), s);
            let mv = {
                let mut coeffs = vec![0.0; basis_count];
                coeffs[mv_idx as usize] = mv_c;
                MultiVector::from_coefficients(config.clone(), coeffs)
            };

            let left = scalar.geometric_product(&mv);
            let right = mv.geometric_product(&scalar);

            for i in 0..basis_count {
                let diff = (left.get_coefficient(i) - right.get_coefficient(i)).abs();
                prop_assert!(diff < 1e-8);
            }
        });
    }

    #[test]
    fn test_scalar_multiplication_distributivity() {
        proptest!(|(dim in 2u32..=5, s in -10.0f64..10.0f64,
                    a_idx in 0u64..8, a_c in -5.0f64..5.0f64,
                    b_idx in 0u64..8, b_c in -5.0f64..5.0f64)| {
            let config = test_config(dim);
            let basis_count = config.basis_count();
            let a_idx = a_idx.min(basis_count as u64 - 1);
            let b_idx = b_idx.min(basis_count as u64 - 1);

            let scalar = MultiVector::from_scalar(config.clone(), s);
            let a = {
                let mut coeffs = vec![0.0; basis_count];
                coeffs[a_idx as usize] = a_c;
                MultiVector::from_coefficients(config.clone(), coeffs)
            };
            let b = {
                let mut coeffs = vec![0.0; basis_count];
                coeffs[b_idx as usize] = b_c;
                MultiVector::from_coefficients(config.clone(), coeffs)
            };

            let left = scalar.geometric_product(&a.add(&b));
            let right = scalar.geometric_product(&a).add(&scalar.geometric_product(&b));

            for i in 0..basis_count {
                let diff = (left.get_coefficient(i) - right.get_coefficient(i)).abs();
                prop_assert!(diff < 1e-8);
            }
        });
    }
}

// ============================================================================
// OUTER PRODUCT PROPERTIES
// ============================================================================

mod outer_product_tests {
    use super::*;

    #[test]
    fn test_outer_product_antisymmetry_vectors() {
        proptest!(|(dim in 2u32..=5, i in 0u32..4, j in 0u32..4,
                    a_c in -5.0f64..5.0f64, b_c in -5.0f64..5.0f64)| {
            let config = test_config(dim);
            let i = i.min(dim - 1);
            let j = j.min(dim - 1);

            let a = MultiVector::basis_vector(config.clone(), i).scale(a_c);
            let b = MultiVector::basis_vector(config.clone(), j).scale(b_c);

            let ab = a.outer_product(&b);
            let neg_ba = b.outer_product(&a).neg();

            for k in 0..config.basis_count() {
                let diff = (ab.get_coefficient(k) - neg_ba.get_coefficient(k)).abs();
                prop_assert!(diff < 1e-8);
            }
        });
    }

    #[test]
    fn test_outer_product_associativity() {
        proptest!(|(dim in 3u32..=5, i in 0u32..3, j in 0u32..3, k in 0u32..3)| {
            let config = test_config(dim);
            let i = i.min(dim - 1);
            let j = j.min(dim - 1);
            let k = k.min(dim - 1);

            let a = MultiVector::basis_vector(config.clone(), i);
            let b = MultiVector::basis_vector(config.clone(), j);
            let c = MultiVector::basis_vector(config.clone(), k);

            let left = a.outer_product(&b).outer_product(&c);
            let right = a.outer_product(&b.outer_product(&c));

            for idx in 0..config.basis_count() {
                let diff = (left.get_coefficient(idx) - right.get_coefficient(idx)).abs();
                prop_assert!(diff < 1e-8);
            }
        });
    }

    #[test]
    fn test_outer_product_self_wedge_zero() {
        proptest!(|(dim in 2u32..=5, i in 0u32..4, c in -10.0f64..10.0f64)| {
            let config = test_config(dim);
            let i = i.min(dim - 1);

            let a = MultiVector::basis_vector(config.clone(), i).scale(c);
            let result = a.outer_product(&a);

            prop_assert!(result.is_zero());
        });
    }

    #[test]
    fn test_outer_product_scalar() {
        proptest!(|(dim in 2u32..=5, s in -5.0f64..5.0f64,
                    i in 0u32..4, a_c in -5.0f64..5.0f64)| {
            let config = test_config(dim);
            let i = i.min(dim - 1);

            let scalar = MultiVector::from_scalar(config.clone(), s);
            let a = MultiVector::basis_vector(config.clone(), i).scale(a_c);

            let left = scalar.outer_product(&a);
            let right = a.outer_product(&scalar);
            let expected = a.scale(s);

            for k in 0..config.basis_count() {
                let diff1 = (left.get_coefficient(k) - expected.get_coefficient(k)).abs();
                let diff2 = (right.get_coefficient(k) - expected.get_coefficient(k)).abs();
                prop_assert!(diff1 < 1e-8);
                prop_assert!(diff2 < 1e-8);
            }
        });
    }

    #[test]
    fn test_outer_product_distributivity() {
        proptest!(|(dim in 3u32..=5, i in 0u32..3, a_c in -5.0f64..5.0f64,
                    j in 0u32..3, b_c in -5.0f64..5.0f64,
                    k in 0u32..3, c_c in -5.0f64..5.0f64)| {
            let config = test_config(dim);
            let i = i.min(dim - 1);
            let j = j.min(dim - 1);
            let k = k.min(dim - 1);

            let a = MultiVector::basis_vector(config.clone(), i).scale(a_c);
            let b = MultiVector::basis_vector(config.clone(), j).scale(b_c);
            let c = MultiVector::basis_vector(config.clone(), k).scale(c_c);

            let left = a.outer_product(&b.add(&c));
            let right = a.outer_product(&b).add(&a.outer_product(&c));

            for idx in 0..config.basis_count() {
                let diff = (left.get_coefficient(idx) - right.get_coefficient(idx)).abs();
                prop_assert!(diff < 1e-8);
            }
        });
    }

    #[test]
    fn test_outer_product_grade_sign() {
        proptest!(|(dim in 4u32..=5)| {
            let config = test_config(dim);

            let e1 = MultiVector::basis_vector(config.clone(), 0);
            let e2 = MultiVector::basis_vector(config.clone(), 1);
            let e3 = MultiVector::basis_vector(config.clone(), 2);
            let e12 = e1.outer_product(&e2);

            let left = e12.outer_product(&e3);
            let right = e3.outer_product(&e12);

            for idx in 0..config.basis_count() {
                let diff = (left.get_coefficient(idx) - right.get_coefficient(idx)).abs();
                prop_assert!(diff < 1e-8);
            }
        });
    }
}

// ============================================================================
// INNER PRODUCT PROPERTIES
// ============================================================================

mod inner_product_tests {
    use super::*;

    #[test]
    fn test_inner_product_orthogonality() {
        proptest!(|(dim in 2u32..=5, i in 0u32..4, j in 0u32..4)| {
            let config = test_config(dim);
            let i = i.min(dim - 1);
            let j = j.min(dim - 1);

            let e_i = MultiVector::basis_vector(config.clone(), i);
            let e_j = MultiVector::basis_vector(config.clone(), j);

            let result = e_i.left_inner(&e_j);

            if i != j {
                prop_assert!(result.is_zero());
            } else {
                prop_assert!((result.scalar_part() - 1.0).abs() < 1e-8);
            }
        });
    }

    #[test]
    fn test_inner_product_scalar() {
        proptest!(|(dim in 2u32..=5, s in -5.0f64..5.0f64,
                    i in 0u32..4, a_c in -5.0f64..5.0f64)| {
            let config = test_config(dim);
            let i = i.min(dim - 1);

            let scalar = MultiVector::from_scalar(config.clone(), s);
            let a = MultiVector::basis_vector(config.clone(), i).scale(a_c);

            // s⌋a = s*a (scalar contracts to give scaled vector)
            let result = scalar.left_inner(&a);
            let expected = a.scale(s);
            for k in 0..config.basis_count() {
                let diff = (result.get_coefficient(k) - expected.get_coefficient(k)).abs();
                prop_assert!(diff < 1e-8);
            }

            // a⌋s = 0 (vector cannot contract into scalar)
            let result2 = a.left_inner(&scalar);
            prop_assert!(result2.is_zero());
        });
    }

    #[test]
    fn test_left_inner_distributivity() {
        proptest!(|(dim in 3u32..=5, i in 0u32..3, a_c in -3.0f64..3.0f64,
                    j in 0u32..3, b_c in -3.0f64..3.0f64,
                    k in 0u32..3, c_c in -3.0f64..3.0f64)| {
            let config = test_config(dim);
            let i = i.min(dim - 1);
            let j = j.min(dim - 1);
            let k = k.min(dim - 1);

            let a = MultiVector::basis_vector(config.clone(), i).scale(a_c);
            let b = MultiVector::basis_vector(config.clone(), j).scale(b_c);
            let c = MultiVector::basis_vector(config.clone(), k).scale(c_c);

            let left = a.left_inner(&b.add(&c));
            let right = a.left_inner(&b).add(&a.left_inner(&c));

            for idx in 0..config.basis_count() {
                let diff = (left.get_coefficient(idx) - right.get_coefficient(idx)).abs();
                prop_assert!(diff < 1e-8);
            }
        });
    }

    #[test]
    fn test_right_inner_distributivity() {
        proptest!(|(dim in 3u32..=5, i in 0u32..3, a_c in -3.0f64..3.0f64,
                    j in 0u32..3, b_c in -3.0f64..3.0f64,
                    k in 0u32..3, c_c in -3.0f64..3.0f64)| {
            let config = test_config(dim);
            let i = i.min(dim - 1);
            let j = j.min(dim - 1);
            let k = k.min(dim - 1);

            let a = MultiVector::basis_vector(config.clone(), i).scale(a_c);
            let b = MultiVector::basis_vector(config.clone(), j).scale(b_c);
            let c = MultiVector::basis_vector(config.clone(), k).scale(c_c);

            let left = a.add(&b).right_inner(&c);
            let right = a.right_inner(&c).add(&b.right_inner(&c));

            for idx in 0..config.basis_count() {
                let diff = (left.get_coefficient(idx) - right.get_coefficient(idx)).abs();
                prop_assert!(diff < 1e-8);
            }
        });
    }

    #[test]
    fn test_inner_product_grade_projection() {
        proptest!(|(dim in 3u32..=5, i in 0u32..3)| {
            let config = test_config(dim);
            let i = i.min(dim - 1);

            let e = MultiVector::basis_vector(config.clone(), i);
            let result = e.left_inner(&e);

            let scalar_part = result.scalar_part();
            prop_assert!((scalar_part - 1.0).abs() < 1e-8);

            for idx in 1..config.basis_count() {
                let coef = result.get_coefficient(idx);
                prop_assert!(coef.abs() < 1e-10);
            }
        });
    }
}

// ============================================================================
// INVOLUTION PROPERTIES
// ============================================================================

mod involution_tests {
    use super::*;

    #[test]
    fn test_grade_involution_formula() {
        proptest!(|(dim in 2u32..=5, scalar_c in -5.0f64..5.0f64,
                    vec_c in -5.0f64..5.0f64, i in 0u32..4)| {
            let config = test_config(dim);
            let i = i.min(dim - 1);

            let scalar = MultiVector::from_scalar(config.clone(), scalar_c);
            let vec = MultiVector::basis_vector(config.clone(), i).scale(vec_c);
            let mv = scalar.add(&vec);

            let result = mv.grade_involution();

            prop_assert!((result.scalar_part() - scalar_c).abs() < 1e-8);

            let vec_coef = result.get_coefficient((1 << i) as usize);
            prop_assert!((vec_coef + vec_c).abs() < 1e-8);
        });
    }

    #[test]
    fn test_grade_involution_involution() {
        proptest!(|(dim in 2u32..=5, idx in 0u64..8, c in -5.0f64..5.0f64)| {
            let config = test_config(dim);
            let basis_count = config.basis_count();
            let idx = idx.min(basis_count as u64 - 1);

            let mut coeffs = vec![0.0; basis_count];
            coeffs[idx as usize] = c;
            let mv = MultiVector::from_coefficients(config.clone(), coeffs);

            let double_involution = mv.grade_involution().grade_involution();

            for k in 0..basis_count {
                let diff = (mv.get_coefficient(k) - double_involution.get_coefficient(k)).abs();
                prop_assert!(diff < 1e-8);
            }
        });
    }

    #[test]
    fn test_reversion_vector() {
        proptest!(|(dim in 2u32..=5, i in 0u32..4, c in -5.0f64..5.0f64)| {
            let config = test_config(dim);
            let i = i.min(dim - 1);

            let v = MultiVector::basis_vector(config.clone(), i).scale(c);
            let result = v.reversion();

            for k in 0..config.basis_count() {
                let diff = (v.get_coefficient(k) - result.get_coefficient(k)).abs();
                prop_assert!(diff < 1e-8);
            }
        });
    }

    #[test]
    fn test_reversion_bivector() {
        proptest!(|(dim in 3u32..=5, i in 0u32..3, j in 0u32..3, c in -3.0f64..3.0f64)| {
            let config = test_config(dim);
            let i = i.min(dim - 1);
            let j = j.min(dim - 1);

            if i != j {
                let e_i = MultiVector::basis_vector(config.clone(), i);
                let e_j = MultiVector::basis_vector(config.clone(), j);
                let bivector = e_i.outer_product(&e_j).scale(c);

                let result = bivector.reversion();

                for k in 0..config.basis_count() {
                    let expected = -bivector.get_coefficient(k);
                    let diff = (result.get_coefficient(k) - expected).abs();
                    prop_assert!(diff < 1e-8);
                }
            }
        });
    }

    #[test]
    fn test_reversion_involution() {
        proptest!(|(dim in 2u32..=5, idx in 0u64..8, c in -5.0f64..5.0f64)| {
            let config = test_config(dim);
            let basis_count = config.basis_count();
            let idx = idx.min(basis_count as u64 - 1);

            let mut coeffs = vec![0.0; basis_count];
            coeffs[idx as usize] = c;
            let mv = MultiVector::from_coefficients(config.clone(), coeffs);

            let double_reversion = mv.reversion().reversion();

            for k in 0..basis_count {
                let diff = (mv.get_coefficient(k) - double_reversion.get_coefficient(k)).abs();
                prop_assert!(diff < 1e-8);
            }
        });
    }

    #[test]
    fn test_clifford_conjugate_formula() {
        proptest!(|(dim in 2u32..=5, idx in 0u64..8, c in -5.0f64..5.0f64)| {
            let config = test_config(dim);
            let basis_count = config.basis_count();
            let idx = idx.min(basis_count as u64 - 1);

            let mut coeffs = vec![0.0; basis_count];
            coeffs[idx as usize] = c;
            let mv = MultiVector::from_coefficients(config.clone(), coeffs);

            let clifford = mv.clifford_conjugate();
            let involution_then_reversion = mv.grade_involution().reversion();
            let reversion_then_involution = mv.reversion().grade_involution();

            for k in 0..basis_count {
                let diff1 = (clifford.get_coefficient(k) - involution_then_reversion.get_coefficient(k)).abs();
                let diff2 = (clifford.get_coefficient(k) - reversion_then_involution.get_coefficient(k)).abs();
                prop_assert!(diff1 < 1e-8);
                prop_assert!(diff2 < 1e-8);
            }
        });
    }

    #[test]
    fn test_clifford_conjugate_involution() {
        proptest!(|(dim in 2u32..=5, idx in 0u64..8, c in -5.0f64..5.0f64)| {
            let config = test_config(dim);
            let basis_count = config.basis_count();
            let idx = idx.min(basis_count as u64 - 1);

            let mut coeffs = vec![0.0; basis_count];
            coeffs[idx as usize] = c;
            let mv = MultiVector::from_coefficients(config.clone(), coeffs);

            let double_clifford = mv.clifford_conjugate().clifford_conjugate();

            for k in 0..basis_count {
                let diff = (mv.get_coefficient(k) - double_clifford.get_coefficient(k)).abs();
                prop_assert!(diff < 1e-8);
            }
        });
    }

    #[test]
    fn test_product_reversion() {
        proptest!(|(dim in 2u32..=5, a_idx in 0u64..8, a_c in -3.0f64..3.0f64,
                    b_idx in 0u64..8, b_c in -3.0f64..3.0f64)| {
            let config = test_config(dim);
            let basis_count = config.basis_count();
            let a_idx = a_idx.min(basis_count as u64 - 1);
            let b_idx = b_idx.min(basis_count as u64 - 1);

            let a = {
                let mut coeffs = vec![0.0; basis_count];
                coeffs[a_idx as usize] = a_c;
                MultiVector::from_coefficients(config.clone(), coeffs)
            };
            let b = {
                let mut coeffs = vec![0.0; basis_count];
                coeffs[b_idx as usize] = b_c;
                MultiVector::from_coefficients(config.clone(), coeffs)
            };

            let ab_rev = a.geometric_product(&b).reversion();
            let b_rev_a_rev = b.reversion().geometric_product(&a.reversion());

            for k in 0..basis_count {
                let diff = (ab_rev.get_coefficient(k) - b_rev_a_rev.get_coefficient(k)).abs();
                prop_assert!(diff < 1e-8);
            }
        });
    }

    #[test]
    fn test_product_involution() {
        proptest!(|(dim in 2u32..=5, a_idx in 0u64..8, a_c in -3.0f64..3.0f64,
                    b_idx in 0u64..8, b_c in -3.0f64..3.0f64)| {
            let config = test_config(dim);
            let basis_count = config.basis_count();
            let a_idx = a_idx.min(basis_count as u64 - 1);
            let b_idx = b_idx.min(basis_count as u64 - 1);

            let a = {
                let mut coeffs = vec![0.0; basis_count];
                coeffs[a_idx as usize] = a_c;
                MultiVector::from_coefficients(config.clone(), coeffs)
            };
            let b = {
                let mut coeffs = vec![0.0; basis_count];
                coeffs[b_idx as usize] = b_c;
                MultiVector::from_coefficients(config.clone(), coeffs)
            };

            let ab_inv = a.geometric_product(&b).grade_involution();
            let a_inv_b_inv = a.grade_involution().geometric_product(&b.grade_involution());

            for k in 0..basis_count {
                let diff = (ab_inv.get_coefficient(k) - a_inv_b_inv.get_coefficient(k)).abs();
                prop_assert!(diff < 1e-8);
            }
        });
    }

    #[test]
    fn test_norm_invariance() {
        proptest!(|(dim in 2u32..=5, scalar_c in -3.0f64..3.0f64,
                    vec_c in -3.0f64..3.0f64, i in 0u32..4)| {
            let config = test_config(dim);
            let i = i.min(dim - 1);

            let scalar = MultiVector::from_scalar(config.clone(), scalar_c);
            let vec = MultiVector::basis_vector(config.clone(), i).scale(vec_c);
            let mv = scalar.add(&vec);

            let norm_orig = mv.norm_squared();
            let norm_star = mv.grade_involution().norm_squared();
            let norm_dagger = mv.reversion().norm_squared();
            let norm_clifford = mv.clifford_conjugate().norm_squared();

            prop_assert!((norm_orig - norm_star).abs() < 1e-8);
            prop_assert!((norm_orig - norm_dagger).abs() < 1e-8);
            prop_assert!((norm_orig - norm_clifford).abs() < 1e-8);
        });
    }
}

// ============================================================================
// DUAL PROPERTIES
// ============================================================================

mod dual_tests {
    use super::*;

    #[test]
    fn test_dual_inverse_dual_relationship() {
        proptest!(|(dim in 3u32..=5, idx in 0u64..8, c in -3.0f64..3.0f64)| {
            let config = test_config(dim);
            let basis_count = config.basis_count();
            let idx = idx.min(basis_count as u64 - 1);

            let mut coeffs = vec![0.0; basis_count];
            coeffs[idx as usize] = c;
            let mv = MultiVector::from_coefficients(config.clone(), coeffs);

            let dual = mv.dual();
            let inverse_dual = dual.inverse_dual();

            let n = dim as usize;
            let grade = config.grade_of_index(idx as usize) as usize;

            // Sign formula: inverse_dual(dual(A_r)) = (-1)^{n(n-1)/2 + r(r-1)/2} * A_r
            let i_squared_sign = (n * (n - 1) / 2) % 2;
            let reversion_sign = if grade == 0 { 0 } else { (grade * (grade - 1) / 2) % 2 };
            let total_sign_exp = i_squared_sign + reversion_sign;
            let expected_sign = if total_sign_exp % 2 == 0 { 1.0 } else { -1.0 };

            let expected = mv.scale(expected_sign);

            for k in 0..basis_count {
                let diff = (inverse_dual.get_coefficient(k) - expected.get_coefficient(k)).abs();
                prop_assert!(diff < 1e-8);
            }
        });
    }

    #[test]
    fn test_dual_scalar() {
        proptest!(|(dim in 3u32..=5, c in -5.0f64..5.0f64)| {
            let config = test_config(dim);
            let scalar = MultiVector::from_scalar(config.clone(), c);

            let dual = scalar.dual();

            let volume_idx = (1u64 << dim) - 1;
            let coef = dual.get_coefficient(volume_idx as usize);
            prop_assert!(coef.abs() > 0.0);
        });
    }

    #[test]
    fn test_dual_volume() {
        proptest!(|(dim in 3u32..=5)| {
            let config = test_config(dim);
            let volume_idx = (1u64 << dim) - 1;

            let mut coeffs = vec![0.0; config.basis_count()];
            coeffs[volume_idx as usize] = 1.0;
            let volume = MultiVector::from_coefficients(config.clone(), coeffs);

            let dual = volume.dual();

            let scalar_part = dual.scalar_part().abs();
            prop_assert!(scalar_part > 0.0);

            for k in 1..config.basis_count() {
                let coef = dual.get_coefficient(k);
                prop_assert!(coef.abs() < 1e-8);
            }
        });
    }

    #[test]
    fn test_dual_grade_change() {
        proptest!(|(dim in 3u32..=5, i in 0u32..4, c in -3.0f64..3.0f64)| {
            let config = test_config(dim);
            let i = i.min(dim - 1);

            let vec = MultiVector::basis_vector(config.clone(), i).scale(c);
            let dual = vec.dual();

            let n = dim as usize;

            let mut non_zero_grades = vec![];
            for k in 0..config.basis_count() {
                let coef = dual.get_coefficient(k);
                if coef.abs() > 1e-10 {
                    let grade = config.grade_of_index(k);
                    non_zero_grades.push(grade);
                }
            }

            for grade in non_zero_grades {
                prop_assert!(grade as usize == n - 1);
            }
        });
    }
}

// ============================================================================
// SCALAR PRODUCT PROPERTIES
// ============================================================================

mod scalar_product_tests {
    use super::*;

    #[test]
    fn test_scalar_product_symmetry() {
        proptest!(|(dim in 2u32..=5, a_idx in 0u64..8, a_c in -3.0f64..3.0f64,
                    b_idx in 0u64..8, b_c in -3.0f64..3.0f64)| {
            let config = test_config(dim);
            let basis_count = config.basis_count();
            let a_idx = a_idx.min(basis_count as u64 - 1);
            let b_idx = b_idx.min(basis_count as u64 - 1);

            let a = {
                let mut coeffs = vec![0.0; basis_count];
                coeffs[a_idx as usize] = a_c;
                MultiVector::from_coefficients(config.clone(), coeffs)
            };
            let b = {
                let mut coeffs = vec![0.0; basis_count];
                coeffs[b_idx as usize] = b_c;
                MultiVector::from_coefficients(config.clone(), coeffs)
            };

            let ab = a.scalar_product(&b);
            let ba = b.scalar_product(&a);

            let diff = (ab - ba).abs();
            prop_assert!(diff < 1e-8);
        });
    }

    #[test]
    fn test_scalar_product_involution_invariance() {
        proptest!(|(dim in 2u32..=5, a_idx in 0u64..8, a_c in -3.0f64..3.0f64,
                    b_idx in 0u64..8, b_c in -3.0f64..3.0f64)| {
            let config = test_config(dim);
            let basis_count = config.basis_count();
            let a_idx = a_idx.min(basis_count as u64 - 1);
            let b_idx = b_idx.min(basis_count as u64 - 1);

            let a = {
                let mut coeffs = vec![0.0; basis_count];
                coeffs[a_idx as usize] = a_c;
                MultiVector::from_coefficients(config.clone(), coeffs)
            };
            let b = {
                let mut coeffs = vec![0.0; basis_count];
                coeffs[b_idx as usize] = b_c;
                MultiVector::from_coefficients(config.clone(), coeffs)
            };

            let ab = a.scalar_product(&b);
            let a_star_b_star = a.grade_involution().scalar_product(&b.grade_involution());
            let a_dagger_b_dagger = a.reversion().scalar_product(&b.reversion());

            prop_assert!((ab - a_star_b_star).abs() < 1e-8);
            prop_assert!((ab - a_dagger_b_dagger).abs() < 1e-8);
        });
    }

    #[test]
    fn test_orthogonal_blades_scalar_product() {
        proptest!(|(dim in 3u32..=5, i in 0u32..4, j in 0u32..4)| {
            let config = test_config(dim);
            let i = i.min(dim - 1);
            let j = j.min(dim - 1);

            if i != j {
                let e_i = MultiVector::basis_vector(config.clone(), i);
                let e_j = MultiVector::basis_vector(config.clone(), j);

                let sp = e_i.scalar_product(&e_j);
                prop_assert!(sp.abs() < 1e-10);
            }
        });
    }
}

// ============================================================================
// INVERSE PROPERTIES
// ============================================================================

mod inverse_tests {
    use super::*;

    #[test]
    fn test_vector_inverse() {
        proptest!(|(dim in 2u32..=5, i in 0u32..4, c in 1.0f64..5.0f64)| {
            let config = test_config(dim);
            let i = i.min(dim - 1);

            let v = MultiVector::basis_vector(config.clone(), i).scale(c);

            let v_inv = v.inverse().expect("Non-zero vector should have inverse");

            let product = v.geometric_product(&v_inv);
            let scalar = product.scalar_part();

            prop_assert!((scalar - 1.0).abs() < 1e-8);
        });
    }

    #[test]
    fn test_inverse_of_inverse() {
        proptest!(|(dim in 2u32..=5, i in 0u32..4, c in 1.0f64..5.0f64)| {
            let config = test_config(dim);
            let i = i.min(dim - 1);

            let v = MultiVector::basis_vector(config.clone(), i).scale(c);

            let v_inv = v.inverse().expect("Vector should have inverse");
            let v_inv_inv = v_inv.inverse().expect("Inverse should have inverse");

            for k in 0..config.basis_count() {
                let diff = (v.get_coefficient(k) - v_inv_inv.get_coefficient(k)).abs();
                prop_assert!(diff < 1e-8);
            }
        });
    }

    #[test]
    fn test_bivector_inverse() {
        proptest!(|(dim in 3u32..=5, i in 0u32..3, j in 0u32..3, c in 1.0f64..3.0f64)| {
            let config = test_config(dim);
            let i = i.min(dim - 1);
            let j = j.min(dim - 1);

            if i != j {
                let e_i = MultiVector::basis_vector(config.clone(), i);
                let e_j = MultiVector::basis_vector(config.clone(), j);
                let bivector = e_i.outer_product(&e_j).scale(c);

                let b_inv = bivector.inverse().expect("Non-zero bivector should have inverse");

                let product = bivector.geometric_product(&b_inv);
                let scalar = product.scalar_part();

                prop_assert!((scalar - 1.0).abs() < 1e-8);
            }
        });
    }

    #[test]
    fn test_scalar_inverse() {
        proptest!(|(dim in 2u32..=5, c in 1.0f64..10.0f64)| {
            let config = test_config(dim);
            let scalar = MultiVector::from_scalar(config.clone(), c);

            let inv = scalar.inverse().expect("Non-zero scalar should have inverse");

            prop_assert!((inv.scalar_part() - 1.0/c).abs() < 1e-8);
        });
    }
}

// ============================================================================
// GEOMETRIC INTERPRETATION TESTS
// ============================================================================

mod geometry_tests {
    use super::*;

    #[test]
    fn test_reflection() {
        proptest!(|(dim in 2u32..=5, i in 0u32..4, j in 0u32..4)| {
            let config = test_config(dim);
            let i = i.min(dim - 1);
            let j = j.min(dim - 1);

            if i != j {
                let e_i = MultiVector::basis_vector(config.clone(), i);
                let e_j = MultiVector::basis_vector(config.clone(), j);

                let reflected = e_j.reflect_in(&e_i).expect("Reflection should succeed");

                let neg_e_j = e_j.neg();
                for k in 0..config.basis_count() {
                    let diff = (reflected.get_coefficient(k) - neg_e_j.get_coefficient(k)).abs();
                    prop_assert!(diff < 1e-8);
                }
            }
        });
    }

    #[test]
    fn test_double_reflection() {
        proptest!(|(dim in 3u32..=5)| {
            let config = test_config(dim);

            let e1 = MultiVector::basis_vector(config.clone(), 0);
            let e2 = MultiVector::basis_vector(config.clone(), 1);

            let reflected_once = e2.reflect_in(&e1).expect("First reflection should succeed");
            let reflected_twice = reflected_once.reflect_in(&e1).expect("Second reflection should succeed");

            for k in 0..config.basis_count() {
                let diff = (e2.get_coefficient(k) - reflected_twice.get_coefficient(k)).abs();
                prop_assert!(diff < 1e-8);
            }
        });
    }
}

#[test]
fn debug_inner_product_scalar() {
    let config = test_config(3);

    let scalar = MultiVector::from_scalar(config.clone(), 2.0);
    let vec = MultiVector::basis_vector(config.clone(), 0).scale(3.0);

    let result = scalar.left_inner(&vec);
    println!("s⌋a result scalar_part: {}", result.scalar_part());
    for i in 0..config.basis_count() {
        let c = result.get_coefficient(i);
        if c.abs() > 1e-10 {
            println!("  coefficient {}: {}", i, c);
        }
    }

    let result2 = vec.left_inner(&scalar);
    println!("a⌋s result scalar_part: {}", result2.scalar_part());
    for i in 0..config.basis_count() {
        let c = result2.get_coefficient(i);
        if c.abs() > 1e-10 {
            println!("  coefficient {}: {}", i, c);
        }
    }
}

#[test]
fn debug_dual_relationship() {
    let config = test_config(3);

    // Test scalar (idx=0)
    let scalar = MultiVector::from_scalar(config.clone(), 1.0);
    let dual_scalar = scalar.dual();
    let inv_dual = dual_scalar.inverse_dual();
    println!("Scalar: idx=0");
    println!("  dual(1):");
    for i in 0..config.basis_count() {
        let c = dual_scalar.get_coefficient(i);
        if c.abs() > 1e-10 {
            println!("    coeff[{}] = {}", i, c);
        }
    }
    println!("  inverse_dual(dual(1)):");
    for i in 0..config.basis_count() {
        let c = inv_dual.get_coefficient(i);
        if c.abs() > 1e-10 {
            println!("    coeff[{}] = {}", i, c);
        }
    }

    // Test vector (idx=1, e1)
    let e1 = MultiVector::basis_vector(config.clone(), 0);
    let dual_e1 = e1.dual();
    let inv_dual_e1 = dual_e1.inverse_dual();
    println!("\nVector e1: idx=1");
    println!("  dual(e1):");
    for i in 0..config.basis_count() {
        let c = dual_e1.get_coefficient(i);
        if c.abs() > 1e-10 {
            println!("    coeff[{}] = {}", i, c);
        }
    }
    println!("  inverse_dual(dual(e1)):");
    for i in 0..config.basis_count() {
        let c = inv_dual_e1.get_coefficient(i);
        if c.abs() > 1e-10 {
            println!("    coeff[{}] = {}", i, c);
        }
    }

    // Check I^2 for 3D
    let I = config.volume_element();
    let I_sq = I.geometric_product(&I);
    println!("\nI^2 (volume element squared):");
    println!("  scalar part: {}", I_sq.scalar_part());
}
