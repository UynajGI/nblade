//! 几何积基准测试
//! 包含低维（SIMD/顺序执行）和高维（并行执行）的性能对比

use criterion::{black_box, criterion_group, criterion_main, Criterion};
use nblade::{AlgebraConfig, MultiVector};
use std::sync::Arc;

fn geometric_product_2d(c: &mut Criterion) {
    let config = Arc::new(AlgebraConfig::euclidean(2));
    let e1 = MultiVector::basis_vector(config.clone(), 0);
    let e2 = MultiVector::basis_vector(config.clone(), 1);

    c.bench_function("geometric_product_2d", |b| {
        b.iter(|| {
            let _result = e1.geometric_product(&e2);
        })
    });
}

fn geometric_product_3d(c: &mut Criterion) {
    let config = Arc::new(AlgebraConfig::euclidean(3));
    let e1 = MultiVector::basis_vector(config.clone(), 0);
    let e2 = MultiVector::basis_vector(config.clone(), 1);

    c.bench_function("geometric_product_3d", |b| {
        b.iter(|| {
            let _result = e1.geometric_product(&e2);
        })
    });
}

fn geometric_product_4d(c: &mut Criterion) {
    let config = Arc::new(AlgebraConfig::euclidean(4));
    let e1 = MultiVector::basis_vector(config.clone(), 0);
    let e2 = MultiVector::basis_vector(config.clone(), 1);

    c.bench_function("geometric_product_4d", |b| {
        b.iter(|| {
            let _result = e1.geometric_product(&e2);
        })
    });
}

fn geometric_product_full_3d(c: &mut Criterion) {
    let config = Arc::new(AlgebraConfig::euclidean(3));
    let coeffs_a: Vec<f64> = (0..8).map(|i| (i as f64 + 1.0) * 0.1).collect();
    let coeffs_b: Vec<f64> = (0..8).map(|i| (i as f64 + 2.0) * 0.1).collect();
    let a = MultiVector::from_coefficients(config.clone(), coeffs_a);
    let b_vec = MultiVector::from_coefficients(config.clone(), coeffs_b);

    c.bench_function("geometric_product_full_3d", |b| {
        b.iter(|| {
            let _result = a.geometric_product(&b_vec);
        })
    });
}

fn geometric_product_full_4d(c: &mut Criterion) {
    let config = Arc::new(AlgebraConfig::euclidean(4));
    let coeffs_a: Vec<f64> = (0..16).map(|i| (i as f64 + 1.0) * 0.1).collect();
    let coeffs_b: Vec<f64> = (0..16).map(|i| (i as f64 + 2.0) * 0.1).collect();
    let a = MultiVector::from_coefficients(config.clone(), coeffs_a);
    let b_vec = MultiVector::from_coefficients(config.clone(), coeffs_b);

    c.bench_function("geometric_product_full_4d", |b| {
        b.iter(|| {
            let _result = a.geometric_product(&b_vec);
        })
    });
}

fn geometric_product_5d(c: &mut Criterion) {
    let config = Arc::new(AlgebraConfig::euclidean(5));
    let e1 = MultiVector::basis_vector(config.clone(), 0);
    let e2 = MultiVector::basis_vector(config.clone(), 1);

    c.bench_function("geometric_product_5d", |b| {
        b.iter(|| {
            let _result = e1.geometric_product(&e2);
        })
    });
}

fn geometric_product_6d(c: &mut Criterion) {
    let config = Arc::new(AlgebraConfig::euclidean(6));
    let e1 = MultiVector::basis_vector(config.clone(), 0);
    let e2 = MultiVector::basis_vector(config.clone(), 1);

    c.bench_function("geometric_product_6d", |b| {
        b.iter(|| {
            let _result = e1.geometric_product(&e2);
        })
    });
}

fn geometric_product_8d(c: &mut Criterion) {
    let config = Arc::new(AlgebraConfig::euclidean(8));
    let e1 = MultiVector::basis_vector(config.clone(), 0);
    let e2 = MultiVector::basis_vector(config.clone(), 1);

    c.bench_function("geometric_product_8d", |b| {
        b.iter(|| {
            let _result = e1.geometric_product(&e2);
        })
    });
}

fn geometric_product_10d(c: &mut Criterion) {
    let config = Arc::new(AlgebraConfig::euclidean(10));
    let e1 = MultiVector::basis_vector(config.clone(), 0);
    let e2 = MultiVector::basis_vector(config.clone(), 1);

    c.bench_function("geometric_product_10d", |b| {
        b.iter(|| {
            let _result = e1.geometric_product(&e2);
        })
    });
}

fn geometric_product_full_6d(c: &mut Criterion) {
    let config = Arc::new(AlgebraConfig::euclidean(6));
    let coeffs_a: Vec<f64> = (0..64).map(|i| (i as f64 + 1.0) * 0.1).collect();
    let coeffs_b: Vec<f64> = (0..64).map(|i| (i as f64 + 2.0) * 0.1).collect();
    let a = MultiVector::from_coefficients(config.clone(), coeffs_a);
    let b_vec = MultiVector::from_coefficients(config.clone(), coeffs_b);

    c.bench_function("geometric_product_full_6d", |bencher| {
        bencher.iter(|| {
            let _result = a.geometric_product(&b_vec);
        })
    });
}

fn geometric_product_full_8d(c: &mut Criterion) {
    let config = Arc::new(AlgebraConfig::euclidean(8));
    let coeffs_a: Vec<f64> = (0..256).map(|i| (i as f64 + 1.0) * 0.01).collect();
    let coeffs_b: Vec<f64> = (0..256).map(|i| (i as f64 + 2.0) * 0.01).collect();
    let a = MultiVector::from_coefficients(config.clone(), coeffs_a);
    let b_vec = MultiVector::from_coefficients(config.clone(), coeffs_b);

    c.bench_function("geometric_product_full_8d", |bencher| {
        bencher.iter(|| {
            let _result = a.geometric_product(&b_vec);
        })
    });
}

criterion_group!(
    benches,
    geometric_product_2d,
    geometric_product_3d,
    geometric_product_4d,
    geometric_product_5d,
    geometric_product_6d,
    geometric_product_8d,
    geometric_product_10d,
    geometric_product_full_3d,
    geometric_product_full_4d,
    geometric_product_full_6d,
    geometric_product_full_8d,
);
criterion_main!(benches);
