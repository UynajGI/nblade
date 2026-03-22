//! 缓冲区池性能基准测试

use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use nblade::{AlgebraConfig, MultiVector};
use std::sync::Arc;

#[cfg(feature = "pool")]
use nblade::multivector::pool::BufferPool;

fn geometric_product_3d_pooled(c: &mut Criterion) {
    let config = Arc::new(AlgebraConfig::euclidean(3));
    let e1 = MultiVector::basis_vector(config.clone(), 0);
    let e2 = MultiVector::basis_vector(config.clone(), 1);

    #[cfg(feature = "pool")]
    BufferPool::clear();

    c.bench_function("geometric_product_3d_pooled", |b| {
        b.iter(|| {
            let _result = e1.geometric_product(&e2);
        })
    });
}

fn geometric_product_5d_pooled(c: &mut Criterion) {
    let config = Arc::new(AlgebraConfig::euclidean(5));
    let e1 = MultiVector::basis_vector(config.clone(), 0);
    let e2 = MultiVector::basis_vector(config.clone(), 1);

    #[cfg(feature = "pool")]
    BufferPool::clear();

    c.bench_function("geometric_product_5d_pooled", |b| {
        b.iter(|| {
            let _result = e1.geometric_product(&e2);
        })
    });
}

fn geometric_product_chain_3d(c: &mut Criterion) {
    let config = Arc::new(AlgebraConfig::euclidean(3));
    let e1 = MultiVector::basis_vector(config.clone(), 0);
    let e2 = MultiVector::basis_vector(config.clone(), 1);
    let e3 = MultiVector::basis_vector(config.clone(), 2);

    #[cfg(feature = "pool")]
    BufferPool::clear();

    c.bench_function("geometric_product_chain_3d", |b| {
        b.iter(|| {
            let r1 = e1.geometric_product(&e2);
            let r2 = r1.geometric_product(&e3);
            let r3 = r2.geometric_product(&e1);
            black_box(r3);
        })
    });
}

fn geometric_product_many_operations(c: &mut Criterion) {
    let mut group = c.benchmark_group("many_operations");

    for dim in [3, 4, 5].iter() {
        let config = Arc::new(AlgebraConfig::euclidean(*dim));
        let vectors: Vec<_> = (0..*dim as u32)
            .map(|i| MultiVector::basis_vector(config.clone(), i))
            .collect();

        #[cfg(feature = "pool")]
        BufferPool::clear();

        group.bench_with_input(BenchmarkId::new("dim", dim), dim, |b, _| {
            b.iter(|| {
                let mut result = vectors[0].clone();
                for i in 1..vectors.len() {
                    result = result.geometric_product(&vectors[i]);
                }
                black_box(result);
            })
        });
    }

    group.finish();
}

#[cfg(feature = "pool")]
fn pool_buffer_operations(c: &mut Criterion) {
    BufferPool::clear();

    c.bench_function("pool_get_return_8", |b| {
        b.iter(|| {
            let buf = BufferPool::get_buffer(8);
            BufferPool::return_buffer(buf);
        })
    });

    c.bench_function("pool_get_return_64", |b| {
        b.iter(|| {
            let buf = BufferPool::get_buffer(64);
            BufferPool::return_buffer(buf);
        })
    });

    c.bench_function("pool_get_return_256", |b| {
        b.iter(|| {
            let buf = BufferPool::get_buffer(256);
            BufferPool::return_buffer(buf);
        })
    });
}

#[cfg(feature = "pool")]
criterion_group!(
    benches,
    geometric_product_3d_pooled,
    geometric_product_5d_pooled,
    geometric_product_chain_3d,
    geometric_product_many_operations,
    pool_buffer_operations,
);

#[cfg(not(feature = "pool"))]
criterion_group!(
    benches,
    geometric_product_3d_pooled,
    geometric_product_5d_pooled,
    geometric_product_chain_3d,
    geometric_product_many_operations,
);

criterion_main!(benches);
