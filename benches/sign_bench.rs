use criterion::{black_box, criterion_group, criterion_main, Criterion};
use pure_dsa::Algorithm; // adjust if your crate name is different

fn bench_sign_verify(c: &mut Criterion) {
    let msg = b"benchmark message";

    for alg in &[Algorithm::Mode2, Algorithm::Mode3, Algorithm::Mode5] {
        let name = match alg {
            Algorithm::Mode2 => "Mode2",
            Algorithm::Mode3 => "Mode3",
            Algorithm::Mode5 => "Mode5",
        };

        c.bench_function(&format!("{name} keygen"), |b| {
            b.iter(|| {
                black_box(alg.generate());
            });
        });

        let keypair = alg.generate();
        let public = keypair.public();

        c.bench_function(&format!("{name} sign"), |b| {
            b.iter(|| {
                let _ = keypair.sign(black_box(msg));
            });
        });

        let signature = keypair.sign(msg);

        c.bench_function(&format!("{name} verify"), |b| {
            b.iter(|| {
                let _ = alg.verify(black_box(&signature), black_box(msg), black_box(public));
            });
        });
    }
}

criterion_group!(benches, bench_sign_verify);
criterion_main!(benches);
