//! Классический пример из статьи Kettner, Mehlhorn, Pion, Schirra, Yap
//! «Classroom examples of robustness problems in geometric computations»:
//! точки около прямой y = x, сдвинутые на доли ulp. Наивная формула в f64
//! путает знаки ориентации, точный предикат — ни разу.

use yadro_core::Vec2;
use yadro_core::predicates::{Sign, orient2d};

#[test]
fn near_collinear_grid_of_ulps() {
    let q = Vec2::new(12.0, 12.0);
    let r = Vec2::new(24.0, 24.0);
    let ulp = f64::EPSILON / 2.0; // шаг представления чисел около 0,5: 2^−53
    let mut naive_wrong = 0;
    for i in 0..256u32 {
        for j in 0..256u32 {
            let p = Vec2::new(0.5 + f64::from(i) * ulp, 0.5 + f64::from(j) * ulp);
            // Точно: (q − p) × (r − p) = 12·2^−53·(j − i).
            let exact = Sign::of_i128(i128::from(j) - i128::from(i));
            assert_eq!(orient2d(p, q, r), exact, "i={i} j={j}");
            let naive = (q.x - p.x) * (r.y - p.y) - (q.y - p.y) * (r.x - p.x);
            if Sign::of_f64(naive) != exact {
                naive_wrong += 1;
            }
        }
    }
    println!("наивная формула ошиблась в {naive_wrong} точках из 65536");
    assert!(
        naive_wrong > 0,
        "пример должен показывать ошибки наивной формулы"
    );
}
