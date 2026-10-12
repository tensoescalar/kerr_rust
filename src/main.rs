use num_complex::Complex;
use nalgebra::DVector;
use plotters::prelude::*;
use std::f64::consts::PI;

type C64 = Complex<f64>;

// Estado squeezed vacío analítico en la base de Fock
fn get_squeezed_vector(d: usize, r: f64) -> DVector<C64> {
    let mut state = DVector::<C64>::zeros(d);
    let cosh_r = r.cosh();
    let tanh_r = r.tanh();
    let prefactor = 1.0 / cosh_r.sqrt();

    for n in (0..d).step_by(2) {
        let k = n / 2;

        let mut fact_n = 1.0;
        for i in 1..=n {
            fact_n *= i as f64;
        }

        let mut fact_k = 1.0;
        for i in 1..=k {
            fact_k *= i as f64;
        }

        let sign = if k % 2 == 1 { -1.0 } else { 1.0 };
        let coeff = prefactor
            * sign
            * (fact_n.sqrt() / (fact_k * 2.0_f64.powi(k as i32)))
            * tanh_r.powi(k as i32);

        state[n] = C64::new(coeff, 0.0);
    }

    // Normalización (versión segura sin errores de tipo)
    let norm = state.norm();
    if norm > 1e-12 {
        for i in 0..state.len() {
            state[i] /= C64::new(norm, 0.0);
        }
    }

    state
}

// Aplica rotación real: c_n -> c_n * exp(-i * theta * n)
fn apply_rotation(state: &DVector<C64>, theta: f64) -> DVector<C64> {
    let mut new_state = state.clone();
    for n in 0..state.len() {
        let phase = C64::new(0.0, -theta * (n as f64));
        new_state[n] *= phase.exp();
    }
    new_state
}

// Aplica Kerr real: c_n -> c_n * exp(-i * chi * n * (n-1))
fn apply_kerr(state: &DVector<C64>, chi: f64) -> DVector<C64> {
    let mut new_state = state.clone();
    for n in 0..state.len() {
        let n_f = n as f64;
        let phase = C64::new(0.0, -chi * n_f * (n_f - 1.0));
        new_state[n] *= phase.exp();
    }
    new_state
}

// Función de Husimi Q
fn husimi_q(state: &DVector<C64>, alpha_re: f64, alpha_im: f64) -> f64 {
    let d = state.len();
    let mut overlap = C64::new(0.0, 0.0);
    let alpha = C64::new(alpha_re, alpha_im);
    let norm_factor = (-0.5 * alpha.norm_sqr()).exp();

    let mut fact = 1.0;
    let mut alpha_pow = C64::new(1.0, 0.0);

    for n in 0..d {
        if n > 0 {
            fact *= n as f64;
            alpha_pow *= alpha;
        }
        let coeff = state[n] * norm_factor * alpha_pow / fact.sqrt();
        overlap += coeff;
    }
    overlap.norm_sqr() / PI
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("--- SIMULADOR FOTONICO CV REAL ---");
    println!("Squeezing analitico + Rotacion real + Kerr real (n(n-1))");

    let d: usize = 40;
    let grid: usize = 55;
    let range: f64 = 4.8;

    let r = 0.85;
    let theta = PI / 4.0;
    let chi = 0.25;

    // 1. Estado inicial: vacio
    let mut state0 = DVector::<C64>::zeros(d);
    state0[0] = C64::new(1.0, 0.0);

    // 2. Squeezing analitico
    let state1 = get_squeezed_vector(d, r);

    // 3. Rotacion real
    let state2 = apply_rotation(&state1, theta);

    // 4. Kerr real
    let state3 = apply_kerr(&state2, chi);

    let states = [&state0, &state1, &state2, &state3];

    println!("Estados generados correctamente.");
    println!("Generando Husimi Q...");

    // Calcular los 4 mapas
    let mut maps = vec![vec![vec![0.0; grid]; grid]; 4];

    for (stage, st) in states.iter().enumerate() {
        for i in 0..grid {
            for j in 0..grid {
                let x = -range + (i as f64) * (2.0 * range) / ((grid - 1) as f64);
                let y = -range + (j as f64) * (2.0 * range) / ((grid - 1) as f64);
                maps[stage][i][j] = husimi_q(st, x, y);
            }
        }
    }

    // Visualizacion
    let root_path = "kerr_evolution_real_4_stages.png";
    let root = BitMapBackend::new(root_path, (1300, 1000)).into_drawing_area();
    root.fill(&WHITE)?;

    let titles = [
        "1. Estado Inicial (Vacio)",
        "2. Squeezing real",
        "3. Rotacion real (45 deg)",
        "4. Kerr real (Banana)",
    ];

    let areas = root.split_evenly((2, 2));

    for idx in 0..4 {
        let mut chart = ChartBuilder::on(&areas[idx])
            .caption(titles[idx], ("sans-serif", 18).into_font())
            .margin(12)
            .x_label_area_size(28)
            .y_label_area_size(28)
            .build_cartesian_2d(-range..range, -range..range)?;

        chart.configure_mesh().draw()?;

        let map_ref = &maps[idx];
        let max_val = map_ref
            .iter()
            .flatten()
            .cloned()
            .fold(0.0_f64, f64::max);

        chart.draw_series((0..grid).flat_map(move |i| {
            (0..grid).map(move |j| {
                let p = if max_val > 1e-12 {
                    map_ref[i][j] / max_val
                } else {
                    0.0
                };
                let x = -range + (i as f64) * (2.0 * range) / ((grid - 1) as f64);
                let y = -range + (j as f64) * (2.0 * range) / ((grid - 1) as f64);
                let color = HSLColor(240.0 / 360.0 - p * 0.75, 1.0, 0.45);
                Circle::new((x, y), (p * 5.0).max(0.6), color.filled())
            })
        }))?;
    }

    root.present()?;
    println!("Imagen guardada en: {}", root_path);
    println!("--- Listo ---");

    Ok(())
}
