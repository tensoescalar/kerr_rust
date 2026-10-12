# Continuous-Variable (CV) Photonic Quantum Simulator in Rust

Simulador numérico de alta eficiencia en Rust para óptica cuántica de variables continuas (CV), modelando la evolución en el espacio de fases bidimensional $(\alpha_r, \alpha_i)$ en la base de Fock truncada.

---

## Características Físicas y Pipeline de Simulación

La plataforma fotónica ofrece ventajas críticas frente a sistemas superconductores, destacando su operación a temperatura ambiente y su alta resistencia a la decoherencia. Este repositorio implementa rigurosamente cuatro etapas físicas clave:

1. **Estado Inicial (Vacío):** Simetría de mínima incertidumbre centrada en el origen del espacio de fases.
2. **Squeezing Analítico:** Compresión de fluctuaciones cuánticas en una cuadratura específica mediante la representación analítica exacta en la base de Fock truncada.
3. **Rotación Real:** Evolución unitaria exacta mediante operadores diagonales de fase ($c_n \rightarrow c_n e^{-i\theta n}$).
4. **Efecto Kerr Óptico No Lineal:** Simulación del término de tercer orden ($n(n-1)$), generando la característica distorsión geométrica en forma de "banana" en el espacio de fases.

---

## Requisitos y Dependencias (`Cargo.toml`)

El proyecto utiliza las siguientes dependencias optimizadas para Rust:
- `nalgebra`: Manejo eficiente de vectores de estado complejos.
- `plotters`: Generación de gráficos de alta resolución.
- `num-complex`: Soporte avanzado para números complejos.

---

## Compilación y Ejecución

Para compilar y ejecutar el simulador en modo de alto rendimiento (`release`):

```bash
cargo run --release

