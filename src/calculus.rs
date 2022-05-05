// Integration through Simpson's rule
pub fn integrate<F: Fn(f64) -> f64>(f: F, a: f64, b: f64, precision: i32) -> f64 {
    let delta_x = (b - a) / precision as f64;
    let mut sum = f(a) + f(b);

    for i in 1..precision {
        let t = a + i as f64 * delta_x;

        if i % 2 == 0 {
            sum += 2.0 * f(t);
        } else {
            sum += 4.0 * f(t);
        }
    }

    sum * delta_x / 3.0
}

#[cfg(test)]
mod tests {

    use super::*;
    use std::f64::consts::*;

    #[test]
    fn integrate_test() {
        // Test on sine function
        let sine_integral = integrate(f64::sin, 0.0, PI, 350);
        let sine_integral = (sine_integral * 10000.0).round() / 10000.0;
        assert_eq!(sine_integral, 2.0);

        // Test on parabula
        let parabula = |x: f64| -> f64 { x * x };
        let parabula_integral = integrate(parabula, 0.0, 1.0, 256);
        let parabula_integral = (parabula_integral * 10000.0).round() / 10000.0;
        assert_eq!(parabula_integral, 0.3333);
    }
}
