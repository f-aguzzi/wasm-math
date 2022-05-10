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

// Naïve differentiation algorithm
pub fn differentiate<F: Fn(f64) -> f64>(f: F, a: f64) -> f64 {
    let delta_x = 0.0000000001;

    (f(a + delta_x) - f(a)) / delta_x
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

    #[test]
    fn differentiate_test() {
        let parabula = |x: f64| -> f64 { x * x };

        let parabula_derivative_1 = differentiate(parabula, 0.0);
        let parabula_derivative_1 = (parabula_derivative_1 * 10000.0).round() / 10000.0;
        assert_eq!(parabula_derivative_1, 0.0);

        let parabula_derivative_2 = differentiate(parabula, 1.0);
        let parabula_derivative_2 = (parabula_derivative_2 * 10000.0).round() / 10000.0;
        assert_eq!(parabula_derivative_2, 2.0);

        let parabula_derivative_3 = differentiate(parabula, -2.5);
        let parabula_derivative_3 = (parabula_derivative_3 * 10000.0).round() / 10000.0;
        assert_eq!(parabula_derivative_3, -5.0);
    }
}
