pub fn maximum<F: Fn(f64) -> f64>(f: F, left: f64, right: f64, precision: f64) -> f64 {
    let mut max = f(left);
    let mut cur_position = left;

    while cur_position < right {
        max = f64::max(max, f(cur_position));
        cur_position += precision;
    }

    max
}

pub fn fzero<F: Fn(f64) -> f64>(f: F, left: f64, right: f64) -> f64 {
    let precision = 2048;

    let step = (right - left) / precision as f64;
    let mut samples: Vec<f64> = Vec::new();
    let mut cur_position = left;

    while cur_position < right {
        samples.push(cur_position);
        cur_position += step;
    }

    let mut index = 0;
    while index < precision - 2 {
        if f(samples[index]).signum() != f(samples[index + 1]).signum() {
            break;
        }
        index += 1;
    }

    let find_min_max = |a: f64, b: f64| -> (f64, f64) {
        if f(a) < f(b) {
            (a, b)
        } else {
            (b, a)
        }
    };

    let min_max = find_min_max(left + step * index as f64, left + step * (index + 1) as f64);
    let mut min = min_max.0;
    let mut max = min_max.1;

    let mut mid = (max + min) * 0.5;

    for _i in 0..(precision / 8) {
        let midpoint = f(mid);

        if midpoint == 0.0 {
            break;
        } else if midpoint > 0.0 {
            max = mid;
        } else {
            min = mid;
        }

        mid = (max + min) * 0.5;
    }

    mid
}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn maximum_test() {
        let max = maximum(f64::cos, -1.0, 1.0, 0.0001);

        // Approximate the result to 5 digits
        let max = (max * 10000.0).round() / 10000.0;

        // Check if the result is correct
        assert_eq!(max, 1.0);
    }

    #[test]
    fn fzero_test() {
        // SINE FUNCTION
        // Look for the first zero in the -3.0 - +3.0 range
        let zero = fzero(f64::sin, -3.0, 3.0);
        let sin_zero = f64::sin(zero);

        // Approximate the result to 5 digits
        let sin_zero = (sin_zero * 10000.0).round() / 10000.0;

        assert_eq!(sin_zero, 0.0);

        // PARABULA
        let parabula = |x: f64| -> f64 { x * x - 3.0 * x + 1.0 };
        let zero = fzero(parabula, 0.0, 1.0);
        let parabula_zero = parabula(zero);
        assert_eq!(parabula_zero, 0.0);
    }
}
