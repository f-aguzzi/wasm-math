use wasm_bindgen::prelude::*;
use fast_inv_sqrt::InvSqrt64;
use std::f64::consts::*;

#[wasm_bindgen]
pub fn normpdf(x: f64, mu: f64, sigma: f64) -> f64 {
    // f64::exp(-0.5 * f64::powi((x - mu) / sigma, 2)) / (sigma * f64::sqrt(2.0 * PI))
    s_normpdf(x - mu) / sigma
}

#[wasm_bindgen]
pub fn s_normpdf(x: f64) -> f64 {
    // Using the Quake fast inverse square root
    (2.0 * PI as f64).inv_sqrt64() * f64::exp(-0.5 * f64::powi(x, 2))
}

#[wasm_bindgen]
pub fn normcdf(x: f64, mu: f64, sigma: f64) -> f64 {
    s_normcdf((x - mu) / sigma)
}

#[wasm_bindgen]
pub fn s_normcdf(x: f64) -> f64 {
    // Using the Zelen & Severo's approximate algorithm
    let t: f64 = 1.0 / (1.0 + 0.2316419 * x);
    1.0 - s_normpdf(x) * (0.319381530 * t - 0.356563782 * f64::powi(t, 2) + 1.781477937 * f64::powi(t, 3) - 1.821255978 * f64::powi(t, 4) + 1.330274429 * f64::powi(t, 5))
}

#[wasm_bindgen]
pub fn s_norminv(x: f64) -> f64 {
    // Using Shore's approximate algorithm
    if x < 0.5 {
        -5.5556 * (1.0 - f64::powf(x / (1.0 - x), 0.1186))
    } else {
        5.5556 * (1.0 - f64::powf((1.0 - x) / x, 0.1186))
    }
}

#[wasm_bindgen]
pub fn norminv(x: f64, mu: f64, sigma: f64) -> f64 {
    mu + s_norminv(x) * sigma
}


#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn normpdf_test() {
        // Does the normal probability function compute the correct result?
        assert_eq!(normpdf(1.0, 1.0, 1.0), 0.39831039613785035);
    }

    #[test]
    fn s_normpdf_test() {
        // Does the normal probability function compute the correct result?

        //  Approximate the result 
        // (two different square root algorithms, so slightly different results are expected)
        let s_norm: String = s_normpdf(1.0).to_string().chars().take(5).collect();
        let norm: String = normpdf(1.0, 0.0, 1.0).to_string().chars().take(5).collect();

        assert_eq!(s_norm, norm);
    }

    #[test]
    fn normcdf_test() {
        //  Approximate the result 
        // (fast square root + approximate formula, so slightly different results are expected)
        let normcdf: String = normcdf(0.5, 1.0, 2.0).to_string().chars().take(5).collect();

        // Does the cumulative probability function compute the correct result?
        assert_eq!(normcdf, "0.402");
    }

    #[test]
    fn s_normcdf_test() {

        //  Approximate the result 
        // (fast square root + approximate formula, so slightly different results are expected)
        let normcdf: String = s_normcdf(0.0).to_string().chars().take(5).collect();

        // Does the cumulative probability function compute the correct result?
        assert_eq!(normcdf, "0.500");
    }

    #[test]
    fn s_norminv_test() {
        //  Approximate the result 
        // (fast square root + approximate formula, so slightly different results are expected)
        let norminv: String = s_norminv(0.3).to_string().chars().take(6).collect();

        // Does the inverse probability function compute the correct result?
        assert_eq!(s_norminv(0.5), 0.0);
        assert_eq!(norminv, "-0.531");
    }

    #[test]
    fn norminv_test() {
        // Approximate the results
        let norminv: String = norminv(0.3, 4.0, 2.0).to_string().chars().take(5).collect();

        // Does the inverse probability function compute the correct result?
        assert_eq!(norminv, "2.937")
    }
}