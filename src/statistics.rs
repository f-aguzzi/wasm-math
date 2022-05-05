use crate::calculus::integrate;
use crate::optimizers::fzero;
use std::f64::consts::*;
use wasm_bindgen::prelude::*;

/**
 *  ------------------------------------------------------
 *  NORMAL DISTRIBUTION
 *  Gaussian distribitions (normalized and non-normalized)
 *  Density, Cumulate and Quantile functions
 *  ------------------------------------------------------
 */

 // Gaussian distribution
#[wasm_bindgen]
pub fn normpdf(x: f64, mu: f64, sigma: f64) -> f64 {
    // f64::exp(-0.5 * f64::powi((x - mu) / sigma, 2)) / (sigma * f64::sqrt(2.0 * PI))
    s_normpdf(x - mu) / sigma
}

// Normal distribution
#[wasm_bindgen]
pub fn s_normpdf(x: f64) -> f64 {
    f64::exp(-0.5 * f64::powi(x, 2)) / f64::sqrt(2.0 * PI)
}

// Gaussian cumulate distribution
#[wasm_bindgen]
pub fn normcdf(x: f64, mu: f64, sigma: f64) -> f64 {
    s_normcdf((x - mu) / sigma)
}

// Standard cumulate distribution
#[wasm_bindgen]
pub fn s_normcdf(x: f64) -> f64 {
    // Using the Zelen & Severo's approximate algorithm
    let t: f64 = 1.0 / (1.0 + 0.2316419 * x);
    1.0 - s_normpdf(x)
        * (0.319381530 * t - 0.356563782 * f64::powi(t, 2) + 1.781477937 * f64::powi(t, 3)
            - 1.821255978 * f64::powi(t, 4)
            + 1.330274429 * f64::powi(t, 5))
}

// Standard distribution quantile
#[wasm_bindgen]
pub fn s_norminv(x: f64) -> f64 {
    // Using Shore's approximate algorithm
    if x < 0.5 {
        -5.5556 * (1.0 - f64::powf(x / (1.0 - x), 0.1186))
    } else {
        5.5556 * (1.0 - f64::powf((1.0 - x) / x, 0.1186))
    }
}

// Gaussian quantile
#[wasm_bindgen]
pub fn norminv(x: f64, mu: f64, sigma: f64) -> f64 {
    mu + s_norminv(x) * sigma
}

/**
 *  --------------------------------------------------------------
 *  GAMMA DISTRIBUTION
 *  Gamma functions (complete, lower incomplete, upper incomplete)
 *  Gamma distribution (pdf, cdf, quantile)
 *  --------------------------------------------------------------
 */

// Gamma function
#[wasm_bindgen]
pub fn gamma(x: f64) -> f64 {
    // coefficients from Chebyshev polynomials
    let q = vec![
        75122.6331530,
        80916.6278952,
        36308.2951477,
        8687.24529705,
        1168.92649479,
        83.8676043424,
        2.50662827511,
    ];

    let mut num: f64 = 0.0;
    let mut den: f64 = 1.0;

    for n in 0..=6 {
        num += q[n] * x.powi(n as i32);
        den *= x + (n as f64);
    }

    // Lanczos's approximation
    num * den.recip() * (x + 5.5).powf(x + 0.5) * f64::exp(-(x + 5.5))
}

// Lower incomplete gamma function
#[wasm_bindgen]
pub fn lowincgamma(s: f64, x: f64) -> f64 {
    gamma(s) * gammapdf(x, s, 1.0)
}

// Upper incomplete gamma function
#[wasm_bindgen]
pub fn uppincgamma(s: f64, x: f64) -> f64 {
    gamma(s) * (1.0 - gammapdf(x, s, 1.0))
}

// Regularized incomplete gamma function
#[wasm_bindgen]
pub fn regincgamma(s: f64, x: f64) -> f64 {
    gamma(s) * (1.0 - gammapdf(x, s, 1.0))
}

// Gamma distribution
#[wasm_bindgen]
pub fn gammapdf(x: f64, a: f64, b: f64) -> f64 {
    if x < 0.0 {
        return f64::NAN;
    }
    x.powf(a - 1.0) * f64::exp(-b * x) * b.powf(a) / gamma(a)
}

// Cumulate gamma distribution
#[wasm_bindgen]
pub fn gammacdf(x: f64, a: f64, b: f64) -> f64 {
    unimplemented!()
}

// Ggamma distribution quantile
#[wasm_bindgen]
pub fn gammainv(x: f64, a: f64, b: f64) -> f64 {
    unimplemented!()
}


/**
 *  -------------------------------------------------
 *  BETA DISTRIBUTION
 *  Beta function (complete, incomplete, regularized)
 *  Beta distribution (pdf, cdf, quantile)
 *  -------------------------------------------------
 */

// Beta function
#[wasm_bindgen]
pub fn beta(x: f64, y: f64) -> f64 {
    gamma(x) * gamma(y) / gamma(x + y)
}

// Incomplete Beta Function
#[wasm_bindgen]
pub fn incbet(x: f64, a: f64, b: f64) -> f64 {
    // I had to give up and use numerical integration
    let beta_deriv = |t: f64| -> f64 { t.powf(a - 1.0) * (1.0 - t).powf(b - 1.0) };
    integrate(beta_deriv, 0.0, x, 256)
}

// Regularized incomplete beta function
#[wasm_bindgen]
pub fn regincbet(x: f64, a: f64, b: f64) -> f64 {
    incbet(x, a, b) / beta(a, b)
}

// Beta distribution
#[wasm_bindgen]
pub fn betapdf(x: f64, a: f64, b: f64) -> f64 {
    unimplemented!()
}

// Cumulative beta distribution
#[wasm_bindgen]
pub fn betacdf(x: f64, a: f64, b: f64) -> f64 {
    unimplemented!()
}

// Beta distribution quantile
#[wasm_bindgen]
pub fn betainv(x: f64, a: f64, b: f64) -> f64 {
    unimplemented!()
}

/**
 *  ----------------------------------------
 *  STUDENT'S T DISTRIBUTION
 *  t distribitions (normalized)
 *  Density, Cumulate and Quantile functions
 *  ----------------------------------------
 */

// t distribution
#[wasm_bindgen]
pub fn tpdf(x: f64, v: f64) -> f64 {
    gamma((v + 1.0) * 0.5)
        * (v * PI).sqrt().recip()
        * gamma(v * 0.5).recip()
        * (1.0 + x * x * v.recip()).powf(-(v + 1.0) * 0.5)
}

// Cumulate t distribution
#[wasm_bindgen]
pub fn tcdf(x: f64, v: f64) -> f64 {
    if x == 0.0 {
        return 0.5;
    }

    let tdist2t = |t: f64, v: f64| -> f64 { 1.0 - regincbet(v / (v + t * t), v * 0.5, 0.5) };
    let tdist1t = |t: f64, v: f64| -> f64 { 1.0 - (1.0 - tdist2t(t, v)) * 0.5 };

    match x > 0.0 {
        true => tdist1t(x, v),
        false => 1.0 - tdist1t(-x, v),
    }
}

// t distribution quantile
#[wasm_bindgen]
pub fn tinv(x: f64, v: f64) -> f64 {
    let to_maximize = |tval: f64| -> f64 { x - tcdf(tval, v) };

    if x < i8::MIN as f64 {
        f64::NEG_INFINITY
    } else if x > i8::MAX as f64 {
        f64::INFINITY
    } else if x > 0.000000001 && x < 0.999999999 {
        fzero(to_maximize, -10.0, 10.0)
    } else {
        fzero(to_maximize, i8::MIN as f64, i8::MAX as f64)
    }
}

/**
*   ----------------------------------------
*   CHI SQUARED DISTRIBUTION
*   X^2 distribition
*   Density, Cumulate and Quantile functions
*   ----------------------------------------
*/

#[wasm_bindgen]
pub fn chi2pdf(x: f64, k: f64) -> f64 {
    if x > 0.0 {
        x.powf(0.5 * k - 1.0) * f64::exp(-0.5 * x) * (-0.5 * k).exp2() / gamma(0.5 * k)
    } else {
        0.0
    }
}

#[wasm_bindgen]
pub fn chi2cdf(x: f64, k: f64) -> f64 {
    gamma(x * 0.5).recip() * lowincgamma(k * 0.5, x * 0.5)
}

#[wasm_bindgen]
pub fn chi2inv(x: f64, k: f64) -> f64 {
    x + (4.0 * (k - 1.0)).recip() * x * x
}

/**
 * |----------------------------|
 * |----------------------------|
 * |         UNIT TESTS         |
 * |----------------------------|
 * |----------------------------|
 */

#[cfg(test)]
mod tests {

    use super::*;

    /**
     *  ------------------------------
     *  NORMAL DISTRIBUTION UNIT TESTS
     *  ------------------------------
     */

    #[test]
    fn normpdf_test() {
        // Does the normal probability function compute the correct result?
        assert_eq!(normpdf(1.0, 1.0, 1.0), 0.3989422804014327);
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
        // Does the cumulative probability function compute the correct result?
        assert_eq!(normcdf(0.5, 1.0, 2.0), 0.4012880670269864);
    }

    #[test]
    fn s_normcdf_test() {
        // Does the cumulative probability function compute the correct result?
        assert_eq!(s_normcdf(0.0), 0.5000000005248086);
    }

    #[test]
    fn s_norminv_test() {
        // Does the inverse probability function compute the correct result?
        assert_eq!(s_norminv(0.5), 0.0);
        assert_eq!(s_norminv(0.3), -0.531145444833719);
    }

    #[test]
    fn norminv_test() {
        // Does the inverse probability function compute the correct result?
        assert_eq!(norminv(0.3, 4.0, 2.0), 2.937709110332562)
    }

    /**
     *  -----------------------------
     *  GAMMA DISTRIBUTION UNIT TESTS
     *  -----------------------------
     */

    #[test]
    fn gamma_test() {
        // Extract the result
        let gamma = gamma(2.3);

        // Round to 5 digits
        let gamma = (gamma * 10000.0).round() / 10000.0;

        // Compare with MATLAB result
        assert_eq!(gamma, 1.1667);
    }

    #[test]
    fn lowincgamma_test() {
        // Extract the result
        let lowincgamma = lowincgamma(1.44, 3.5);

        // Round to 5 digits
        let lowincgamma = (lowincgamma * 10000.0).round() / 10000.0;

        // Compare with MATLAB result
        assert_eq!(lowincgamma, 0.1041);
    }

    #[test]
    fn uppincgamma_test() {
        // Extract the result
        let lowincgamma = uppincgamma(1.44, 3.5);

        // Round to 5 digits
        let lowincgamma = (lowincgamma * 10000.0).round() / 10000.0;

        // Compare with MATLAB result
        assert_eq!(lowincgamma, 0.8959);
    }

    #[test]
    fn gammapdf_test() {
        // Compute the result
        let pdf = gammapdf(0.93, 3.0, 7.0);

        // Round to 5 digits
        let pdf = (pdf * 10000.0).round() / 10000.0;

        // Compare with MATLAB result
        assert_eq!(pdf, 0.0011);
    }

    #[test]
    fn gammacdf_test() {
        // Compute the result
        let cdf = gammacdf(1.23, 2.0, 2.0);

        // Round to 5 digits
        let cdf = (cdf * 10000.0).round() / 10000.0;

        // Compare with MATLAB result
        assert_eq!(cdf, 0.1269);
    }

    #[test]
    fn gammainv_test() {
        // Compute the result
        let inv = gammainv(0.68, 3.0, 3.0);

        // Round to 5 digits
        let inv = (inv * 10000.0).round() / 10000.0;

        // Compare with MATLAB result
        assert_eq!(inv, 10.5138);
    }

    /**
     *  -----------------------
     *  BETA DISTRIBUTION TESTS
     *  -----------------------
     */

    #[test]
    fn beta_test() {
        // Extract the result
        let beta = beta(3.0, 2.0);

        // Round to 5 digits
        let beta = (beta * 10000.0).round() / 10000.0;

        // Compare with MATLAB result
        assert_eq!(beta, 0.0833);
    }

    #[test]
    fn incbet_test() {
        // Extract the result
        let incbet = incbet(0.3, 2.0, 6.0);

        // Round to 4 digits
        let incbet = (incbet * 1000.0).round() / 1000.0;

        // Compare with MATLAB result
        assert_eq!(incbet, 0.0160);
    }

    #[test]
    fn regincbet_test() {
        // Extract the result
        let regincbet = regincbet(0.3, 11.0, 7.0);

        // Round to 5 digits
        let regincbet = (regincbet * 10000.0).round() / 10000.0;

        // Compare with MATLAB result
        assert_eq!(regincbet, 0.0032);
    }

    #[test]
    fn betapdf_test() {
        // Compute the result
        let pdf = betapdf(0.01, 2.0, 5.0);

        // Round to 5 digits
        let pdf = (pdf * 10000.0).round() / 10000.0;

        // Compare with MATLAB result
        assert_eq!(pdf, 0.2882);
    }

    #[test]
    fn betacdf_test() {
        // Compute the result
        let cdf = betacdf(0.6, 9.0, 4.0);

        // Round to 5 digits
        let cdf = (cdf * 10000.0).round() / 10000.0;

        // Compare with MATLAB result
        assert_eq!(cdf, 0.2253);
    }

    #[test]
    fn betainv_test() {
        // Compute the result
        let inv = betainv(0.72, 3.0, 4.0);

        // Round to 5 digits
        let inv = (inv * 10000.0).round() / 10000.0;

        // Compare with MATLAB result
        assert_eq!(inv, 0.5354);
    }


    //  --------------------
    //  T DISTRIBUTION TESTS
    //  --------------------

    #[test]
    fn tpdf_test() {
        // Extract the result
        let tpdf = tpdf(0.8, 11.0);

        // Round to 5 digits
        let tpdf = (tpdf * 10000.0).round() / 10000.0;

        // Compare with MATLAB result
        assert_eq!(tpdf, 0.2778);
    }

    #[test]
    fn tcdf_test() {
        // Extract the first result
        let tcdf_1 = tcdf(3.2, 11.0);

        // Round to 5 digits
        let tcdf_1 = (tcdf_1 * 10000.0).round() / 10000.0;

        // Compare with MATLAB result
        assert_eq!(tcdf_1, 0.9958);

        // Extract second the result
        let tcdf_2 = tcdf(-3.2, 11.0);
        let tcdf_2 = (tcdf_2 * 10000.0).round() / 10000.0;
        assert_eq!(tcdf_2, 0.0042);

        // Third result
        let tcdf_3 = tcdf(0.0, 11.0);
        let tcdf_3 = (tcdf_3 * 10000.0).round() / 10000.0;
        assert_eq!(tcdf_3, 0.5);
    }

    #[test]
    fn tinv_test() {
        // Extract the result
        let tinv = tinv(0.4, 8.0);

        // Truncate to 5 digits
        let tinv = (tinv * 10000.0).trunc() / 10000.0;

        // Compare with MATLAB result
        assert_eq!(tinv, -0.2619);
    }

    // ------------------------
    // CHI^2 DISTRIBUTION TESTS
    // ------------------------

    #[test]
    fn chi2pdf_test() {
        // Compute the result
        let pdf = chi2pdf(0.56, 5.0);

        // Round to 5 digits
        let pdf = (pdf * 10000.0).round() / 10000.0;

        // Check against MATLAB result
        assert_eq!(pdf, 0.0421);
    }

    #[test]
    fn chi2cdf_test() {
        // Compute the result
        let cdf = chi2cdf(0.42, 3.0);

        // Round to 5 digits
        let cdf = (cdf * 10000.0).round() / 10000.0;

        // Check against MATLAB result
        assert_eq!(cdf, 0.0639);
    }

    #[test]
    fn chi2inv_test() {
        // Compute the result
        let inv = chi2inv(0.076, 4.0);

        // Round to 5 digits
        let inv = (inv * 10000.0).round() / 10000.0;

        // Check against MATLAB result
        assert_eq!(inv, 0.9039);
    }
}
