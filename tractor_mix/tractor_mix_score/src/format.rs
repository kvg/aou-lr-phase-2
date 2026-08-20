/// Match R `round(x, digits)`.
pub fn round_digits(x: f64, digits: i32) -> f64 {
    if !x.is_finite() {
        return x;
    }
    let factor = 10_f64.powi(digits);
    (x * factor).round() / factor
}

/// Match R `signif(x, digits)` for finite values.
pub fn signif_digits(x: f64, digits: i32) -> f64 {
    if !x.is_finite() {
        return x;
    }
    if x == 0.0 {
        return 0.0;
    }
    let exp = x.abs().log10().floor() as i32;
    let factor = 10_f64.powi(digits - 1 - exp);
    (x * factor).round() / factor
}

/// Upper tail chi-square CDF: P(X > x) for X ~ chi2(df).
pub fn pchisq_upper(x: f64, df: usize) -> f64 {
    if !x.is_finite() || x < 0.0 {
        return f64::NAN;
    }
    if df == 0 {
        return f64::NAN;
    }
    1.0 - chi2_cdf(x, df)
}

fn chi2_cdf(x: f64, df: usize) -> f64 {
    if x <= 0.0 {
        return 0.0;
    }
    regularized_gamma(df as f64 / 2.0, x / 2.0)
}

/// Lower regularized incomplete gamma P(a, x) = gamma(a, x) / Gamma(a).
fn regularized_gamma(a: f64, x: f64) -> f64 {
    if x <= 0.0 {
        return 0.0;
    }
    if x < a + 1.0 {
        series_gamma(a, x)
    } else {
        1.0 - continued_fraction_gamma(a, x)
    }
}

fn series_gamma(a: f64, x: f64) -> f64 {
    let mut ap = a;
    let mut sum = 1.0 / a;
    let mut del = sum;
    for _ in 0..200 {
        ap += 1.0;
        del *= x / ap;
        sum += del;
        if del.abs() < sum.abs() * 3e-7 {
            break;
        }
    }
    sum * (-x + a * x.ln() - ln_gamma(a)).exp()
}

fn continued_fraction_gamma(a: f64, x: f64) -> f64 {
    let mut b = x + 1.0 - a;
    let mut c = 1.0 / f64::MIN_POSITIVE;
    let mut d = 1.0 / b;
    let mut h = d;
    for i in 1..=200 {
        let an = -i as f64 * (i as f64 - a);
        b += 2.0;
        d = an * d + b;
        if d.abs() < f64::MIN_POSITIVE {
            d = f64::MIN_POSITIVE;
        }
        c = b + an / c;
        if c.abs() < f64::MIN_POSITIVE {
            c = f64::MIN_POSITIVE;
        }
        d = 1.0 / d;
        let del = d * c;
        h *= del;
        if (del - 1.0).abs() < 3e-7 {
            break;
        }
    }
    (-x + a * x.ln() - ln_gamma(a)).exp() * h
}

fn ln_gamma(z: f64) -> f64 {
    let coeffs = [
        76.18009172947146,
        -86.50532032941677,
        24.01409824083091,
        -1.231739572450155,
        0.001208650973866179,
        -0.000005395049384995,
    ];
    let x = z;
    let mut y = z;
    let mut tmp = x + 5.5;
    tmp -= (x + 0.5) * tmp.ln();
    let mut ser = 1.000000000190015;
    for c in coeffs.iter() {
        y += 1.0;
        ser += c / y;
    }
    -tmp + (2.5066282746310005 * ser / x).ln()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn round_matches_r_style() {
        assert_eq!(round_digits(1.234567, 5), 1.23457);
        assert!(round_digits(f64::NAN, 5).is_nan());
    }

    #[test]
    fn signif_matches_r_style() {
        let v = signif_digits(0.001234567, 5);
        assert!((v - 0.0012346).abs() < 1e-10);
    }

    #[test]
    fn pchisq_upper_sanity() {
        let p = pchisq_upper(3.841, 1);
        assert!((p - 0.05).abs() < 0.01);
    }
}
