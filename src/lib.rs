
fn sqrt(n: u32) -> u32 {
    let mut x = n;
    let mut c: u32 = 0;

    let mut d: u32 = 1 << 30;

    while d < n {
        d = d >> 2;
    }

    while d != 0 {
        if x >= c + d {      
            x = x - (c + d);        
            c = (c >> 1) + d;  
        } else {
            c = c >> 1;
        }
        d = d >> 2;
    }

    return c;

}

fn float_sqrt (n: f64) -> f64 {

    const MANTISSA: u64 = 0xFFFFFFFFFFFFFF;
    const EXPONENT: u64 = 0xFFF0000000000000;

    let man = n.to_bits() & EXPONENT;
    let exp = ((n.to_bits() >> 52 ) as i64).wrapping_sub(1023);

    let mut new_exp = (exp as i64).wrapping_sub(1023);
    println!("{}", new_exp);

    new_exp = new_exp >> 2;
    let mut final_exp = new_exp.wrapping_add(1023) as u64;
    final_exp = final_exp & EXPONENT;

    let result = f64::from_bits(man | (final_exp  << 52) );

    return result;
}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn it_works() {
        assert_eq!(16, sqrt(256));
    }

    #[test]
    fn float_sqrt_test() {
        assert_eq!(16.0, float_sqrt(256.0));
    }
}