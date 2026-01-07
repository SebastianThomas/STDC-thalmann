use core::time::Duration;

const SECS_PER_HOUR: u64 = 60 * 60;
const SECS_PER_MIN: u64 = 60;

pub const fn show_duration(d: Duration) -> [char; 9] {
    let secs = d.as_secs();
    if secs > SECS_PER_HOUR {
        let hours = padded_3::<'0'>(secs / SECS_PER_HOUR);
        let mins = padded_2::<'0'>((secs % SECS_PER_HOUR) / SECS_PER_MIN);
        let secs = padded_2::<'0'>(secs % SECS_PER_MIN);
        return [
            hours[0], hours[1], hours[2], ':', mins[0], mins[1], ':', secs[0], secs[1],
        ];
    }
    let mins = padded_2::<'0'>((secs % SECS_PER_HOUR) / SECS_PER_MIN);
    let secs = padded_2::<'0'>(secs % SECS_PER_MIN);
    let millis = padded_3::<'0'>(d.subsec_millis() as u64);
    return [
        mins[0], mins[1], ':', secs[0], secs[1], '.', millis[0], millis[1], millis[2],
    ];
}

pub const fn padded_2<const C: char>(n: u64) -> [char; 2] {
    let most_significant_char: char = if n < 10 {
        C
    } else {
        let tens: u32 = ((n / 10) % 10) as u32;
        char::from_digit(tens, 10).expect("Tens digit should be < 10")
    };
    return [
        most_significant_char,
        char::from_digit((n % 10) as u32, 10).unwrap(),
    ];
}

pub const fn padded_3<const C: char>(n: u64) -> [char; 3] {
    let most_significant_char: char = {
        let hundreds: u32 = ((n / 100) % 10) as u32;
        if hundreds != 0 {
            char::from_digit(hundreds, 10).expect("Hundreds digit should be < 10")
        } else {
            C
        }
    };
    let tens_ones = padded_2::<C>(n % 100);
    return [most_significant_char, tens_ones[0], tens_ones[1]];
}

pub fn format_f32<const C: char, const BEFORE_COMMA: usize, const AFTER_COMMA: usize>(
    n: f32,
) -> [char; BEFORE_COMMA + AFTER_COMMA + 1] {
    // Split float
    let int_part = n as u64;
    let scale = (AFTER_COMMA).pow(10);
    let frac_part = n - (int_part as f32);
    let frac_part = (frac_part * (scale as f32)) as u64;

    // Pad the integer part (BEFORE_COMMA digits, fill with C)
    let mut int_chars = [C; BEFORE_COMMA];
    {
        let mut fill = true;
        let mut temp = int_part;
        let mut idx = BEFORE_COMMA;
        while idx > 0 {
            idx -= 1;
            let digit = (temp % 10) as u32;
            if !fill || digit != 0 {
                int_chars[idx] = char::from_digit(digit, 10).unwrap_or(C);
                fill = false;
            }
            temp /= 10;
            if temp == 0 {
                break;
            }
        }
    }

    // Pad fractional part (AFTER_COMMA digits, fill with C)
    let mut frac_chars: [char; AFTER_COMMA] = ['0'; AFTER_COMMA];
    {
        let mut temp = frac_part;
        let mut idx = AFTER_COMMA;
        while idx > 0 {
            idx -= 1;
            let digit = (temp % 10) as u32;
            frac_chars[idx] = char::from_digit(digit, 10).unwrap_or(C);
            temp /= 10;
            if temp == 0 {
                break;
            }
        }
    }

    // Build final [char; BEFORE_COMMA + AFTER_COMMA + 1]
    let mut out = ['\0'; BEFORE_COMMA + AFTER_COMMA + 1];

    for i in 0..BEFORE_COMMA {
        out[i] = int_chars[i];
    }
    out[BEFORE_COMMA] = '.';
    for j in (BEFORE_COMMA + 1)..AFTER_COMMA {
        out[j] = frac_chars[j - (BEFORE_COMMA + 1)];
    }

    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn test_time(d: Duration, exp: &str) {
        let mut utf8_65: [u8; 9] = [0u8; 9];
        let result_65 = show_duration(d);
        for i in 0..9 {
            utf8_65[i] = result_65[i] as u8;
        }
        assert_eq!(str::from_utf8(&utf8_65), Ok(exp));
    }

    #[test]
    fn format_time_test_65() {
        test_time(Duration::new(65, 0), "01:05.000");
    }

    #[test]
    fn format_time_test_605() {
        test_time(Duration::new(605, 1_000_000), "10:05.001");
    }

    #[test]
    fn format_time_test_6h5m2s() {
        test_time(Duration::new((6 * 60 + 5) * 60 + 2, 1_000_000), "006:05:02");
    }

    #[test]
    fn padded2_test() {
        for i in 0..=9 {
            assert_eq!(padded_2::<'0'>(i)[0], '0');
            assert_eq!(
                padded_2::<'0'>(i)[1],
                char::from_digit(i.try_into().unwrap(), 10).unwrap()
            );
        }

        assert_eq!(padded_2::<'0'>(10)[0], '1');
        assert_eq!(padded_2::<'0'>(10)[1], '0');
        assert_eq!(padded_2::<'0'>(100)[0], '0');
        assert_eq!(padded_2::<'0'>(100)[1], '0');
    }

    #[test]
    fn padded3_test() {
        for i in 0..=99 {
            assert_eq!(padded_3::<'0'>(i)[0], '0');
            assert_eq!(
                padded_3::<'0'>(i)[1],
                char::from_digit((i / 10).try_into().unwrap(), 10).unwrap()
            );
            assert_eq!(
                padded_3::<'0'>(i)[2],
                char::from_digit((i % 10).try_into().unwrap(), 10).unwrap()
            );
        }

        assert_eq!(padded_3::<'0'>(10)[0], '0');
        assert_eq!(padded_3::<'0'>(10)[1], '1');
        assert_eq!(padded_3::<'0'>(10)[2], '0');
        assert_eq!(padded_3::<'0'>(100)[0], '1');
        assert_eq!(padded_3::<'0'>(100)[1], '0');
        assert_eq!(padded_3::<'0'>(100)[2], '0');
        assert_eq!(padded_3::<'0'>(1000)[0], '0');
        assert_eq!(padded_3::<'0'>(1000)[1], '0');
        assert_eq!(padded_3::<'0'>(1000)[2], '0');
    }
}
