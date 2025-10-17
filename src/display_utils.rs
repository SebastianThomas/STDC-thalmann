use core::time::Duration;

const SECS_PER_HOUR: u64 = 60 * 60;
const SECS_PER_MIN: u64 = 60;
const ZERO_CHAR: char = char::from_digit(0, 10).unwrap();

pub fn show_duration(d: Duration) -> [char; 9] {
    let secs = d.as_secs();
    if secs > SECS_PER_HOUR {
        let hours = padded_3(secs / SECS_PER_HOUR);
        let mins = padded_2((secs % SECS_PER_HOUR) / SECS_PER_MIN);
        let secs = padded_2(secs % SECS_PER_MIN);
        return [
            hours[0], hours[1], hours[2], ':', mins[0], mins[1], ':', secs[0], secs[1],
        ];
    }
    let mins = padded_2((secs % SECS_PER_HOUR) / SECS_PER_MIN);
    let secs = padded_2(secs % SECS_PER_MIN);
    let millis = padded_3(d.subsec_millis() as u64);
    return [
        ':', mins[0], mins[1], ':', secs[0], secs[1], millis[0], millis[1], millis[2],
    ];
}

const fn padded_2(n: u64) -> [char; 2] {
    let most_significant_char: char = if n < 10 {
        ZERO_CHAR
    } else {
        let tens: u32 = ((n / 10) % 10) as u32;
        char::from_digit(tens, 10).expect("Tens digit should be < 10")
    };
    return [
        most_significant_char,
        char::from_digit((n % 10) as u32, 10).unwrap(),
    ];
}

const fn padded_3(n: u64) -> [char; 3] {
    let most_significant_char: char = {
        let hundreds: u32 = ((n / 100) % 10) as u32;
        char::from_digit(hundreds, 10).expect("Hundreds digit should be < 10")
    };
    let tens_ones = padded_2(n % 100);
    return [most_significant_char, tens_ones[0], tens_ones[1]];
}
