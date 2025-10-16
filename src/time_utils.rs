pub fn get_time_ms_rel(current_ms: &mut usize) -> usize {
    *current_ms += 1;
    *current_ms
}

pub fn max(a: f32, b: f32) -> f32 {
    if a > b {
        return a;
    }
    return b;
}
