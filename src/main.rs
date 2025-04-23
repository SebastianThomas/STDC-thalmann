const depth_incr: f64 = 3.0;

#[no_mangle]
pub fn thalmann(loading: Vec<f64>, maxDepth: f64, beta: Vec<Vec<f64>>) -> bool {
    let n: u64 = beta::len();
    let M: Vec<Vec<f64>> = Vec::new();
    for i in 0..n {
        for depth in 0..(maxDepth) {
            M[i][depth] = f(beta[i], depth);
        }
    }
    
    let mut Ddelta = ceil(maxDepth / depth_incr);
    let p = loading;
    // Reach first Stop
    while checkNoMpttExceeded(n, Ddelta, p, M) {
        Ddelta -= 1;
    }

    loop {
        let mut max = -1;
        for i in 0..n {
            max = max(max, p[i] - M[i][next_depth]);
        }
        if max >= 0 {
            Ddelta -= 1;
        }
    }
}

#[no_mangle]
pub fn checkNoMpttExceeded(n: u64, dep_delta: u64, p: Vec<f64>, M: Vec<Vec<f64>>) -> bool {
    assert(p::len() == n)
    assert(M::len() == n);
    for i in 0..n {
        if (p[i] - M[i][dep_delta] > 0) {
            return false;
        }
    }
    return true;
}

// If you use `main()`, declare it as `pub` to see it in the output:
pub fn main() { 
    thalmann();
}
