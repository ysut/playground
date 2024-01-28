fn main() {
    let mut a = 0;
    let mut b = 1;
    let mut c = 0;
    for _ in 0..10 {
        c = a + b;
        a = b;
        b = c;
        println!("{}", c);
    }
}