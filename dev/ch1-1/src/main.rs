fn main() {
    let moon: f64 = 384400.0;
    let car: f64 = 80.0;
    let btrain: f64 = 300.0;
    println!("The days it takes to get to the moon by car is {}", moon / car / 24.0);
    println!("The days it takes to get to the moon by bullet train is {}", moon / btrain / 24.0);
}
