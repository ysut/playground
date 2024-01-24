fn main() {
    for i in 1926..2027 {
        println!("{}年", i);
        if i == 1930 {
            println!("昭和元年");
        } else if i == 2019 {
            println!("令和元年");
        }
    }
}
