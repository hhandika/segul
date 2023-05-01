use segul::cli;

#[cfg(not(tarpaulin_include))]
fn main() {
    cli::parse_cli();
}
