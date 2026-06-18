fn main() {
    if let Err(err) = gatk_faster_rnaseq_rust::interval_tools::run(std::env::args().collect()) {
        eprintln!("error: {}", err.message);
        std::process::exit(err.exit_code);
    }
}
