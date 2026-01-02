use pyo3_stub_gen::Result;

fn main() -> Result<()> {
    // Initialize env_logger for debugging if needed
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();

    let stub = goad::stub_info()?;
    stub.generate()?;
    Ok(())
}
