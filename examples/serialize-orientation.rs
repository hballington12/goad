use std::{fs::File, io::Write};

use pbt::orientation::{Config, OrientationScheme};

fn main() -> std::io::Result<()> {
    // Create an example OrientationScheme
    let orientation_scheme = OrientationScheme::Uniform { num_orients: 100 };

    // Create a Config struct with the orientation scheme
    let config = Config {
        orientation: orientation_scheme,
    };

    // Serialize the Config struct to a TOML string
    let toml_string = toml::to_string(&config).unwrap();

    // Write the TOML string to a file
    let mut file = File::create("config.toml")?;
    file.write_all(toml_string.as_bytes())?;

    println!("Configuration written to config.toml");

    Ok(())
}
