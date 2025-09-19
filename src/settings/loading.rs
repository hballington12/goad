use anyhow::Result;
use config::{Config, Environment, File};
use std::env;
use std::path::PathBuf;

use super::{cli, validation, Settings};

pub fn load_default_config() -> Result<Settings> {
    let goad_dir = retrieve_project_root()?;
    let default_config_file = goad_dir.join("config/default.toml");

    let settings: Config = Config::builder()
        .add_source(File::from(default_config_file).required(true))
        .build()
        .unwrap_or_else(|err| {
            eprintln!("Error loading configuration: {}", err);
            std::process::exit(1);
        });

    let mut config: Settings = settings.try_deserialize().unwrap_or_else(|err| {
        eprintln!("Error deserializing configuration: {}", err);
        std::process::exit(1);
    });

    validation::validate_config(&mut config);

    Ok(config)
}

pub fn load_config() -> Result<Settings> {
    load_config_with_cli(true)
}

pub fn load_config_with_cli(apply_cli_updates: bool) -> Result<Settings> {
    let config_file = get_config_file()?;

    let settings: Config = Config::builder()
        .add_source(File::from(config_file).required(true))
        .add_source(Environment::with_prefix("goad"))
        .build()
        .unwrap_or_else(|err| {
            eprintln!("Error loading configuration: {}", err);
            std::process::exit(1);
        });

    let mut config: Settings = settings.try_deserialize().unwrap_or_else(|err| {
        eprintln!("Error deserializing configuration: {}", err);
        std::process::exit(1);
    });

    if apply_cli_updates {
        cli::update_settings_from_cli(&mut config);
    }

    validation::validate_config(&mut config);

    Ok(config)
}

fn get_config_file() -> Result<PathBuf, anyhow::Error> {
    let current_dir_config = std::env::current_dir()
        .map(|dir| dir.join("local.toml"))
        .unwrap();
    let config_file = if current_dir_config.exists() {
        current_dir_config
    } else {
        let goad_dir = retrieve_project_root()?;
        let default_config_file = goad_dir.join("config/default.toml");
        let local_config = goad_dir.join("config/local.toml");

        if local_config.exists() {
            println!("Using local configuration: {:?}", local_config);
            local_config
        } else {
            println!("Using default configuration: {:?}", default_config_file);
            default_config_file
        }
    };
    Ok(config_file)
}

fn retrieve_project_root() -> Result<std::path::PathBuf> {
    if let Ok(manifest_dir) = env::var("CARGO_MANIFEST_DIR") {
        Ok(std::path::PathBuf::from(manifest_dir))
    } else if let Ok(path) = env::var("GOAD_ROOT_DIR") {
        Ok(std::path::PathBuf::from(path))
    } else {
        let exe_path = env::current_exe().expect("Failed to get current executable path");
        let mut current_dir = exe_path
            .parent()
            .expect("Failed to get executable directory")
            .to_path_buf();
        let mut found = false;

        while !found && current_dir.parent().is_some() {
            if current_dir.join("config").is_dir() {
                found = true;
            } else {
                current_dir = current_dir.parent().unwrap().to_path_buf();
            }
        }

        if found {
            Ok(current_dir)
        } else {
            Err(anyhow::anyhow!("Could not find project root directory"))
        }
    }
}
