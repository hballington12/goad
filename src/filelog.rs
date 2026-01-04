//! File-based logging to avoid conflicts with indicatif progress bars.
//!
//! Usage:
//! ```ignore
//! use goad::logging;
//!
//! // Initialize with output directory
//! logging::init(&output_dir);
//! // Prints: "Logs will be written to: /path/to/output/goad.log"
//!
//! // Now all log macros write to file instead of stderr
//! log::info!("This goes to the log file");
//! ```

use std::fs::File;
use std::io::Write;
use std::path::Path;
use std::sync::Mutex;

use log::{Level, LevelFilter, Log, Metadata, Record};

struct FileLogger {
    file: Mutex<File>,
}

impl Log for FileLogger {
    fn enabled(&self, metadata: &Metadata) -> bool {
        metadata.level() <= Level::Info
    }

    fn log(&self, record: &Record) {
        if self.enabled(record.metadata()) {
            if let Ok(mut file) = self.file.lock() {
                let timestamp = chrono::Local::now().format("%Y-%m-%d %H:%M:%S");
                let _ = writeln!(
                    file,
                    "[{}] [{}] [{}] {}",
                    timestamp,
                    record.level(),
                    record.target(),
                    record.args()
                );
            }
        }
    }

    fn flush(&self) {
        if let Ok(mut file) = self.file.lock() {
            let _ = file.flush();
        }
    }
}

/// Initialize file-based logging to the specified directory.
///
/// Creates a `goad.log` file in the output directory and redirects
/// all log output there. Prints the log file location to stdout
/// before any progress bars start.
///
/// Returns the path to the log file.
pub fn init(output_dir: &Path) -> std::io::Result<std::path::PathBuf> {
    // Ensure directory exists
    std::fs::create_dir_all(output_dir)?;

    let log_path = output_dir.join("goad.log");

    // Print location before progress bars start
    println!("Logs: {}", log_path.display());

    let file = File::create(&log_path)?;

    let logger = FileLogger {
        file: Mutex::new(file),
    };

    // Set as global logger (ignore error if already set)
    let _ = log::set_boxed_logger(Box::new(logger));
    log::set_max_level(LevelFilter::Info);

    Ok(log_path)
}

/// Initialize logging to stderr (default behavior, for CLI tools without progress bars).
pub fn init_stderr() {
    env_logger::init();
}
