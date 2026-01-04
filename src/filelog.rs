//! File-based logging to avoid conflicts with indicatif progress bars.
//!
//! Usage:
//! ```ignore
//! use goad::filelog;
//!
//! // Initialize with output directory
//! filelog::init(&output_dir);
//! // Logs: "Logs will be written to: /path/to/output/goad.log"
//!
//! // Now all log macros write to file instead of stderr
//! log::info!("This goes to the log file");
//! ```

use std::fs::File;
use std::io::{self, Stderr, Write};
use std::path::Path;
use std::sync::{Mutex, OnceLock};

use log::{Level, LevelFilter, Log, Metadata, Record};

enum LogOutput {
    Stderr(Stderr),
    File(File),
}

impl Write for LogOutput {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        match self {
            LogOutput::Stderr(s) => s.write(buf),
            LogOutput::File(f) => f.write(buf),
        }
    }

    fn flush(&mut self) -> io::Result<()> {
        match self {
            LogOutput::Stderr(s) => s.flush(),
            LogOutput::File(f) => f.flush(),
        }
    }
}

struct SwitchableLogger {
    output: Mutex<LogOutput>,
}

impl Log for SwitchableLogger {
    fn enabled(&self, metadata: &Metadata) -> bool {
        metadata.level() <= Level::Info
    }

    fn log(&self, record: &Record) {
        if self.enabled(record.metadata()) {
            if let Ok(mut output) = self.output.lock() {
                let timestamp = chrono::Local::now().format("%Y-%m-%d %H:%M:%S");
                let _ = writeln!(
                    output,
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
        if let Ok(mut output) = self.output.lock() {
            let _ = output.flush();
        }
    }
}

static LOGGER: OnceLock<SwitchableLogger> = OnceLock::new();

/// Initialize file-based logging to the specified directory.
///
/// First logs to stderr to announce the log file location, then switches
/// to file-based logging. This allows the "Logs will be written to:" message
/// to appear in the terminal before progress bars start.
///
/// Returns the path to the log file.
pub fn init(output_dir: &Path) -> std::io::Result<std::path::PathBuf> {
    // Ensure directory exists
    std::fs::create_dir_all(output_dir)?;

    let log_path = output_dir.join("goad.log");

    // Initialize logger with stderr output first (if not already initialized)
    let logger = LOGGER.get_or_init(|| SwitchableLogger {
        output: Mutex::new(LogOutput::Stderr(io::stderr())),
    });

    // Set as global logger (ignore error if already set)
    let _ = log::set_logger(logger);
    log::set_max_level(LevelFilter::Info);

    // Log the path to stderr
    log::info!("Logs will be written to: {}", log_path.display());

    // Now switch to file output
    let file = File::create(&log_path)?;
    if let Ok(mut output) = logger.output.lock() {
        *output = LogOutput::File(file);
    }

    Ok(log_path)
}
