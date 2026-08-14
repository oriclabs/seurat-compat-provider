//! GPL-3.0-only process provider for strict Seurat numeric boundaries.
//!
//! This executable is deliberately separate from BioLang. It invokes a pinned
//! R package environment and exchanges only generic BLMATF64 files with the
//! MIT-licensed caller; it is never linked into `bl`.

use std::env;
use std::fs;
use std::io::Write;
use std::path::PathBuf;
use std::process::{Command, ExitCode, Stdio};
use std::time::{SystemTime, UNIX_EPOCH};

const VERSION: &str = env!("CARGO_PKG_VERSION");
const PROVIDER_R: &str = include_str!("../r/provider.R");

fn usage() {
    eprintln!(
        "bl-seurat-provider {VERSION} (GPL-3.0-only)\n\n\
         Usage:\n\
           bl-seurat-provider doctor\n\
           bl-seurat-provider self-test\n\
           bl-seurat-provider cca LEFT.f64 RIGHT.f64 OUTPUT_DIR [DIMS] [SEED] [MAX_FEATURES]\n\
           bl-seurat-provider pca INPUT.f64 OUTPUT_DIR [COMPONENTS] [SEED] [CENTER]\n\n\
         Environment:\n\
           BIOLANG_RSCRIPT                 Rscript executable (default: Rscript)\n\
           BL_SEURAT_ALLOW_VERSION_MISMATCH=true  allow unpinned validation\n\n\
         Required R packages: Seurat 5.5.1, irlba 2.3.7, RcppAnnoy 0.0.23"
    );
}

fn temporary_script() -> Result<PathBuf, String> {
    let nonce = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map_err(|error| error.to_string())?
        .as_nanos();
    let path = env::temp_dir().join(format!(
        "bl-seurat-provider-{}-{nonce}.R",
        std::process::id()
    ));
    fs::write(&path, PROVIDER_R)
        .map_err(|error| format!("cannot create temporary provider script: {error}"))?;
    Ok(path)
}

fn run_provider(arguments: &[String]) -> Result<(), String> {
    let script = temporary_script()?;
    let rscript = env::var("BIOLANG_RSCRIPT").unwrap_or_else(|_| "Rscript".to_string());
    let status = Command::new(&rscript)
        .arg("--vanilla")
        .arg(&script)
        .args(arguments)
        .env("BL_SEURAT_PROVIDER_VERSION", VERSION)
        .stdin(Stdio::null())
        .stdout(Stdio::inherit())
        .stderr(Stdio::inherit())
        .status()
        .map_err(|error| format!("cannot start {rscript}: {error}"));
    let _ = fs::remove_file(&script);
    let status = status?;
    if status.success() {
        Ok(())
    } else {
        Err(format!(
            "provider R process exited with {}",
            status
                .code()
                .map_or_else(|| "no status".to_string(), |v| v.to_string())
        ))
    }
}

fn write_smoke_matrix(path: &PathBuf, batch_shift: f64) -> Result<(), String> {
    let rows = 12_u64;
    let columns = 20_u64;
    let mut file = fs::File::create(path).map_err(|error| error.to_string())?;
    file.write_all(b"BLMATF64")
        .and_then(|_| file.write_all(&rows.to_le_bytes()))
        .and_then(|_| file.write_all(&columns.to_le_bytes()))
        .map_err(|error| error.to_string())?;
    for row in 0..rows {
        for column in 0..columns {
            let group = if row < rows / 2 { 2.0 } else { -2.0 };
            let signal = if column < 8 { group } else { -group * 0.5 };
            let value = signal + batch_shift + ((row * 17 + column * 11 + 3) % 13) as f64 / 20.0;
            file.write_all(&value.to_le_bytes())
                .map_err(|error| error.to_string())?;
        }
    }
    Ok(())
}

fn self_test() -> Result<(), String> {
    let nonce = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map_err(|error| error.to_string())?
        .as_nanos();
    let directory = env::temp_dir().join(format!(
        "bl-seurat-provider-self-test-{}-{nonce}",
        std::process::id()
    ));
    fs::create_dir(&directory).map_err(|error| error.to_string())?;
    let left = directory.join("left.f64");
    let right = directory.join("right.f64");
    let cca_output = directory.join("cca");
    let pca_output = directory.join("pca");
    let result = (|| {
        write_smoke_matrix(&left, 0.0)?;
        write_smoke_matrix(&right, 0.4)?;
        run_provider(&[
            "cca".to_string(),
            left.to_string_lossy().into_owned(),
            right.to_string_lossy().into_owned(),
            cca_output.to_string_lossy().into_owned(),
            "3".to_string(),
            "42".to_string(),
            "12".to_string(),
        ])?;
        run_provider(&[
            "pca".to_string(),
            left.to_string_lossy().into_owned(),
            pca_output.to_string_lossy().into_owned(),
            "3".to_string(),
            "42".to_string(),
            "false".to_string(),
        ])?;
        for path in [
            cca_output.join("left-embedding.f64"),
            cca_output.join("right-embedding.f64"),
            cca_output.join("weight-reduction.f64"),
            cca_output.join("filter-features.csv"),
            pca_output.join("scores.f64"),
            pca_output.join("loadings.f64"),
        ] {
            if fs::metadata(&path).map(|value| value.len()).unwrap_or(0) == 0 {
                return Err(format!(
                    "self-test output missing or empty: {}",
                    path.display()
                ));
            }
        }
        Ok(())
    })();
    let _ = fs::remove_dir_all(&directory);
    result?;
    println!("BL_SEURAT_PROVIDER_SELF_TEST_OK");
    Ok(())
}

fn validate_command(arguments: &[String]) -> Result<(), String> {
    match arguments.first().map(String::as_str) {
        Some("doctor") if arguments.len() == 1 => Ok(()),
        Some("self-test") if arguments.len() == 1 => Ok(()),
        Some("cca") if (4..=7).contains(&arguments.len()) => Ok(()),
        Some("pca") if (3..=6).contains(&arguments.len()) => Ok(()),
        Some("--version" | "-V") if arguments.len() == 1 => Ok(()),
        Some("--help" | "-h") if arguments.len() == 1 => Ok(()),
        _ => Err("invalid arguments".to_string()),
    }
}

fn main() -> ExitCode {
    let arguments: Vec<String> = env::args().skip(1).collect();
    if validate_command(&arguments).is_err() {
        usage();
        return ExitCode::from(2);
    }
    match arguments.first().map(String::as_str) {
        Some("--version" | "-V") => {
            println!("bl-seurat-provider {VERSION}");
            ExitCode::SUCCESS
        }
        Some("--help" | "-h") => {
            usage();
            ExitCode::SUCCESS
        }
        Some("self-test") => match self_test() {
            Ok(()) => ExitCode::SUCCESS,
            Err(error) => {
                eprintln!("bl-seurat-provider: {error}");
                ExitCode::FAILURE
            }
        },
        _ => match run_provider(&arguments) {
            Ok(()) => ExitCode::SUCCESS,
            Err(error) => {
                eprintln!("bl-seurat-provider: {error}");
                ExitCode::FAILURE
            }
        },
    }
}

#[cfg(test)]
mod tests {
    use super::validate_command;

    fn args(values: &[&str]) -> Vec<String> {
        values.iter().map(|value| (*value).to_string()).collect()
    }

    #[test]
    fn accepts_documented_commands() {
        assert!(validate_command(&args(&["doctor"])).is_ok());
        assert!(validate_command(&args(&["self-test"])).is_ok());
        assert!(validate_command(&args(&["cca", "a", "b", "out"])).is_ok());
        assert!(validate_command(&args(&["pca", "a", "out", "50", "42", "false"])).is_ok());
    }

    #[test]
    fn rejects_incomplete_commands() {
        assert!(validate_command(&args(&[])).is_err());
        assert!(validate_command(&args(&["cca", "a"])).is_err());
        assert!(validate_command(&args(&["pca", "a"])).is_err());
        assert!(validate_command(&args(&["unknown"])).is_err());
    }
}
