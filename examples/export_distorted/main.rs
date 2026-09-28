//! Export a single distorted realisation of a geometry as OBJ.
//!
//! The reviewer-requested irregular particle must be one fixed realisation
//! shared by every method, whereas a normal GOAD run draws a fresh
//! realisation per orientation. This tool applies `Geom::distort` once, with
//! an explicit seed, and writes the result.
//!
//! Usage:
//!     cargo run --release --example export_distorted -- <in.obj> <sigma> <seed> <out.obj>

use goad::geom::Geom;
use nalgebra::Complex;

fn main() -> anyhow::Result<()> {
    let args: Vec<String> = std::env::args().collect();
    let (input, sigma, seed, output) = (
        &args[1],
        args[2].parse::<f32>()?,
        args[3].parse::<u64>()?,
        &args[4],
    );
    let mut geoms = Geom::load(input, vec![Complex::new(1.31, 0.0)])?;
    let geom = &mut geoms[0];
    geom.distort(sigma, Some(seed))?;
    geom.write_obj(output)?;
    println!("wrote {output} (sigma {sigma}, seed {seed})");
    Ok(())
}
