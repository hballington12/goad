use std::sync::{
    atomic::{AtomicBool, AtomicUsize},
    Mutex,
};

use anyhow::Result;

use rayon::iter::{IntoParallelIterator, IntoParallelRefIterator, ParallelIterator};

use crate::{
    geom::Geom,
    orientation::{Euler, Orientations, Scheme},
    problem::{init_geom, Problem},
    result::{GOComponent, Results},
    settings::Settings,
};

pub struct Convergence {
    pub geom: Geom,
    pub orientations: Orientations, // shared (maximum) orientations
    pub settings: Settings,         // runtime settings
    pub result: Results,            // mean result of the problems
    pub result_var: Results,        // result variance
}

impl Convergence {
    /// Creates a new `MultiProblem` from optional `Geom` and `Settings`.
    /// If settings not provided, loads from config file.
    /// If geom not provided, loads from file using settings.geom_name.
    pub fn new(
        geom: Option<Geom>,
        settings: Option<Settings>,
        max_orientations: usize,
    ) -> anyhow::Result<Self> {
        let settings = settings
            .unwrap_or_else(|| crate::settings::load_config().expect("Failed to load config"));
        let mut geom = match geom {
            Some(g) => g,
            None => Geom::from_file(&settings.geom_name).map_err(|e| {
                anyhow::anyhow!(
                    "Failed to load geometry file '{}': {}\n\
                    Hint: This may be caused by degenerate faces (zero cross product), \
                    faces that are too small, or non-planar geometry. \
                    Please check and fix the geometry file.",
                    settings.geom_name,
                    e
                )
            })?,
        };

        init_geom(&settings, &mut geom);

        let orientation_scheme = Scheme::Uniform {
            num_orients: max_orientations,
        };
        let orientations = Orientations::generate(&orientation_scheme, settings.seed);
        let bins = &settings.binning.scheme.generate();

        let result = Results::new_empty(&bins);
        let result_var = Results::new_empty(&bins);

        Ok(Self {
            geom,
            orientations,
            settings,
            result,
            result_var,
        })
    }

    // 1. Update function signature to return Result (Error fix E0277)
    pub fn run(&mut self) {
        // let problem_base = Problem::new(Some(self.geom.clone()), Some(self.settings.clone()));

        // let bins = self.result.bins();

        // self.result = self
        //     .orientations
        //     .eulers
        //     .par_iter()
        //     .try_fold(
        //         // 2. Identity must return the raw type 'Results', not 'Ok(Results)'
        //         || Results::new_empty(&bins),
        //         |accum, (a, b, g)| {
        //             let mut problem = problem_base.clone();
        //             let euler = Euler::new(*a, *b, *g);

        //             // This '?' works because the closure returns Result<Results, E>
        //             let b: () = problem.run(Some(&euler))?;

        //             Ok(accum + problem.result)
        //         },
        //     )
        //     .try_reduce(
        //         // 3. Identity for reduce must also return 'Results'
        //         || Results::new_empty(&bins),
        //         // 4. Inputs 'a' and 'b' are already unwrapped 'Results'.
        //         //    No '?' needed on a or b. Return 'Ok' to satisfy signature.
        //         |a, b| Ok(a + b),
        //     )?; // 5. The final '?' propagates the error out of 'fn run'

        let bytes = 0..22_u8;
        let sum = bytes
            .into_par_iter()
            .try_fold(|| 0_u32, |a: u32, b: u8| a.checked_add(b as u32))
            .try_reduce(|| 0, u32::checked_add); // Ok(())

        let problem_base = Problem::new(Some(self.geom.clone()), Some(self.settings.clone()));
        let bins = &self.result.bins();
        let null_results = Results::new_empty(bins);

        let result = self
            .orientations
            .eulers
            .par_iter()
            .try_fold(
                || null_results.clone(),
                |mut accum, (a, b, g)| {
                    let mut problem = problem_base.clone();
                    let euler = Euler::new(*a, *b, *g);

                    if let Err(err) = problem.run(Some(&euler)) {
                        eprintln!("Error running problem (will skip this iteration): {}", err);
                    }
                    println!(
                        "asymmetry individual is {:?}, current is {:?}",
                        problem.result.params.asymmetry(&GOComponent::Total),
                        accum.params.asymmetry(&GOComponent::Total)
                    );

                    accum = accum + problem.result;
                    // println!("results areis {:?}", accum);

                    if dummy_condition() {
                        Err(accum)
                    } else {
                        Ok(accum)
                    }
                },
            )
            .try_reduce(
                || null_results.clone(),
                |mut accum, item| {
                    // stop on global count > 50
                    println!(
                        "accumulating: {:?} | {:?}",
                        accum.params.asymmetry(&GOComponent::Total),
                        item.params.asymmetry(&GOComponent::Total)
                    );
                    accum = accum + item;
                    if dummy_condition() {
                        Err(accum)
                    } else {
                        Ok(accum)
                    }
                },
            );
        if let Ok(result) = result {
            println!(
                "finished: {:?}",
                result.params.asymmetry(&GOComponent::Total)
            );
        }
    }
}

// replace later with a more meaningful condition to terminate
pub fn dummy_condition() -> bool {
    false
}
