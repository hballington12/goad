use std::sync::{
    atomic::{AtomicBool, AtomicUsize},
    Mutex,
};

use anyhow::Result;

use numpy::npyffi::NPY_ITER_GLOBAL_FLAGS;
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
    pub orientations: Orientations,  // shared (maximum) orientations
    pub settings: Settings,          // runtime settings
    pub result: Option<Results>,     // mean result of the problems
    pub result_var: Option<Results>, // result variance
    pub max_orientations: usize,
    is_converged: bool,
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

        let result = None;
        let result_var = Some(Results::new_empty(&bins));

        Ok(Self {
            geom,
            orientations,
            settings,
            result,
            result_var,
            max_orientations,
            is_converged: false,
        })
    }

    pub fn run(&mut self) {
        let problem_base = Problem::new(Some(self.geom.clone()), Some(self.settings.clone()));
        let mut i = 0;
        while !self.is_converged && i < self.max_orientations {
            println!("Iteration {}", i);
            let (alpha, beta, gamma) = self.orientations.eulers[i];
            let mut problem = problem_base.clone();
            if let Err(err) = problem.run(Some(&Euler::new(alpha, beta, gamma))) {
                eprintln!("Error running problem (will skip this iteration): {}", err);
            }

            i += 1;
            if i == 1 {
                self.result = Some(problem.result);
            } else {
                let mm = self.result.take().unwrap(); // remove unwrap
                let ss = self.result_var.take().unwrap();
                let d = problem.result - mm.clone();
                self.result = Some(mm + d.clone() / (i as f32));
                self.result_var = Some(ss + d.clone() * d * (((i - 1) as f32) / (i as f32)));
            }

            let asymmetry = self.result.clone().unwrap().get_asymmetry().unwrap();

            let asymmetry_var = self
                .result_var
                .clone()
                .unwrap()
                .get_asymmetry()
                .unwrap_or_default();
            let asymmetry_sem = (asymmetry_var / ((i - 1) as f32)).sqrt();
            println!("asymmetry = {:?} + {:?}", asymmetry, asymmetry_sem);

            if i > 500 {
                self.is_converged = true;
            }
        }
        println!("done.")
    }
}

// replace later with a more meaningful condition to terminate
pub fn dummy_condition() -> bool {
    false
}
