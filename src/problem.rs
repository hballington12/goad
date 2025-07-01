//! Main simulation orchestrator for electromagnetic scattering calculations.
//!
//! This module contains the [`Problem`] struct which coordinates the complete GOAD
//! simulation pipeline from initial beam generation through near-field propagation
//! to far-field diffraction analysis. It manages the electromagnetic beam lifecycle,
//! handles surface interactions, and orchestrates result collection.
//!
//! The simulation process involves:
//! - Geometry initialization and beam discretization
//! - Near-field beam propagation with reflection/refraction
//! - Beam queue management and power tracking
//! - Far-field diffraction calculations
//! - Result aggregation and analysis
//!
//! # Simulation Pipeline
//!
//! 1. **Initialization**: Geometry setup, scaling, and centering
//! 2. **Illumination**: Initial beam creation and discretization
//! 3. **Propagation**: Near-field beam tracing through particle
//! 4. **Diffraction**: Far-field calculation from exit apertures
//! 5. **Analysis**: Mueller matrix computation and parameter extraction

use crate::{
    beam::{Beam, BeamPropagation, BeamType, BeamVariant},
    bins::{generate_bins, Scheme},
    field::Field,
    geom::{self, Face, Geom},
    helpers::draw_face,
    orientation, output,
    result::{self, Results},
    settings::{load_config, Settings},
};
use macroquad::prelude::*;
use nalgebra::{Complex, Matrix2, Point3, Vector3};
use ndarray::Array2;
use pyo3::prelude::*;
use rayon::prelude::*;

#[cfg(test)]
mod tests {

    use super::*;
    use nalgebra::Complex;

    #[test]
    fn cube_inside_ico() {
        let mut geom = Geom::from_file("./examples/data/cube_inside_ico.obj").unwrap();
        geom.shapes[0].refr_index = Complex {
            // modify the refractive index of the outer shape
            re: 2.0,
            im: 0.1,
        };
        geom.shapes[1].parent_id = Some(0); // set the parent of the second shape to be the first
        assert_ne!(geom.shapes[0].refr_index, geom.shapes[1].refr_index);
        assert_eq!(
            geom.shapes[0],
            geom.shapes[geom.shapes[1].parent_id.unwrap()]
        );

        let mut problem = Problem::new(geom, None);

        problem.propagate_next();
    }
}

/// An electromagnetic scattering simulation.
/// 
/// **Context**: Electromagnetic scattering simulations require coordinating geometry setup,
/// beam generation, near-field propagation through the particle, far-field diffraction
/// calculations, and result analysis. Each simulation tracks [`crate::beam::Beam`] segments as they interact
/// with particle surfaces, undergo reflection and refraction, and escape to produce the
/// scattered field.
/// 
/// **How it Works**: A [`Problem`] manages the simulation through multiple queues and state
/// tracking. The `beam_queue` holds beams awaiting near-field propagation, while `out_beam_queue`
/// and `ext_diff_beam_queue` manage beams ready for far-field diffraction. The `base_geom`
/// preserves the original [`crate::geom::Geom`] for reset operations, while `geom` is the working copy
/// that can be oriented and modified. Results accumulate in the [`crate::result::Results`] field, tracking
/// power flows, Mueller matrices, and derived parameters.
/// 
/// The simulation proceeds through initialization (geometry setup, centering, scaling),
/// illumination (initial beam creation), near-field solving (beam propagation and interaction),
/// and far-field solving ([`crate::diff`] to scattering angles).
#[pyclass]
#[derive(Debug, Clone)] // Added Default derive
pub struct Problem {
    pub base_geom: Geom,                // original [`crate::geom::Geom`] for reset operations
    pub geom: Geom,                     // working [`crate::geom::Geom`] copy for beam tracing
    pub beam_queue: Vec<Beam>,          // [`crate::beam::Beam`] objects awaiting near-field propagation
    pub out_beam_queue: Vec<Beam>,      // beams awaiting [`crate::diff`] calculations
    pub ext_diff_beam_queue: Vec<Beam>, // beams awaiting external diffraction
    pub settings: Settings,             // runtime [`crate::settings::Settings`]
    pub result: Results,                // accumulated [`crate::result::Results`] of the simulation
}

#[pymethods]
impl Problem {
    #[new]
    fn py_new(settings: Settings, geom: Option<Geom>) -> Self {
        match geom {
            // If a geometry is provided, use it
            Some(geom) => Problem::new(geom, Some(settings)),
            None => {
                // If no geometry is provided, load geometry from file
                let geom = geom::Geom::from_file(&settings.geom_name).unwrap();
                Problem::new(geom, Some(settings))
            }
        }
    }

    /// Setter function for the problem settings
    #[setter]
    pub fn set_settings(&mut self, settings: Settings) {
        self.settings = settings;
    }

    /// Getter function for the problem settings
    #[getter]
    pub fn get_settings(&self) -> Settings {
        self.settings.clone()
    }

    /// Getter function for the geometry
    #[getter]
    pub fn get_geom(&self) -> Geom {
        self.geom.clone()
    }

    /// Python interface for solving with discrete orientations.
    pub fn py_solve(&mut self) -> PyResult<()> {
        self.reset();

        // println!("{:#?}", self.settings);

        let euler = match self.settings.orientation.scheme {
            orientation::Scheme::Discrete { ref eulers } => &eulers[0].clone(),
            _ => {
                panic!("Python solve is only supperted for discrete orientation scheme")
            }
        };

        self.run(Some(euler));
        Ok(())
    }

    /// Python interface for printing power statistics.
    pub fn py_print_stats(&self) -> PyResult<()> {
        println!("{}", self.result.powers);
        Ok(())
    }

    // getter function to retrieve Python object containing the mueller matrix
    // convert the Array2 to a list of lists and return
    #[getter]
    pub fn get_mueller(&self) -> Vec<Vec<f32>> {
        collect_mueller(&self.result.mueller)
    }

    // getter function to retrieve Python object containing the asymmetry parameter
    #[getter]
    pub fn get_asymmetry(&self) -> Option<f32> {
        self.result.params.asymettry
    }

    // getter function to retrieve Python object containing the scattering cross section
    #[getter]
    pub fn get_scat_cross(&self) -> Option<f32> {
        self.result.params.scat_cross
    }

    // getter function to retrieve Python object containing the extinction cross section
    #[getter]
    pub fn get_ext_cross(&self) -> Option<f32> {
        self.result.params.ext_cross
    }

    // getter function to retrieve Python object containing the albedo
    #[getter]
    pub fn get_albedo(&self) -> Option<f32> {
        self.result.params.albedo
    }

    // getter function to retrieve Python object containing the 1d mueller matrix
    // convert the Array2 to a list of lists and return
    #[getter]
    pub fn get_mueller_1d(&self) -> Vec<Vec<f32>> {
        if let Some(mueller_1d) = &self.result.mueller_1d {
            let mut mueller_list = Vec::new();
            for row in mueller_1d.outer_iter() {
                let mut row_list = Vec::new();
                for val in row.iter() {
                    row_list.push(*val);
                }
                mueller_list.push(row_list);
            }
            mueller_list
        } else {
            Vec::new()
        }
    }

    /// getter function to retrieve Python object containing the theta values for the 1d mueller matrix
    #[getter]
    pub fn get_theta_1d(&self) -> Vec<f32> {
        if let Some(theta) = &self.result.bins_1d {
            theta.clone()
        } else {
            Vec::new()
        }
    }
}

impl Problem {
    /// Creates a new electromagnetic scattering problem.
    /// 
    /// **Context**: Setting up a scattering simulation requires initializing the [`crate::geom::Geom`]
    /// with appropriate material properties, configuring simulation parameters, and
    /// preparing the [`crate::result::Results`] structures. The problem needs both the original geometry
    /// for reset operations and a working copy for transformations.
    /// 
    /// **How it Works**: Applies initial geometry processing through [`init_geom()`], 
    /// generates angular bins for result collection using [`crate::bins::generate_bins`], and initializes empty result
    /// structures. Creates both `base_geom` (immutable reference) and `geom` (working copy)
    /// to enable multiple runs with different [`crate::orientation`] settings while preserving the
    /// original geometry state.
    /// 
    /// # Example
    /// ```rust
    /// let mut geom = crate::geom::Geom::from_file("./examples/data/hex2.obj").unwrap();
    /// geom.shapes[0].refr_index.re = 1.5;
    /// geom.shapes[0].refr_index.im = 0.0001;
    /// let mut problem = Problem::new(geom, None);
    /// ```
    pub fn new(mut geom: Geom, settings: Option<Settings>) -> Self {
        let settings = settings.unwrap_or_else(|| load_config().expect("Failed to load config"));
        init_geom(&settings, &mut geom);

        let bins = generate_bins(&settings.binning.scheme);
        let solution = Results::new_empty(&bins);

        let problem = Self {
            base_geom: geom.clone(),
            geom,
            beam_queue: vec![],
            out_beam_queue: vec![],
            ext_diff_beam_queue: vec![],
            settings,
            result: solution,
        };

        problem
    }

    /// Resets the problem state for a new simulation run.
    /// 
    /// **Context**: Multi-orientation simulations require running the same base
    /// problem with different particle orientations. Each run must start from
    /// a clean state with empty [`crate::beam::Beam`] queues and fresh [`crate::result::Results`] structures while
    /// preserving the original geometry configuration.
    /// 
    /// **How it Works**: Clears all beam queues, resets the [`crate::result::Results`] structure
    /// to empty while preserving the angular bins, and restores the working
    /// geometry from the preserved `base_geom` copy.
    /// 
    /// # Example
    /// ```rust
    /// println!("Resetting problem...");
    /// problem.reset();
    /// ```
    pub fn reset(&mut self) {
        self.beam_queue.clear();
        self.out_beam_queue.clear();
        self.ext_diff_beam_queue.clear();
        self.result = Results::new_empty(&self.result.bins);
        self.geom.clone_from(&self.base_geom);
    }

    /// Initializes geometry for simulation.
    /// 
    /// **Context**: Raw geometry from files needs processing before electromagnetic
    /// simulation - applying scaling transformations, distortions for roughness studies,
    /// centering for rotation operations, and normalizing to unit scale for numerical
    /// stability in the [`crate::field`] calculations.
    /// 
    /// **How it Works**: Applies optional anisotropic scaling from [`crate::settings::Settings`], applies
    /// optional [`crate::distortion::distort`] for surface roughness simulation, centers the geometry at
    /// the origin, then rescales to unit size. Stores the scaling factor for later
    /// conversion of results back to physical units.
    /// 
    /// # Example
    /// ```rust
    /// // From run() workflow
    /// self.init();
    /// match euler {
    ///     Some(euler) => {
    ///         self.orient(euler);
    ///     }
    ///     None => {}
    /// }
    /// self.illuminate();
    /// ```
    pub fn init(&mut self) {
        // Apply geometry scaling if set
        if let Some(scale) = &self.settings.geom_scale {
            self.geom.vector_scale(scale);
        }
        // Apply distortion if set
        if let Some(distortion) = self.settings.distortion {
            self.geom.distort(distortion, self.settings.seed);
        }
        self.geom.recentre();
        self.settings.scale = self.geom.rescale();
    }

    /// Creates the incident beam for the simulation.
    /// 
    /// **Context**: Electromagnetic scattering simulations require an incident
    /// electromagnetic wave to interact with the particle. The incident [`crate::beam::Beam`]
    /// must fully illuminate the particle geometry to capture all scattering
    /// interactions.
    /// 
    /// **How it Works**: Creates a plane wave [`crate::beam::Beam`] with cross-section large enough
    /// to encompass the particle bounds, propagating along the negative z-axis.
    /// The beam wavelength is scaled according to the geometry normalization
    /// factor to maintain consistent physics.
    /// 
    /// # Example
    /// ```rust
    /// println!("Illuminating problem...");
    /// problem.illuminate();
    /// ```
    pub fn illuminate(&mut self) {
        let scaled_wavelength = self.settings.wavelength * self.settings.scale;

        let beam = basic_initial_beam(
            &self.geom,
            scaled_wavelength,
            self.settings.medium_refr_index,
        );

        self.beam_queue.push(beam);
    }

    /// Creates a problem with a custom initial beam.
    /// 
    /// **Context**: Some simulations require specific incident [`crate::beam::Beam`] configurations
    /// rather than the standard plane wave illumination. This constructor allows
    /// direct specification of the initial electromagnetic [`crate::field::Field`].
    /// 
    /// **How it Works**: Sets up the [`Problem`] structure with the provided [`crate::beam::Beam`]
    /// already in the `beam_queue`, bypassing the standard [`illuminate()`] step.
    pub fn new_with_field(geom: Geom, beam: Beam) -> Self {
        let settings = load_config().expect("Failed to load config");

        let bins = generate_bins(&settings.binning.scheme);
        let solution = Results::new_empty(&bins);

        Self {
            base_geom: geom.clone(),
            geom,
            beam_queue: vec![beam],
            out_beam_queue: vec![],
            ext_diff_beam_queue: vec![],
            settings,
            result: solution,
        }
    }

    /// Computes far-field diffraction from outgoing beams.
    /// 
    /// **Context**: [`crate::beam::Beam`] objects exiting the particle geometry must be converted to
    /// far-field scattering amplitudes at specific angular bins. This involves
    /// Fourier transform calculations in [`crate::diff`] for each beam.
    /// 
    /// **How it Works**: Applies diffraction calculations to each beam in parallel,
    /// accumulates the resulting amplitude matrices, and adds them to the total
    /// far-field amplitude.
    fn diffract_outbeams(
        queue: &mut Vec<Beam>,
        bins: &[(f32, f32)],
        total_ampl_far_field: &mut [Matrix2<Complex<f32>>],
        fov_factor: Option<f32>,
    ) {
        let ampl_far_field = queue
            .par_iter()
            .map(|outbeam| outbeam.diffract(bins, fov_factor))
            .reduce(
                || vec![Matrix2::<Complex<f32>>::zeros(); bins.len()],
                |mut acc, local| {
                    for (a, l) in acc.iter_mut().zip(local) {
                        *a += l;
                    }
                    acc
                },
            );

        // TODO move this outside the function
        for (i, ampl) in ampl_far_field.iter().enumerate() {
            total_ampl_far_field[i] += ampl;
        }
    }

    /// Combines external diffraction and beam diffraction components.
    /// 
    /// **Context**: The total scattered field consists of contributions from
    /// both external diffraction (edge effects) and beam diffraction (transmitted
    /// and reflected beams). These must be coherently combined.
    /// 
    /// **How it Works**: Adds the beam amplitude matrix to the external
    /// diffraction amplitude matrix element-wise to produce the total amplitude.
    fn combine_far(&mut self) {
        self.result.ampl = self.result.ampl_ext.clone();
        for (i, ampl) in self.result.ampl.iter_mut().enumerate() {
            *ampl += self.result.ampl_beam[i];
        }
    }

    /// Solves external diffraction contribution to far-field.
    /// 
    /// **Context**: Geometric edges and surfaces create diffraction effects
    /// independent of transmitted beams. These external diffraction contributions
    /// must be calculated separately.
    /// 
    /// **How it Works**: Processes beams in the external diffraction queue
    /// without field-of-view truncation to capture all diffraction effects.
    pub fn solve_far_ext_diff(&mut self) {
        let fov_factor = None; // don't truncate by field of view for external diffraction
        Self::diffract_outbeams(
            &mut self.ext_diff_beam_queue,
            &self.result.bins,
            &mut self.result.ampl_ext,
            fov_factor,
        );
    }

    /// Solves beam diffraction contribution to far-field.
    /// 
    /// **Context**: Beams exiting the particle after reflection and refraction
    /// interactions contribute to the scattered field through diffraction.
    /// These contributions are typically the dominant scattering mechanism.
    /// 
    /// **How it Works**: Processes beams in the outgoing beam queue with
    /// optional field-of-view truncation for computational efficiency.
    pub fn solve_far_outbeams(&mut self) {
        let fov_factor = self.settings.fov_factor; // truncate by field of view for outbeams
        Self::diffract_outbeams(
            &mut self.out_beam_queue,
            &self.result.bins,
            &mut self.result.ampl_beam,
            fov_factor,
        );
    }
    /// Solves the complete far-field scattering problem.
    /// 
    /// **Context**: Far-field solving requires combining multiple scattering
    /// mechanisms - external diffraction and beam diffraction - to produce
    /// the total scattered amplitude.
    /// 
    /// **How it Works**: Sequentially solves external diffraction and beam
    /// diffraction, then combines the results coherently.
    pub fn solve_far(&mut self) {
        self.solve_far_ext_diff();
        self.solve_far_outbeams();
        self.combine_far();
    }

    /// Executes the complete electromagnetic scattering simulation.
    /// 
    /// **Context**: The scattering simulation involves two distinct phases - near-field
    /// propagation where beams interact with particle surfaces, and far-field diffraction
    /// where outgoing beams are converted to scattering angles. Each phase requires
    /// different algorithms and produces different components of the final solution.
    /// 
    /// **How it Works**: Calls solve_near() to propagate beams through the particle geometry,
    /// then solve_far() to compute diffraction from outgoing beams. Converts the resulting
    /// amplitude matrices to Mueller matrices for different scattering components (total,
    /// beam-only, external diffraction). Attempts to compute 1D angular distributions and
    /// derived scattering parameters when possible.
    /// 
    /// # Example
    /// ```rust
    /// // Traditional workflow
    /// problem.reset();
    /// problem.init();
    /// problem.illuminate();
    /// problem.solve();
    /// problem.writeup();
    /// ```
    pub fn solve(&mut self) {
        self.solve_near();
        self.solve_far();
        self.result.mueller = output::ampl_to_mueller(&self.result.bins, &self.result.ampl);
        self.result.mueller_beam =
            output::ampl_to_mueller(&self.result.bins, &self.result.ampl_beam);
        self.result.mueller_ext = output::ampl_to_mueller(&self.result.bins, &self.result.ampl_ext);

        // try compute 1d mueller
        match self.settings.binning.scheme {
            Scheme::Custom { .. } => {} // 1d mueller not supported for custom bins
            _ => {
                match result::try_mueller_to_1d(&self.result.bins, &self.result.mueller) {
                    Ok((theta, mueller_1d)) => {
                        // let _ = output::write_mueller_1d(&theta, &mueller_1d, "");
                        self.result.bins_1d = Some(theta);
                        self.result.mueller_1d = Some(mueller_1d);
                    }
                    Err(e) => {
                        println!("Failed to compute 1d mueller: {}", e);
                    }
                };
                match result::try_mueller_to_1d(&self.result.bins, &self.result.mueller_beam) {
                    Ok((theta, mueller_1d_beam)) => {
                        // let _ = output::write_mueller_1d(&theta, &mueller_1d_beam, "_beam");
                        self.result.bins_1d = Some(theta);
                        self.result.mueller_1d_beam = Some(mueller_1d_beam);
                    }
                    Err(e) => {
                        println!("Failed to compute 1d mueller (beam): {}", e);
                    }
                };
                match result::try_mueller_to_1d(&self.result.bins, &self.result.mueller_ext) {
                    Ok((theta, mueller_1d_ext)) => {
                        // let _ = output::write_mueller_1d(&theta, &mueller_1d_ext, "_ext");
                        self.result.bins_1d = Some(theta);
                        self.result.mueller_1d_ext = Some(mueller_1d_ext);
                    }
                    Err(e) => {
                        println!("Failed to compute 1d mueller (ext): {}", e);
                    }
                };
            }
        }
    }

    /// Attempts to compute 1D Mueller matrix from 2D results.
    /// 
    /// **Context**: When angular bins have appropriate symmetry, the 2D Mueller
    /// matrix can be reduced to a 1D angular distribution for easier analysis
    /// and comparison with other codes.
    /// 
    /// **How it Works**: Calls the result structure's try_mueller_to_1d method
    /// and handles any errors that occur during the reduction process.
    pub fn try_mueller_to_1d(&mut self) {
        match self.result.try_mueller_to_1d() {
            Ok(()) => {
                // println!("1d mueller computed successfully");
            }
            Err(e) => {
                println!("Failed to compute 1d mueller: {}", e);
            }
        }
    }

    /// Attempts to compute derived scattering parameters.
    /// 
    /// **Context**: The Mueller matrix contains all scattering information, but
    /// specific parameters like scattering cross-section, asymmetry parameter,
    /// and albedo are often needed for analysis and comparison.
    /// 
    /// **How it Works**: Calls the result structure's compute_params method
    /// using the simulation wavelength and handles any computation errors.
    pub fn try_params(&mut self) {
        match self.result.compute_params(self.settings.wavelength) {
            Ok(()) => {
                // println!("Params computed successfully");
            }
            Err(e) => {
                println!("Failed to compute params: {}", e);
            }
        }
    }

    /// Executes a complete scattering simulation with optional orientation.
    /// 
    /// **Context**: The standard workflow for electromagnetic scattering involves
    /// multiple sequential steps that must be performed in the correct order.
    /// This high-level interface handles the entire pipeline from geometry
    /// initialization through result analysis.
    /// 
    /// **How it Works**: Performs geometry initialization, applies optional rotation,
    /// creates the incident beam, solves both near-field and far-field problems,
    /// and computes derived results. This is the primary interface for single-
    /// orientation simulations and the building block for multi-orientation studies.
    /// 
    /// # Example
    /// ```rust
    /// // With orientation
    /// let mut problem = problem_base.clone();
    /// let euler = Euler::new(*a, *b, *g);
    /// problem.run(Some(&euler));
    /// 
    /// // Without orientation
    /// problem.run(None);
    /// ```
    pub fn run(&mut self, euler: Option<&orientation::Euler>) {
        // println!("Running problem...");
        // println!("{:#?}", self.settings);
        self.init();
        match euler {
            Some(euler) => {
                self.orient(euler);
            }
            None => {
                // No rotation
            }
        }
        self.illuminate();
        self.solve();
        self.try_mueller_to_1d();
        self.try_params();
    }

    /// Solves near-field beam propagation through the particle.
    /// 
    /// **Context**: Near-field simulation involves iteratively propagating beams
    /// through the particle geometry until either all beams exit the system or
    /// the remaining beam power falls below a convergence threshold.
    /// 
    /// **How it Works**: Processes beams from the queue using propagate_next()
    /// until the queue is empty or the power cutoff criterion is met. Tracks
    /// truncated power for beams discarded due to the cutoff.
    pub fn solve_near(&mut self) {
        loop {
            if self.beam_queue.len() == 0 {
                break;
            }

            let input_power = self.result.powers.input;
            let output_power = self.result.powers.output;

            if output_power / input_power > self.settings.cutoff {
                // add remaining power in beam queue to missing power due to cutoff
                self.result.powers.trnc_cop += self
                    .beam_queue
                    .iter()
                    .map(|beam| beam.power() / self.settings.scale.powi(2))
                    .sum::<f32>();
                break;
            }

            self.propagate_next();
        }
    }

    /// Writes simulation results to output files.
    /// 
    /// **Context**: Simulation results need to be saved in standard formats
    /// for analysis, visualization, and comparison with other codes or
    /// experimental data.
    /// 
    /// **How it Works**: Writes Mueller matrices for total, beam, and external
    /// diffraction components, plus summary result data to the configured
    /// output directory.
    pub fn writeup(&self) {
        let _ = output::write_mueller(
            &self.result.bins,
            &self.result.mueller,
            "",
            &self.settings.directory,
        );
        let _ = output::write_mueller(
            &self.result.bins,
            &self.result.mueller_beam,
            "_beam",
            &self.settings.directory,
        );
        let _ = output::write_mueller(
            &self.result.bins,
            &self.result.mueller_ext,
            "_ext",
            &self.settings.directory,
        );
        let _ = output::write_result(&self.result, &self.settings.directory);
    }

    /// Processes the next beam from the propagation queue.
    /// 
    /// **Context**: Near-field electromagnetic simulation works by tracking individual
    /// beam segments as they interact with particle surfaces. Each beam can split into
    /// multiple reflected and transmitted beams at interfaces, creating a tree of
    /// interactions that must be processed iteratively until convergence.
    /// 
    /// **How it Works**: Pops the highest-power beam from the queue, propagates it
    /// through the geometry to find surface intersections, and generates output beams
    /// from reflection/refraction. Tracks power flows for different beam types and
    /// adds new beams to appropriate queues. Returns the propagation result for
    /// visualization or analysis.
    /// 
    /// # Example
    /// ```rust
    /// // Interactive propagation
    /// loop {
    ///     if is_key_pressed(KeyCode::Enter) {
    ///         let next_propagation = problem.propagate_next().unwrap();
    ///         println!("number of beams in beam queue: {:?}", problem.beam_queue.len());
    ///     }
    ///     next_frame().await;
    /// }
    /// ```
    pub fn propagate_next(&mut self) -> Option<BeamPropagation> {
        // Try to pop the next beam from the queue
        let Some(mut beam) = self.beam_queue.pop() else {
            println!("No beams left to pop!");
            return None;
        };

        // Compute the outputs by propagating the beam
        let outputs = match &mut beam.type_ {
            BeamType::Default => self.propagate_default(&mut beam),
            BeamType::Initial => self.propagate_initial(&mut beam),
            _ => {
                println!("Unknown beam type, returning empty outputs.");
                Vec::new()
            }
        };

        self.result.powers.absorbed += beam.absorbed_power / self.settings.scale.powi(2);
        self.result.powers.trnc_clip +=
            (beam.clipping_area - beam.csa()) * beam.power() / self.settings.scale.powi(2);

        // Process each output beam
        for output in outputs.iter() {
            let output_power = output.power() / self.settings.scale.powi(2);
            match (&beam.type_, &output.type_) {
                (BeamType::Default, BeamType::Default) => self.insert_beam(output.clone()),
                (BeamType::Default, BeamType::OutGoing) => {
                    self.result.powers.output += output_power;
                    self.insert_outbeam(output.clone());
                }
                (BeamType::Initial, BeamType::Default) => {
                    self.result.powers.input += output_power;
                    self.insert_beam(output.clone());
                }
                (BeamType::Initial, BeamType::ExternalDiff) => {
                    self.result.powers.ext_diff += output_power;
                    self.ext_diff_beam_queue.push(output.clone());
                }
                _ => {}
            }
        }
        Some(BeamPropagation::new(beam, outputs))
    }

    /// Propagates an initial beam without threshold checks.
    /// 
    /// **Context**: Initial beams require different handling than internal beams
    /// since they haven't yet been subject to splitting and attenuation within
    /// the particle geometry.
    /// 
    /// **How it Works**: Directly calls the beam propagation method without
    /// applying power, area, or recursion threshold checks.
    fn propagate_initial(&mut self, beam: &mut Beam) -> Vec<Beam> {
        match beam.propagate(
            &mut self.geom,
            self.settings.medium_refr_index,
            self.settings.beam_area_threshold(),
        ) {
            Ok((outputs, ..)) => outputs,

            Err(_) => Vec::new(),
        }
    }

    /// Propagates a beam with threshold checks for computational efficiency.
    /// 
    /// **Context**: During beam propagation, many beams become small or weak
    /// and contribute negligibly to the final result. Threshold checks prevent
    /// wasting computation on these insignificant beams.
    /// 
    /// **How it Works**: Applies sequential checks for beam power, area, total
    /// internal reflection count, and recursion depth. Beams failing any check
    /// are terminated with their power added to the appropriate truncation category.
    fn propagate_default(&mut self, beam: &mut Beam) -> Vec<Beam> {
        // beam power is below threshold
        if beam.power() < self.settings.beam_power_threshold * self.settings.scale.powi(2) {
            self.result.powers.trnc_energy += beam.power() / self.settings.scale.powi(2);
            return Vec::new();
        }

        // beam area is below threshold
        if beam.face.data().area.unwrap() < self.settings.beam_area_threshold() {
            self.result.powers.trnc_area += beam.power() / self.settings.scale.powi(2);
            return Vec::new();
        }

        // total internal reflection considerations
        if beam.variant == Some(BeamVariant::Tir) {
            if beam.tir_count > self.settings.max_tir {
                self.result.powers.trnc_ref += beam.power() / self.settings.scale.powi(2);
                return Vec::new();
            } else {
                return self.propagate(beam);
            }
        }

        // beam recursion over the maximum
        if beam.rec_count > self.settings.max_rec {
            self.result.powers.trnc_rec += beam.power() / self.settings.scale.powi(2);
            return Vec::new();
        }

        // else, propagate the beam
        self.propagate(beam)
    }

    /// Core beam propagation with error handling.
    /// 
    /// **Context**: Beam propagation through geometry can fail due to numerical
    /// issues or geometric degeneracies. These failures must be handled gracefully
    /// without crashing the simulation.
    /// 
    /// **How it Works**: Calls the beam's propagate method and handles both
    /// successful propagation (returning output beams) and errors (adding
    /// beam power to the clipping error category).
    fn propagate(&mut self, beam: &mut Beam) -> Vec<Beam> {
        match beam.propagate(
            &mut self.geom,
            self.settings.medium_refr_index,
            self.settings.beam_area_threshold(),
        ) {
            Ok((outputs, area_power_loss)) => {
                self.result.powers.trnc_area += area_power_loss / self.settings.scale.powi(2);
                outputs
            }
            Err(_) => {
                self.result.powers.clip_err += beam.power() / self.settings.scale.powi(2);
                Vec::new()
            }
        }
    }

    /// Renders geometry and beam propagation for visualization.
    /// 
    /// **Context**: Interactive debugging and educational visualization require
    /// rendering both the particle geometry and the beam propagation paths.
    /// 
    /// **How it Works**: Draws all faces of all shapes in the geometry, then
    /// calls the propagation's draw method to render beam paths.
    pub fn draw_propagation(&self, propagation: &BeamPropagation) {
        // draw the geometry
        for shape in &self.geom.shapes {
            for face in &shape.faces {
                draw_face(face, GREEN, 4.0);
            }
        }
        propagation.draw();
    }

    /// Inserts a beam into the propagation queue in power-sorted order.
    /// 
    /// **Context**: Processing beams in order of decreasing power improves
    /// convergence by handling the most significant contributions first.
    /// This enables earlier termination when power cutoffs are reached.
    /// 
    /// **How it Works**: Uses binary search to find the insertion position
    /// that maintains ascending power order (since beams are processed by popping).
    pub fn insert_beam(&mut self, beam: Beam) {
        let pos = get_position_by_power(beam.power(), &self.beam_queue, true);
        self.beam_queue.insert(pos, beam);
    }

    /// Inserts an outgoing beam into the far-field queue in power-sorted order.
    /// 
    /// **Context**: Outgoing beams contribute to far-field scattering through
    /// diffraction calculations. Processing them in order of decreasing power
    /// allows for power-based truncation in the far-field solver.
    /// 
    /// **How it Works**: Uses binary search to find the insertion position
    /// that maintains descending power order for sequential processing.
    pub fn insert_outbeam(&mut self, beam: Beam) {
        let pos = get_position_by_power(beam.power(), &self.out_beam_queue, false);
        self.out_beam_queue.insert(pos, beam);
    }

    /// Applies a rotation to the particle geometry.
    /// 
    /// **Context**: Light scattering properties depend strongly on particle orientation
    /// relative to the incident beam. Systematic orientation studies require rotating
    /// the particle through different Euler angles while keeping the incident beam
    /// direction fixed.
    /// 
    /// **How it Works**: Applies the Euler rotation to the working geometry using
    /// the rotation convention specified in the settings. Requires the geometry
    /// to be centered at the origin for rotation around the particle center.
    /// 
    /// # Example
    /// ```rust
    /// pub fn orient(&mut self, euler: &orientation::Euler) {
    ///     if let Err(error) = self
    ///         .geom
    ///         .euler_rotate(euler, self.settings.orientation.euler_convention)
    ///     {
    ///         panic!("Error rotating geometry: {}", error);
    ///     }
    /// }
    /// ```
    pub fn orient(&mut self, euler: &orientation::Euler) {
        if let Err(error) = self
            .geom
            .euler_rotate(euler, self.settings.orientation.euler_convention)
        {
            panic!("Error rotating geometry: {}", error);
        }
    }
}

/// Converts a 2D Mueller matrix to nested vectors for Python interface.
/// 
/// **Context**: Python bindings require converting Rust ndarray structures
/// to Python-compatible nested lists for data exchange.
/// 
/// **How it Works**: Iterates through matrix rows and converts each row
/// to a vector, collecting all rows into a vector of vectors.
pub fn collect_mueller(array2: &Array2<f32>) -> Vec<Vec<f32>> {
    let mut mueller_list = Vec::new();
    for row in array2.outer_iter() {
        let mut row_list = Vec::new();
        for val in row.iter() {
            row_list.push(*val);
        }
        mueller_list.push(row_list);
    }
    mueller_list
}

/// Finds insertion position for power-sorted beam queue.
/// 
/// **Context**: Maintaining power-sorted beam queues requires efficient
/// insertion that preserves the sort order. Binary search provides
/// O(log n) insertion position finding.
/// 
/// **How it Works**: Uses binary search with configurable sort order
/// (ascending for pop-based processing, descending for sequential processing)
/// to find the correct insertion index.
fn get_position_by_power(value: f32, queue: &Vec<Beam>, ascending: bool) -> usize {
    queue
        .binary_search_by(|x| {
            let cmp = x
                .power()
                .partial_cmp(&value)
                .unwrap_or(std::cmp::Ordering::Equal);

            if ascending {
                cmp
            } else {
                cmp.reverse()
            }
        })
        .unwrap_or_else(|e| e)
}

/// Creates a plane wave incident beam that fully illuminates the particle.
/// 
/// **Context**: Electromagnetic scattering simulations require an incident wave
/// that completely illuminates the particle to capture all scattering interactions.
/// The beam must be large enough to encompass the particle bounds with appropriate
/// field amplitude and phase relationships.
/// 
/// **How it Works**: Creates a rectangular beam cross-section slightly larger than
/// the particle bounding box, positioned above the particle along the z-axis.
/// Sets up a plane wave field with linear polarization and propagates the phase
/// backwards to simulate incidence from z=0.
/// 
/// # Example
/// ```rust
/// let scaled_wavelength = self.settings.wavelength * self.settings.scale;
/// let beam = basic_initial_beam(
///     &self.geom,
///     scaled_wavelength,
///     self.settings.medium_refr_index,
/// );
/// self.beam_queue.push(beam);
/// ```
fn basic_initial_beam(geom: &Geom, wavelength: f32, medium_refractive_index: Complex<f32>) -> Beam {
    const FAC: f32 = 1.1; // scale factor to stretch beam to cover geometry
    let bounds = geom.bounds();
    let (min, max) = (bounds.0.map(|v| v * FAC), bounds.1.map(|v| v * FAC));

    let clip_vertices = vec![
        Point3::new(max[0], max[1], max[2]),
        Point3::new(max[0], min[1], max[2]),
        Point3::new(min[0], min[1], max[2]),
        Point3::new(min[0], max[1], max[2]),
    ];

    let mut clip = Face::new_simple(clip_vertices, None, None).unwrap();
    clip.data_mut().area = Some((max[0] - min[0]) * (max[1] - min[1]));
    let mut field = Field::new_identity(Vector3::x(), -Vector3::z()).unwrap();

    // propagate field backwards so its as if the beam comes from z=0
    let dist = bounds.1[2];
    let wavenumber = 2.0 * std::f32::consts::PI / wavelength;
    let arg = -dist * wavenumber * medium_refractive_index.re;
    field.ampl *= Complex::new(arg.cos(), arg.sin());

    let beam = Beam::new_from_field(
        clip,
        -Vector3::z(),
        medium_refractive_index,
        field,
        wavelength,
    );
    beam
}

/// Initializes geometry with material properties from settings.
/// 
/// **Context**: Geometry loaded from files contains only shape information.
/// Electromagnetic simulations require material properties (refractive indices)
/// to be assigned to each shape for reflection/refraction calculations.
/// 
/// **How it Works**: Assigns refractive indices from the settings array to
/// shapes in order, using the first value as default for any shapes beyond
/// the array length. Recenters the geometry as a final preparation step.
/// 
/// # Example
/// ```rust
/// let settings = load_config().expect("Failed to load config");
/// init_geom(&settings, &mut geom);
/// ```
pub fn init_geom(settings: &Settings, geom: &mut Geom) {
    for shape in geom.shapes.iter_mut() {
        shape.refr_index = settings.particle_refr_index[0]; // default refr index is first value
    }
    for (i, refr_index) in settings.particle_refr_index.iter().enumerate() {
        if i >= geom.shapes.len() {
            break;
        }
        geom.shapes[i].refr_index = *refr_index;
    }
    geom.recentre();
}
