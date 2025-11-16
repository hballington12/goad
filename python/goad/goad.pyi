"""Type stubs for goad.goad module"""

class Euler:
    """
    Euler angles representation. All angles in degrees.

    Attributes:
        `alpha`: alpha Euler angle
        `beta`: beta Euler angle
        `gamma`: gamma Euler angle
    """

    alpha: float
    beta: float
    gamma: float

    def __init__(self, alpha: float, beta: float, gamma: float) -> None:
        """
        Create a new Euler angle.

        Args:
            `alpha`: alpha Euler angle in degrees
            `beta`: beta Euler angle in degrees
            `gamma`: gamma Euler angle in degrees
        """
        ...

    def __repr__(self) -> str:
        """Return string representation of Euler angles."""
        ...

class EulerConvention:
    """Euler angle convention for rotation order.

    Specifies the order of rotations when applying Euler angles.
    Conventions are specified as three-letter strings (case-insensitive).

    Valid conventions:
        Proper Euler angles (symmetric):
            - XZX, XYX, YXY, YZY, ZYZ, ZXZ

        Tait-Bryan angles (asymmetric):
            - XZY, XYZ, YXZ, YZX, ZYX, ZXY

    Example:
        convention = EulerConvention('ZYZ')
        convention = EulerConvention('XYZ')
    """

    def __init__(self, order: str) -> None:
        """
        Create a new Euler angle convention.

        Args:
            order: Three-letter rotation order.
                   Valid options: `XZX`, `XYX`, `YXY`, `YZY`, `ZYZ`, `ZXZ`,
                                 `XZY`, `XYZ`, `YXZ`, `YZX`, `ZYX`, `ZXY`

        Raises:
            ValueError: If order is not a valid Euler convention
        """
        ...

class Shape:
    """
    GOAD shape object.
    """

    def __init__(
        self,
        vertices: list[float],
        face_indices: list[list[int]],
        id: int,
        refr_index_re: float,
        refr_index_im: float,
    ) -> None:
        """
        Create a new GOAD shape object.

        Args:
            vertices: List of vertex coordinates.
            face_indices: List of face index arrays. Integers must be positive and correspond to elements of the vertex list.
            id: Unique identifier for the shape.
            refr_index_re: Real part of the refractive index.
            refr_index_im: Imaginary part of the refractive index.

        """
        ...

class Geometry:
    """
    A GOAD geometry object. Represents a collection of one or more `Shape` objects.
    """
    def __init__(
        self,
        shapes: list[Shape],
    ) -> None:
        """
        Create a new GOAD geometry object.

        Args:
            shapes: List of GOAD shape objects.
        """
    ...

class Orientation:
    """
    A GOAD orientation object. Represents the orientation of a `Shape` object.
    """
    @staticmethod
    def uniform(
        num_orients: int, euler_convention: EulerConvention = ...
    ) -> Orientation:
        """
        A random uniform distribution of orientations with the specified number of orientations.

        Args:
            num_orients: Number of orientations.
            euler_convention: Optional Euler convention to use for the orientations.
        """

    @staticmethod
    def discrete(
        eulers: list[Euler], euler_convention: EulerConvention = ...
    ) -> Orientation:
        """
        A discrete distribution of orientations with the specified Euler angles.

        Args:
            eulers: List of Euler angles.
            euler_convention: Optional Euler convention to use for the orientations.
        """

class BinningScheme:
    """
    A GOAD binning scheme. Represents the query points for computing the far field.
    """
    @staticmethod
    def simple(
        num_theta: int,
        num_phi: int,
    ) -> BinningScheme:
        """
        Create a simple binning scheme with uniform theta and phi spacing

        Args:
            num_theta: Number of theta bins.
            num_phi: Number of phi bins.

        Raises:
            ValueError: If num_theta or num_phi is less than 1.
        """
    ...

    @staticmethod
    def interval(
        thetas: list[float],
        theta_spacings: list[float],
        phis: list[float],
        phi_spacings: list[float],
    ) -> BinningScheme:
        """
        Create an interval binning scheme with variable spacing.

        Args:
            thetas: Theta split points
            theta_spacings: Spacing for each theta interval
            phis: Phi split points
            phi_spacings: Spacing for each phi interval

        Examples:
            Theta in steps of 1 from 0 to 90 degrees, and then in steps of 2 from 90 to 180 degrees. Phi in steps of 2 from 0 to 360 degrees.
            >>> BinningScheme.interval(thetas=[0, 90, 180], theta_spacings=[1, 2], phis=[0, 360], phi_spacings=[2])

        Raises:
            Internal panic if the steps do not match the split points. ie. thetas=[0,180],theta_spacings=[7] would result in an error because 180 is not divisible by 7.
        """
        ...

    @staticmethod
    def custom(bins: list[list[list[float]]]) -> BinningScheme:
        """
        Create a custom binning scheme from explicit bin edges.

        Args:
            bins: List of bins, each specified as [[theta_min, theta_max], [phi_min, phi_max]]

        Examples:
            >>> BinningScheme.custom(bins=[[[0, 90], [0, 360]], [[90, 180], [0, 360]]])
            BinningScheme(...)

        """
        ...

    def thetas(self) -> list[float]:
        """
        Returns a list of theta values for each bin in the scheme.
        """
        ...

    def phis(self) -> list[float]:
        """
        Returns a list of phi values for each bin in the scheme.
        """
        ...

class Mapping:
    """
    A GOAD mapping. Represents the method for mapping from the near to the far field.
    """
    def __init__(
        self,
        method: str,
    ) -> None:
        """
        Create a new GOAD mapping object.

        Args:
            method: The mapping method. Valid options are 'go' (Geometric Optics) and 'ad' (Aperture Diffraction).
        """
    ...

class Settings:
    """
    A GOAD settings object. A customisable configuration for the GOAD simulation.
    """
    def __init__(
        self,
        geom_path: str,
        wavelength: float = ...,
        particle_refr_index_re: float = ...,
        particle_refr_index_im: float = ...,
        medium_refr_index_re: float = ...,
        medium_refr_index_im: float = ...,
        orientation: Orientation = ...,
        binning: BinningScheme = ...,
        beam_power_threshold: float = ...,
        beam_area_threshold_fac: float = ...,
        cutoff: float = ...,
        max_rec: int = ...,
        max_tir: int = ...,
        scale: float = ...,
        directory: str = ...,
        mapping: Mapping = ...,
        coherence: bool = ...,
        quiet: bool = ...,
    ) -> None:
        """
        Create a new GOAD settings object.

        Args:

        """
    ...

class MultiProblem:
    """
    A GOAD problem. Multi stands for multi-orientation problem, but you can also use it for particles in fixed orientations.
    """

    def __init__(self, settings: Settings, geom: Geometry = ...) -> None:
        """
        Create a new multi-orientation problem.

        Args:
            settings: A GOAD settings object.
            geom: An optional GOAD geometry object.

        """
        ...

    def solve(self) -> None:
        """Solve the scattering problem."""
        ...

    @property
    def results(self) -> Results:
        """Return the results of the scattering problem."""
        ...

class Results:
    """GOAD scattering results."""

    @property
    def bins(self) -> list[list[float]]:
        """Return the 2d bins [[theta, phi], ...]"""
        ...

    @property
    def bins_1d(self) -> list[float] | None:
        """Return the 1d bins."""
        ...

    @property
    def mueller(self) -> list[list[float]]:
        """Return the Mueller matrix [[s11, s12, ..., s44], ...]"""
        ...

    @property
    def mueller_beam(self) -> list[list[float]]:
        """Return the beam component Mueller matrix [[s11, s12, ..., s44], ...]"""
        ...

    @property
    def mueller_ext(self) -> list[list[float]]:
        """Return the external diffraction Mueller matrix [[s11, s12, ..., s44], ...]"""
        ...

    @property
    def mueller_1d(self) -> list[list[float]]:
        """Return the phi-integrated Mueller matrix [[s11, s12, ..., s44], ...]"""
        ...

    @property
    def mueller_1d_beam(self) -> list[list[float]]:
        """Return the phi-integrated beam Mueller matrix [[s11, s12, ..., s44], ...]"""
        ...

    @property
    def mueller_1d_ext(self) -> list[list[float]]:
        """Return the phi-integrated external diffraction Mueller matrix [[s11, s12, ..., s44], ...]"""
        ...

    @property
    def asymmetry(self) -> float | None:
        """Return the asymmetry parameter."""
        ...

    @property
    def scat_cross(self) -> float | None:
        """Return the scattering cross-section."""
        ...

    @property
    def ext_cross(self) -> float | None:
        """Return the extinction cross-section."""
        ...

    @property
    def albedo(self) -> float | None:
        """Return the single scattering albedo."""
        ...

    @property
    def powers(self) -> dict[str, float]:
        """Return list of power contributions.

        Keys: input, output, absorbed, trnc_ref, trnc_rec, trnc_clip,
        trnc_energy, clip_err, trnc_area, trnc_cop, ext_diff, missing.
        """
        ...

    @property
    def params(self) -> dict[str, float | dict[str, float]]:
        """Return all parameters as a dictionary.

        Deprecated: Use specific parameter properties instead.
        """
        ...
