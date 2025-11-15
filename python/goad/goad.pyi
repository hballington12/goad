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
