"""Oblique torus handling for MCNP writer — surface-side TR card approach."""

import logging
import math

import FreeCAD

logger = logging.getLogger("general_logger")


def build_rotation_matrix(direction: FreeCAD.Vector) -> tuple:
    """Build orthonormal rotation matrix B such that B @ (0,0,1) = direction.

    Parameters
    ----------
    direction : FreeCAD.Vector
        Unit vector of the torus axis in the global frame.
        Must be normalized (length ≈ 1).

    Returns
    -------
    tuple of 9 floats
        (b1, b2, b3, b4, b5, b6, b7, b8, b9) — row-major B matrix.
        Row 1 = local x' in global, Row 2 = local y' in global,
        Row 3 = local z' in global = direction.

    Raises
    ------
    ValueError
        If direction has near-zero length.
    """
    d = FreeCAD.Vector(direction)
    if d.Length < 1e-12:
        raise ValueError(f"Torus axis direction has near-zero length: {direction}")
    d.normalize()

    # z_local = d (the torus axis becomes local z)
    z_loc = d

    # Choose reference vector for cross product to avoid near-parallel case
    if abs(d.z) < 0.9:
        ref = FreeCAD.Vector(0, 0, 1)
    else:
        ref = FreeCAD.Vector(0, 1, 0)

    # x_local = normalize(ref × z_local)
    x_loc = ref.cross(z_loc)
    x_loc.normalize()

    # y_local = z_local × x_local  (right-handed system)
    y_loc = z_loc.cross(x_loc)
    y_loc.normalize()

    # Verify proper rotation (det = +1)
    det = (
        x_loc.x * (y_loc.y * z_loc.z - y_loc.z * z_loc.y)
        - x_loc.y * (y_loc.x * z_loc.z - y_loc.z * z_loc.x)
        + x_loc.z * (y_loc.x * z_loc.y - y_loc.y * z_loc.x)
    )
    if det < 0:
        # Flip x to get proper rotation
        x_loc = x_loc * (-1)
        y_loc = z_loc.cross(x_loc)
        y_loc.normalize()

    return (
        x_loc.x, x_loc.y, x_loc.z,
        y_loc.x, y_loc.y, y_loc.z,
        z_loc.x, z_loc.y, z_loc.z,
    )


class ObliqueTorusRegistry:
    """Registry for oblique torus TR cards during MCNP writing.

    Each registered oblique torus produces one TR card. The registry
    assigns monotonically increasing TR numbers.

    Attributes
    ----------
    entries : list of dict
        Each entry: {
            "surf_id": int,       # original surface ID
            "tr_id": int,         # assigned TR card number
            "origin": (float, float, float),  # torus center in cm (global)
            "cosines": tuple of 9 floats,     # B matrix row-major
        }
    """

    def __init__(self, tr_start: int = 9900):
        """
        Parameters
        ----------
        tr_start : int
            First TR card number to assign. Subsequent ones increment by 1.
        """
        self._next_tr = tr_start
        self.entries = []

    @property
    def next_tr_id(self) -> int:
        return self._next_tr

    def register(self, surf_id: int, direction: FreeCAD.Vector,
                 center_cm: FreeCAD.Vector) -> int:
        """Register an oblique torus and return the assigned TR card number.

        Parameters
        ----------
        surf_id : int
            The GEOUNED surface index for this torus.
        direction : FreeCAD.Vector
            Unit axis of the torus (global frame, already normalized).
        center_cm : FreeCAD.Vector
            Torus center in cm (already converted from mm).

        Returns
        -------
        int
            The TR card number assigned to this surface.
        """
        cosines = build_rotation_matrix(direction)
        tr_id = self._next_tr
        self._next_tr += 1

        entry = {
            "surf_id": surf_id,
            "tr_id": tr_id,
            "origin": (center_cm.x, center_cm.y, center_cm.z),
            "cosines": cosines,
        }
        self.entries.append(entry)

        logger.info(
            f"Oblique torus surface {surf_id}: assigned TR{tr_id}, "
            f"axis=({direction.x:.6f}, {direction.y:.6f}, {direction.z:.6f}), "
            f"center=({center_cm.x:.4f}, {center_cm.y:.4f}, {center_cm.z:.4f}) cm"
        )
        return tr_id

    def format_tr_cards(self, fmt_origin: str = "14.7e",
                        fmt_cosine: str = "20.15e") -> str:
        """Format all registered TR cards as a single MCNP-ready string.

        Parameters
        ----------
        fmt_origin : str
            Python format spec for origin coordinates.
        fmt_cosine : str
            Python format spec for direction cosines (need ≥12 sig digits).

        Returns
        -------
        str
            Multi-line string ready to write into MCNP data block.
            Empty string if no oblique tori registered.
        """
        if not self.entries:
            return ""

        lines = []
        lines.append("C ")
        lines.append("C " + "=" * 58)
        lines.append("C     COORDINATE TRANSFORMATIONS (oblique tori)")
        lines.append("C " + "=" * 58)
        lines.append("C ")

        for entry in self.entries:
            ox, oy, oz = entry["origin"]
            b = entry["cosines"]
            tr_id = entry["tr_id"]

            card = "TR{:<5d}".format(tr_id)
            # Line 1: origin
            line1 = "{} {:{fo}} {:{fo}} {:{fo}}".format(
                card, ox, oy, oz, fo=fmt_origin
            )
            # Line 2-4: cosines (3 per line)
            line2 = "{:8s}{:{fc}} {:{fc}} {:{fc}}".format(
                "", b[0], b[1], b[2], fc=fmt_cosine
            )
            line3 = "{:8s}{:{fc}} {:{fc}} {:{fc}}".format(
                "", b[3], b[4], b[5], fc=fmt_cosine
            )
            line4 = "{:8s}{:{fc}} {:{fc}} {:{fc}} 1".format(
                "", b[6], b[7], b[8], fc=fmt_cosine
            )

            lines.append(f"C  Surface {entry['surf_id']} oblique torus")
            lines.append(line1)
            lines.append(line2)
            lines.append(line3)
            lines.append(line4)

        return "\n".join(lines) + "\n"
