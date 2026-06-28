import math
import sys

import pytest

# Import the data_classes module directly to avoid the FreeCAD-requiring
# package __init__.py chain.  Options and NumericFormat are pure Python.
def _import_data_classes():
    import importlib.util

    spec = importlib.util.spec_from_file_location(
        "geouned.GEOUNED.utils.data_classes",
        "src/geouned/GEOUNED/utils/data_classes.py",
    )
    mod = importlib.util.module_from_spec(spec)
    sys.modules["geouned.GEOUNED.utils.data_classes"] = mod
    spec.loader.exec_module(mod)
    return mod


_data_classes = _import_data_classes()
Options = _data_classes.Options
NumericFormat = _data_classes.NumericFormat

try:
    import FreeCAD
    import Part

    HAS_FREECAD = True
except ImportError:
    HAS_FREECAD = False


# ---------------------------------------------------------------------------
# Options & NumericFormat (no FreeCAD needed)
# ---------------------------------------------------------------------------


class TestObliqueTorusOptions:

    def test_default_strategy_is_tr(self):
        options = Options()
        assert options.obliqueTorusStrategy == "tr"

    def test_set_strategy_gq_fallback(self):
        options = Options(obliqueTorusStrategy="gq_fallback")
        assert options.obliqueTorusStrategy == "gq_fallback"

    def test_set_strategy_tr_explicit(self):
        options = Options(obliqueTorusStrategy="tr")
        assert options.obliqueTorusStrategy == "tr"

    def test_invalid_strategy_raises_value_error(self):
        with pytest.raises(ValueError, match="obliqueTorusStrategy"):
            Options(obliqueTorusStrategy="invalid")

    def test_non_string_strategy_raises_type_error(self):
        with pytest.raises(TypeError, match="obliqueTorusStrategy"):
            Options(obliqueTorusStrategy=123)

    def test_strategy_setter_validation(self):
        options = Options()
        options.obliqueTorusStrategy = "gq_fallback"
        assert options.obliqueTorusStrategy == "gq_fallback"
        options.obliqueTorusStrategy = "tr"
        assert options.obliqueTorusStrategy == "tr"
        with pytest.raises(ValueError):
            options.obliqueTorusStrategy = "bad_value"


class TestNumericFormatTRFields:

    def test_default_tr_o(self):
        nf = NumericFormat()
        assert nf.TR_o == "14.7e"

    def test_default_tr_cos(self):
        nf = NumericFormat()
        assert nf.TR_cos == "20.15e"

    def test_custom_tr_o(self):
        nf = NumericFormat(TR_o="12.5e")
        assert nf.TR_o == "12.5e"

    def test_custom_tr_cos(self):
        nf = NumericFormat(TR_cos="18.12e")
        assert nf.TR_cos == "18.12e"

    def test_tr_o_type_check(self):
        nf = NumericFormat()
        with pytest.raises(TypeError, match="TR_o"):
            nf.TR_o = 123

    def test_tr_cos_type_check(self):
        nf = NumericFormat()
        with pytest.raises(TypeError, match="TR_cos"):
            nf.TR_cos = 456


# ---------------------------------------------------------------------------
# build_rotation_matrix (FreeCAD required)
# ---------------------------------------------------------------------------


@pytest.mark.skipif(not HAS_FREECAD, reason="FreeCAD not available in this environment")
class TestBuildRotationMatrix:

    def test_axis_aligned_x(self):
        from geouned.GEOUNED.write.oblique_torus_tr import build_rotation_matrix

        direction = FreeCAD.Vector(1, 0, 0)
        b = build_rotation_matrix(direction)
        assert len(b) == 9
        self._assert_orthonormal(b)
        self._assert_proper_rotation(b)
        self._assert_row3_matches_direction(b, direction)

    def test_axis_aligned_y(self):
        from geouned.GEOUNED.write.oblique_torus_tr import build_rotation_matrix

        direction = FreeCAD.Vector(0, 1, 0)
        b = build_rotation_matrix(direction)
        self._assert_orthonormal(b)
        self._assert_proper_rotation(b)
        self._assert_row3_matches_direction(b, direction)

    def test_axis_aligned_z(self):
        from geouned.GEOUNED.write.oblique_torus_tr import build_rotation_matrix

        direction = FreeCAD.Vector(0, 0, 1)
        b = build_rotation_matrix(direction)
        self._assert_orthonormal(b)
        self._assert_proper_rotation(b)
        self._assert_row3_matches_direction(b, direction)

    def test_oblique_diagonal(self):
        from geouned.GEOUNED.write.oblique_torus_tr import build_rotation_matrix

        direction = FreeCAD.Vector(1, 1, 1)
        b = build_rotation_matrix(direction)
        self._assert_orthonormal(b)
        self._assert_proper_rotation(b)
        self._assert_row3_matches_direction(b, direction)

    def test_oblique_arbitrary(self):
        from geouned.GEOUNED.write.oblique_torus_tr import build_rotation_matrix

        direction = FreeCAD.Vector(0.577350269, 0.577350269, 0.577350269)
        b = build_rotation_matrix(direction)
        self._assert_orthonormal(b)
        self._assert_proper_rotation(b)
        self._assert_row3_matches_direction(b, direction)

    def test_near_parallel_z_positive(self):
        from geouned.GEOUNED.write.oblique_torus_tr import build_rotation_matrix

        direction = FreeCAD.Vector(0.01, 0.02, 0.9998)
        b = build_rotation_matrix(direction)
        self._assert_orthonormal(b)
        self._assert_proper_rotation(b)
        self._assert_row3_matches_direction(b, direction)

    def test_near_parallel_z_negative(self):
        from geouned.GEOUNED.write.oblique_torus_tr import build_rotation_matrix

        direction = FreeCAD.Vector(0.01, 0.02, -0.9998)
        b = build_rotation_matrix(direction)
        self._assert_orthonormal(b)
        self._assert_proper_rotation(b)
        self._assert_row3_matches_direction(b, direction)

    def test_zero_length_raises_value_error(self):
        from geouned.GEOUNED.write.oblique_torus_tr import build_rotation_matrix

        direction = FreeCAD.Vector(0, 0, 0)
        with pytest.raises(ValueError, match="near-zero length"):
            build_rotation_matrix(direction)

    def test_very_small_direction_raises_value_error(self):
        from geouned.GEOUNED.write.oblique_torus_tr import build_rotation_matrix

        direction = FreeCAD.Vector(1e-14, 0, 0)
        with pytest.raises(ValueError, match="near-zero length"):
            build_rotation_matrix(direction)

    def test_negated_direction_is_still_proper(self):
        from geouned.GEOUNED.write.oblique_torus_tr import build_rotation_matrix

        direction = FreeCAD.Vector(-1, -2, -3)
        b = build_rotation_matrix(direction)
        self._assert_orthonormal(b)
        self._assert_proper_rotation(b)
        self._assert_row3_matches_direction(b, direction)

    # -- helpers ----------------------------------------------------------

    @staticmethod
    def _assert_orthonormal(b):
        def _row(i):
            return FreeCAD.Vector(b[i * 3], b[i * 3 + 1], b[i * 3 + 2])

        r0 = _row(0)
        r1 = _row(1)
        r2 = _row(2)
        assert math.isclose(r0.Length, 1.0, abs_tol=1e-12), f"row 0 not unit: {r0.Length}"
        assert math.isclose(r1.Length, 1.0, abs_tol=1e-12), f"row 1 not unit: {r1.Length}"
        assert math.isclose(r2.Length, 1.0, abs_tol=1e-12), f"row 2 not unit: {r2.Length}"
        assert math.isclose(r0.dot(r1), 0.0, abs_tol=1e-12), f"row 0·1 not orthogonal: {r0.dot(r1)}"
        assert math.isclose(r1.dot(r2), 0.0, abs_tol=1e-12), f"row 1·2 not orthogonal: {r1.dot(r2)}"
        assert math.isclose(r2.dot(r0), 0.0, abs_tol=1e-12), f"row 2·0 not orthogonal: {r2.dot(r0)}"

    @staticmethod
    def _assert_proper_rotation(b):
        det = (
            b[0] * (b[4] * b[8] - b[5] * b[7])
            - b[1] * (b[3] * b[8] - b[5] * b[6])
            + b[2] * (b[3] * b[7] - b[4] * b[6])
        )
        assert math.isclose(det, 1.0, abs_tol=1e-12), f"determinant {det} ≠ 1"

    @staticmethod
    def _assert_row3_matches_direction(b, direction):
        expected = FreeCAD.Vector(direction).normalize()
        row3 = FreeCAD.Vector(b[6], b[7], b[8])
        assert row3.isEqual(expected, 1e-12), f"row3 {row3} ≠ direction {expected}"


# ---------------------------------------------------------------------------
# ObliqueTorusRegistry (FreeCAD required)
# ---------------------------------------------------------------------------


@pytest.mark.skipif(not HAS_FREECAD, reason="FreeCAD not available in this environment")
class TestObliqueTorusRegistry:

    def test_default_constructor(self):
        from geouned.GEOUNED.write.oblique_torus_tr import ObliqueTorusRegistry

        reg = ObliqueTorusRegistry()
        assert reg.next_tr_id == 9900
        assert reg.entries == []

    def test_custom_tr_start(self):
        from geouned.GEOUNED.write.oblique_torus_tr import ObliqueTorusRegistry

        reg = ObliqueTorusRegistry(tr_start=5000)
        assert reg.next_tr_id == 5000

    def test_register_returns_monotonic_ids(self):
        from geouned.GEOUNED.write.oblique_torus_tr import ObliqueTorusRegistry

        reg = ObliqueTorusRegistry(tr_start=100)
        id1 = reg.register(10, FreeCAD.Vector(1, 0, 0), FreeCAD.Vector(0, 0, 0))
        id2 = reg.register(20, FreeCAD.Vector(0, 1, 0), FreeCAD.Vector(1, 2, 3))
        id3 = reg.register(30, FreeCAD.Vector(0, 0, 1), FreeCAD.Vector(4, 5, 6))
        assert id1 == 100
        assert id2 == 101
        assert id3 == 102
        assert reg.next_tr_id == 103

    def test_register_entry_structure(self):
        from geouned.GEOUNED.write.oblique_torus_tr import ObliqueTorusRegistry

        reg = ObliqueTorusRegistry(tr_start=1)
        direction = FreeCAD.Vector(3, 0, 0)
        center = FreeCAD.Vector(10, 20, 30)
        tid = reg.register(99, direction, center)
        assert tid == 1
        entry = reg.entries[0]
        assert entry["surf_id"] == 99
        assert entry["tr_id"] == 1
        assert entry["origin"] == (10.0, 20.0, 30.0)
        assert len(entry["cosines"]) == 9

    def test_format_tr_cards_empty(self):
        from geouned.GEOUNED.write.oblique_torus_tr import ObliqueTorusRegistry

        reg = ObliqueTorusRegistry()
        output = reg.format_tr_cards()
        assert output == ""

    def test_format_tr_cards_single_entry(self):
        from geouned.GEOUNED.write.oblique_torus_tr import ObliqueTorusRegistry

        reg = ObliqueTorusRegistry(tr_start=1)
        reg.register(100, FreeCAD.Vector(0, 0, 1), FreeCAD.Vector(0, 0, 0))
        output = reg.format_tr_cards()
        assert "TR1" in output
        assert "oblique tori" in output
        assert "Surface 100" in output

    def test_format_tr_cards_multiple_entries(self):
        from geouned.GEOUNED.write.oblique_torus_tr import ObliqueTorusRegistry

        reg = ObliqueTorusRegistry(tr_start=1)
        reg.register(100, FreeCAD.Vector(0, 0, 1), FreeCAD.Vector(0, 0, 0))
        reg.register(200, FreeCAD.Vector(1, 0, 0), FreeCAD.Vector(5, 5, 5))
        reg.register(300, FreeCAD.Vector(0, 1, 0), FreeCAD.Vector(-1, -1, -1))
        output = reg.format_tr_cards()
        assert "TR1" in output
        assert "TR2" in output
        assert "TR3" in output
        assert "Surface 100" in output
        assert "Surface 200" in output
        assert "Surface 300" in output

    def test_format_tr_cards_numbering_continues_after_registration(self):
        from geouned.GEOUNED.write.oblique_torus_tr import ObliqueTorusRegistry

        reg = ObliqueTorusRegistry(tr_start=1)
        reg.register(10, FreeCAD.Vector(0, 0, 1), FreeCAD.Vector(0, 0, 0))
        output = reg.format_tr_cards()
        assert "TR1" in output
        assert "TR2" not in output

    def test_register_cosines_match_build_rotation_matrix(self):
        from geouned.GEOUNED.write.oblique_torus_tr import ObliqueTorusRegistry, build_rotation_matrix

        direction = FreeCAD.Vector(1, 2, 3)
        reg = ObliqueTorusRegistry(tr_start=1)
        reg.register(50, direction, FreeCAD.Vector(0, 0, 0))
        expected = build_rotation_matrix(direction)
        actual = reg.entries[0]["cosines"]
        for i in range(9):
            assert math.isclose(actual[i], expected[i], abs_tol=1e-12), f"index {i}: {actual[i]} ≠ {expected[i]}"


# ---------------------------------------------------------------------------
# GQ quadric fallback — oblique torus → circumscribed cylinder (FreeCAD required)
# ---------------------------------------------------------------------------


def _make_torus_surface(axis, center, major_r, minor_r):
    axis = FreeCAD.Vector(axis)
    axis.normalize()
    p = Part.makeTorus(major_r, minor_r, FreeCAD.Vector(0, 0, 0), axis, 360)
    p.translate(center)
    face = p.Faces[0]
    from geouned.GEOUNED.utils.geometry_gu import getSurface

    surf = getSurface(face)
    surf.MajorRadius = major_r
    surf.MinorRadius = minor_r
    return surf


@pytest.mark.skipif(not HAS_FREECAD, reason="FreeCAD not available in this environment")
class TestObliqueTorusGQFallback:

    def test_oblique_torus_omc_surface_is_quadric_not_none(self):
        from geouned.GEOUNED.utils.data_classes import NumericFormat, Tolerances
        from geouned.GEOUNED.write.functions import open_mc_surface

        surf = _make_torus_surface(FreeCAD.Vector(1, 1, 1), FreeCAD.Vector(0, 0, 0), 4, 1)
        omc_type, coeffs = open_mc_surface("Torus", surf, Tolerances(), NumericFormat())
        assert omc_type == "quadric"
        assert coeffs

    def test_oblique_torus_omc_coeffs_ten_tokens(self):
        from geouned.GEOUNED.utils.data_classes import NumericFormat, Tolerances
        from geouned.GEOUNED.write.functions import open_mc_surface

        surf = _make_torus_surface(FreeCAD.Vector(1, 1, 1), FreeCAD.Vector(0, 0, 0), 4, 1)
        _, coeffs = open_mc_surface("Torus", surf, Tolerances(), NumericFormat())
        assert len(coeffs.split()) == 10

    def test_oblique_torus_mcnp_gq_fallback(self):
        from geouned.GEOUNED.utils.data_classes import NumericFormat, Options, Tolerances
        from geouned.GEOUNED.write.functions import mcnp_surface

        surf = _make_torus_surface(FreeCAD.Vector(1, 1, 1), FreeCAD.Vector(0, 0, 0), 4, 1)
        options = Options(obliqueTorusStrategy="gq_fallback")
        result = mcnp_surface(10, "Torus", surf, options, Tolerances(), NumericFormat())
        assert "GQ" in result

    def test_oblique_torus_mcnp_tr_strategy(self):
        from geouned.GEOUNED.utils.data_classes import NumericFormat, Options, Tolerances
        from geouned.GEOUNED.write.functions import mcnp_surface
        from geouned.GEOUNED.write.oblique_torus_tr import ObliqueTorusRegistry

        surf = _make_torus_surface(FreeCAD.Vector(1, 1, 1), FreeCAD.Vector(0, 0, 0), 4, 1)
        options = Options(obliqueTorusStrategy="tr")
        reg = ObliqueTorusRegistry(tr_start=9900)
        result = mcnp_surface(10, "Torus", surf, options, Tolerances(), NumericFormat(), reg)
        assert "TZ" in result
        assert "9900" in result
        assert len(reg.entries) == 1

    def test_oblique_torus_mcnp_tr_strategy_no_registry_falls_back_to_gq(self):
        from geouned.GEOUNED.utils.data_classes import NumericFormat, Options, Tolerances
        from geouned.GEOUNED.write.functions import mcnp_surface

        surf = _make_torus_surface(FreeCAD.Vector(1, 1, 1), FreeCAD.Vector(0, 0, 0), 4, 1)
        options = Options(obliqueTorusStrategy="tr")
        result = mcnp_surface(10, "Torus", surf, options, Tolerances(), NumericFormat())
        assert "GQ" in result

    def test_axis_aligned_torus_still_works_omc(self):
        from geouned.GEOUNED.utils.data_classes import NumericFormat, Tolerances
        from geouned.GEOUNED.write.functions import open_mc_surface

        for axis in (FreeCAD.Vector(1, 0, 0), FreeCAD.Vector(0, 1, 0), FreeCAD.Vector(0, 0, 1)):
            surf = _make_torus_surface(axis, FreeCAD.Vector(0, 0, 0), 4, 1)
            omc_type, coeffs = open_mc_surface("Torus", surf, Tolerances(), NumericFormat())
            assert omc_type in ("x-torus", "y-torus", "z-torus")
            assert coeffs
            assert "quadric" not in omc_type

    def test_axis_aligned_torus_still_works_mcnp(self):
        from geouned.GEOUNED.utils.data_classes import NumericFormat, Options, Tolerances
        from geouned.GEOUNED.write.functions import mcnp_surface

        for axis, card in (
            (FreeCAD.Vector(1, 0, 0), "TX"),
            (FreeCAD.Vector(0, 1, 0), "TY"),
            (FreeCAD.Vector(0, 0, 1), "TZ"),
        ):
            surf = _make_torus_surface(axis, FreeCAD.Vector(0, 0, 0), 4, 1)
            result = mcnp_surface(10, "Torus", surf, Options(), Tolerances(), NumericFormat())
            assert card in result
            assert "GQ" not in result
