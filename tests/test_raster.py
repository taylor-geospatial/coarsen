import numpy as np
import pytest
import shapely

import coarsen


def banks(n=20):
    upper, lower = [(0, 0)], [(0, 0)]
    for i in range(n):
        upper.extend([(i, i + 1), (i + 1, i + 1)])
        lower.extend([(i + 1, i), (i + 1, i + 1)])
    return [
        shapely.Polygon([*upper, (n, n + 5), (-5, n + 5), (-5, 0)]),
        shapely.Polygon([*lower, (n + 5, n), (n + 5, -5), (0, -5)]),
    ]


@pytest.mark.parametrize('reverse', [False, True])
def test_paired_raster_banks(reverse):
    g = banks()
    if reverse:
        g = [shapely.reverse(g[1]), g[0]]
    result = coarsen.coverage_simplify(g, 0.1, raster_pixel_size=1, threads=2)
    assert shapely.is_valid(result).all()
    assert result[0].intersection(result[1]).area == 0
    lane = shapely.union_all(result)
    for x in range(3, 17):
        section = shapely.LineString([(x, x - 3), (x, x + 3)])
        assert section.difference(lane).length == pytest.approx(1)
    assert sum(shapely.get_num_coordinates(result)) < 30


def test_protect_third_parcel():
    g = [*banks(), shapely.box(9, 9, 10, 10)]
    aa, bb, third = coarsen.coverage_simplify(g, 0.1, raster_pixel_size=1)
    assert third.equals(g[2])
    assert aa.intersection(third).area == bb.intersection(third).area == 0


@pytest.mark.parametrize('pixel', [0, -1, np.nan, np.inf])
def test_invalid_pixel_size(pixel):
    with pytest.raises(ValueError, match='raster_pixel_size'):
        coarsen.coverage_simplify(banks(), 1, raster_pixel_size=pixel)


def test_true_edges_and_lanes():
    geoms = [shapely.box(0, 0, 10, 10), shapely.box(10, 0, 20, 10), shapely.box(21, 0, 30, 10)]
    expected = coarsen.coverage_simplify(geoms, 2)
    result = coarsen.coverage_simplify(geoms, 2, raster_pixel_size=1)
    np.testing.assert_array_equal(shapely.to_wkb(result), shapely.to_wkb(expected))
