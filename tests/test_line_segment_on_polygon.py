from polygon import Polygon, Point, is_line_segment_in_polygon, regular_polygon
import numpy as np
import sys, math, random

sys.path.insert(1, '../')

def test_regular_polygon():
    n = random.randint(3, 10)
    polygon = regular_polygon(n)
    
    assert is_line_segment_in_polygon([Point(0, 0), Point(0.1, 0.1)], polygon)
    assert not is_line_segment_in_polygon([Point(0, 0), Point(1, 1)], polygon)
    
def test_non_convex():
    polygon = Polygon(5, [Point(0, 0), Point(2, 0), Point(2, 2), Point(1, 1), Point(0, 2)])
    assert not is_line_segment_in_polygon([Point(0, 2), Point(2, 2)], polygon)
    assert not is_line_segment_in_polygon([Point(0, 2), Point(1, 2)], polygon)
    assert is_line_segment_in_polygon([Point(0, 2), Point(1, 1)], polygon)
    
    # Need to fix this
    assert is_line_segment_in_polygon([Point(0, 2), Point(1.1, 0.9)], polygon)
    assert is_line_segment_in_polygon([Point(1.1, 0.9), Point(0, 2)], polygon)
    
def test_rectangle():
    polygon = Polygon(6, [Point(0, 0), Point(1, 0), Point(2, 0), Point(2, 1), Point(1, 1), Point(0, 1)])
    assert is_line_segment_in_polygon([Point(0, 0), Point(2, 1)], polygon)
    assert is_line_segment_in_polygon([Point(0, 0), Point(2, 0)], polygon)