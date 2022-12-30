from polygon import Polygon, Point, is_line_segment_in_polygon, regular_polygon
import numpy as np
import sys, math, random

sys.path.insert(1, '../')

def test_regular_polygon():
    n = random.randint(3, 10)
    polygon = regular_polygon(n)
    
    assert is_line_segment_in_polygon([Point(0, 0), Point(0.1, 0.1)], polygon)
    assert not is_line_segment_in_polygon([Point(0, 0), Point(1, 1)], polygon)
    