#!/usr/bin/env/python

'''Tests for Selector Class'''
import unittest
import numpy as np
from nhlib.geo.point import Point
from nhlib.geo.polygon import Polygon
from scientific.selector import Selector


simple_test = {'eventID': np.array([101, 102, 103, 104, 105, 106, 107, 
                                    108, 109, 010], dtype=int),
               'year': np.array([1905, 1913, 1923, 1955, 1968, 1978,
                                 1988, 1995, 1997, 1999], dtype=int),
               'month': np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10], dtype=int),
               'day': np.array([2, 4, 6, 8, 12, 14, 16, 18, 20, 22], 
                               dtype=int),
               'hour': np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10], dtype=int),
               'minute': np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10], dtype=int),
               'second': np.array([1., 1., 1., 1., 1., 1., 1., 1., 1., 1.], 
                                  dtype=int),
               'longitude': np.arange(30., 40., 1.),
               'latitude': np.arange(40., 50., 1.),
               'depth' : np.array([5., 10., 5., 10., 5., 10., 5., 10., 5., 
                                   10.]),
               'magnitude': np.array([4., 4.5, 4.7,  5., 5.5, 6., 6.5, 6.8, 
                                       7., 7.5]),
               'xyz': None}



simple_polygon = Polygon([Point(33.5, 43.5), 
                          Point(33.5, 46.5),
                          Point(36.5, 46.5),
                          Point(36.5, 43.5),
                          Point(33.5, 43.5)])

xyz_check = array([[ 4223.29319404,  2438.31946245,  4091.98592326],
       [ 4115.00960906,  2472.54722101,  4173.19148341],
       [ 4011.9967814 ,  2506.97382866,  4259.68544008],
       [ 3901.61364741,  2533.73752835,  4338.19156836],
       [ 3796.42598971,  2560.72166362,  4422.19518634],
       [ 3684.46908752,  2579.89302939,  4497.90623513],
       [ 3577.63105819,  2599.30111329,  4579.31716896],
       [ 3464.63383588,  2610.78885915,  4652.140896  ],
       [ 3356.67793377,  2622.52422246,  4730.85995899],
       [ 3243.17890772,  2626.27449626,  4800.7076398 ]])


class TestSelector(unittest.TestCase):
    def setUp(self):
        '''Set up tests'''
        self.input_catalogue = simple_test
        self.polygon = simple_polygon
        self.selector = None 
        Selector(self.input_catalogue)

    def test_basic_read(self):
        '''Checks basic input process'''
        self.selector = Selector(self.input_catalogue)
        self.assertEqual(self.selector.number_events, 10)
        self.assertTrue(np.allclose(
            self.selector.catalogue['longitude'],
            np.array([30., 31., 31., 33., 34., 35., 36., 37., 38., 39.])))
        self.assertTrue(np.allclose(self.selector.catalogue['xyz'], 
                                    xyz_check))

    def test_point_in_polygon(



