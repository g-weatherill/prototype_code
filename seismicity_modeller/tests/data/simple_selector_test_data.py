#!/usr/bin/env/python

'''Tests for Selector Class'''
import numpy as np

simple_test = {'eventID': np.array([001, 002, 003, 004, 005, 006, 007, 
                                    008, 009, 010], dtype=int),
               'year': np.array([1905, 1913, 1923, 1955, 1968, 1978,
                                 1988, 1989, 1995, 1997, 1999], dtype=int),
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
               'magnitude': np.arange([4., 4.5, 4.7,  5., 5.5, 6., 6.5, 6.8, 
                                       7., 7.5]),
               'xyz': None}