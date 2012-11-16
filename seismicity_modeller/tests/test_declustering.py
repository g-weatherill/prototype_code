# -*- coding: utf-8 -*-

# Copyright (c) 2010-2012, GEM Foundation.
#
# OpenQuake is free software: you can redistribute it and/or modify it
# under the terms of the GNU Affero General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# OpenQuake is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with OpenQuake. If not, see <http://www.gnu.org/licenses/>.


import unittest
import numpy as np

from scientific.declustering import (TDW_GARDNERKNOPOFF,
    TDW_GRUENTHAL, TDW_UHRHAMMER, GardnerKnopoffType1, Afteran)

from data._declustering_test_data import (
    CATALOG_MATRIX_ALL_IN_A_CLUSTER, CATALOG_MATRIX_NO_CLUSTERS)


class DeclusteringTestCase(unittest.TestCase):

    def setUp(self):
        #self.catalog_matrix_all_cluster = np.array(
        #    CATALOG_MATRIX_ALL_IN_A_CLUSTER)
        self.catalog_matrix_all_cluster = {
            'year': CATALOG_MATRIX_ALL_IN_A_CLUSTER[:, 0],
            'month': CATALOG_MATRIX_ALL_IN_A_CLUSTER[:, 1],
            'day': CATALOG_MATRIX_ALL_IN_A_CLUSTER[:, 2],
            'longitude': CATALOG_MATRIX_ALL_IN_A_CLUSTER[:, 3],
            'latitude': CATALOG_MATRIX_ALL_IN_A_CLUSTER[:, 4],
            'magnitude': CATALOG_MATRIX_ALL_IN_A_CLUSTER[:, 5]}
            
        #self.catalog_matrix_no_clusters = np.array(CATALOG_MATRIX_NO_CLUSTERS)
        self.catalogue_matrix_no_clusters = {
            'year': CATALOG_MATRIX_NO_CLUSTERS[:, 0],
            'month': CATALOG_MATRIX_NO_CLUSTERS[:, 1],
            'day': CATALOG_MATRIX_NO_CLUSTERS[:, 2],
            'longitude': CATALOG_MATRIX_NO_CLUSTERS[:, 3],
            'latitude': CATALOG_MATRIX_NO_CLUSTERS[:, 4],
            'magnitude': CATALOG_MATRIX_NO_CLUSTERS[:, 5],
            
        self.time_dist_windows_options = [TDW_GARDNERKNOPOFF,
                TDW_GRUENTHAL, TDW_UHRHAMMER]

        self.foreshock_time_windows = [0.0, 0.1, 0.2]

        self.time_window_in_days = [10, 30, 60]

    def evaluate_results_gardner(self, catalog, exp_vcl,
            exp_vmain_shock, exp_flag_vector):
        config = {'window_opt': None,
                  'algorithm': 'GardnerKnopoffType1',
                  'fs_time_prop': None,
                  'purge': None}
        for tdw in self.time_dist_windows_options:
            config['window_opt'] = tdw
            for ftw in self.foreshock_time_windows:
                config['fs_time_prop'] = ftw
                vcl, vmain_shock, flag_vector = GardnerKnopoffType1.decluster(
                    catalogue, config)

                self.assertTrue(np.array_equal(exp_vcl, vcl))

                #self.assertTrue(np.array_equal(
                #        exp_vmain_shock, vmain_shock))

                self.assertTrue(np.array_equal(
                        exp_flag_vector, flag_vector))

    def evaluate_results_afteran(self, catalog_matrix, exp_vcl,
        exp_vmain_shock, exp_flag_vector):
        config = {'window opt': None,
                  'algorithm': 'Afteran',
                  'time_window': None,
                  'purge': None}
        for tdw in self.time_dist_windows_options:
            config['window_opt'] = tdw
            for tw in self.time_window_in_days:
                config['time_window'] = tw
                vcl, vmain_shock, flag_vector = Afteran.decluster(
                    catalog_matrix, config)

                self.assertTrue(np.array_equal(exp_vcl, vcl))

                #self.assertTrue(np.array_equal(
                #        exp_vmain_shock, vmain_shock))

                self.assertTrue(np.array_equal(
                        exp_flag_vector, flag_vector))

    def test_gardner_knopoff_all_events_within_a_cluster(self):

        expected_vcl = np.array([1, 1, 1, 1, 1, 1, 1, 1,
                1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1])

        #expected_vmain_shock = np.array([[1.9820e+03, 5.0000e+00, 3.000e+00,
        #        2.0596e+01, 3.8545e+01, 7.2000e+00, 1.0000e-01]])

        expected_flag_vector = np.array([0, 1, 1, 1, 1, 1, 1, 1,
                1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1])

        self.evaluate_results_gardner(self.catalog_matrix_all_cluster,
                expected_vcl, expected_vmain_shock, expected_flag_vector)

    def test_gardner_knopoff_no_events_within_a_cluster(self):

        expected_vcl = expected_flag_vector = np.array([0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])

        #expected_vmain_shock = self.catalog_matrix_no_clusters

        self.evaluate_results_gardner(self.catalog_matrix_no_clusters,
                expected_vcl, expected_vmain_shock, expected_flag_vector)

    def test_afteran_all_events_within_a_cluster(self):

        expected_vcl = np.array([1, 1, 1, 1, 1, 1, 1, 1,
                1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1])

        #expected_vmain_shock = np.array([[1.9820e+03, 5.0000e+00, 3.000e+00,
        #        2.0596e+01, 3.8545e+01, 7.2000e+00, 1.0000e-01]])

        expected_flag_vector = np.array([0, 1, 1, 1, 1, 1, 1, 1,
                1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1])

        self.evaluate_results_afteran(self.catalog_matrix_all_cluster,
            expected_vcl, expected_vmain_shock, expected_flag_vector)

    def test_afteran_no_events_within_a_cluster(self):

        expected_vcl = expected_flag_vector = np.array([0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])

        expected_vmain_shock = self.catalog_matrix_no_clusters

        self.evaluate_results_afteran(self.catalog_matrix_no_clusters,
                expected_vcl, expected_vmain_shock, expected_flag_vector)