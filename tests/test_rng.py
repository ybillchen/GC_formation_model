# Licensed under BSD-3-Clause License - see LICENSE
"""RNG isolation and assignment reproducibility regression tests."""
from concurrent.futures import ProcessPoolExecutor
from copy import deepcopy
import importlib
import multiprocessing as mp
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
from GC_formation_model import astro_utils

assignment = importlib.import_module('GC_formation_model.assign')


def scatter_sequence(halo):
    rng = np.random.default_rng(np.random.SeedSequence([7, halo, 0x534D484D]))
    return [astro_utils.SMHM(1e11, 1, scatter=True, rng=rng) for _ in range(8)]


class RNGTests(unittest.TestCase):
    def test_deterministic_calls_consume_no_randomness(self):
        global_before = np.random.get_state()
        rng = np.random.default_rng(123)
        before = deepcopy(rng.bit_generator.state)
        with patch.object(np.random, 'normal', side_effect=AssertionError('global RNG used')):
            for function in [astro_utils.SMHMparameters, astro_utils.SMHMparameters2, astro_utils.SMHMparameters3]:
                self.assertEqual(function(1, rng=rng)[-1], 0)
                self.assertEqual(len(function(1)), 6)
            for k, mdef in [(False, 'm200'), (True, 'm200'), (True, 'mvir')]:
                self.assertEqual(astro_utils.SMHM(1e11, 1, k=k, mdef=mdef, rng=rng),
                                 astro_utils.SMHM(1e11, 1, k=k, mdef=mdef))
        self.assertEqual(before, rng.bit_generator.state)
        global_after = np.random.get_state()
        for a, b in zip(global_before, global_after):
            np.testing.assert_array_equal(a, b)

    def test_scatter_requires_explicit_rng(self):
        for function in [astro_utils.SMHMparameters, astro_utils.SMHMparameters2, astro_utils.SMHMparameters3]:
            with self.assertRaisesRegex(ValueError, 'explicit rng'):
                function(1, scatter=True)
        for k, mdef in [(False, 'm200'), (True, 'm200'), (True, 'mvir')]:
            with self.assertRaisesRegex(ValueError, 'explicit rng'):
                astro_utils.SMHM(1e11, 1, k=k, mdef=mdef, scatter=True)

    def test_scatter_is_one_explicit_draw(self):
        for k, mdef in [(False, 'm200'), (True, 'm200'), (True, 'mvir')]:
            rng = np.random.default_rng(123)
            expected_rng = np.random.default_rng(123)
            xi = expected_rng.normal(0, .218 - .023 * (1 / 2 - 1))
            median = astro_utils.SMHM(1e11, 1, k=k, mdef=mdef)
            with patch.object(np.random, 'normal', side_effect=AssertionError('global RNG used')):
                actual = astro_utils.SMHM(1e11, 1, k=k, mdef=mdef, scatter=True, rng=rng)
            self.assertAlmostEqual(np.log10(actual / median), xi)
            self.assertEqual(rng.bit_generator.state, expected_rng.bit_generator.state)

    def test_scatter_independent_of_order_and_process(self):
        serial = {halo: scatter_sequence(halo) for halo in [42, 84]}
        reversed_run = {halo: scatter_sequence(halo) for halo in [84, 42]}
        self.assertEqual(serial, reversed_run)
        with ProcessPoolExecutor(max_workers=2, mp_context=mp.get_context('spawn')) as pool:
            parallel = dict(zip([84, 42], pool.map(scatter_sequence, [84, 42])))
        self.assertEqual(serial, parallel)
        self.assertNotEqual(serial[42], serial[84])

    def test_collisionless_assignment_order_and_subset(self):
        for peaks in [False, True]:
            with self.subTest(peaks=peaks):
                expected = self.assignment_run([42, 84], peaks)
                self.assertEqual(expected, self.assignment_run([84, 42], peaks))
                self.assertEqual(expected[42], self.assignment_run([42], peaks)[42])

    def assignment_run(self, halos, peaks):
        coords = np.random.default_rng(333).normal(0, .2, (32, 3))
        def load_tree(base, halo):
            return dict(SnapNum=np.array([1]), SubfindID=np.array([halo]),
                        SubhaloPos=np.zeros((1, 3)), SubhaloMass=np.array([1.]), ScaleRad=np.array([1.]))
        def load_halo(base, root, halo, snapshot, kind, fields):
            if kind == 'dm':
                return dict(count=32, Coordinates=coords.copy(),
                            ParticleIDs=np.arange(root * 1000, root * 1000 + 32), Masses=np.ones(32))
            return dict(count=20, Coordinates=np.vstack([np.zeros((1, 3)), coords[:19]]),
                        Density=np.r_[1e12, np.ones(19)])
        class Cosmo:
            h = 1
            def cosmicTime(self, z, units=None):
                return 10 - z
        with tempfile.TemporaryDirectory() as directory:
            name = Path(directory) / 'catalog.txt'
            rows = [[halo, 10, 8, 10, 8, 5, 1, -1, 1, halo, 1] for halo in halos for _ in range(3)]
            np.savetxt(name, rows)
            np.savetxt(Path(directory) / 'catalog_offset_root.txt',
                       [[halo, 3*i, 3*i+3, i, i+1] for i, halo in enumerate(halos)], fmt='%d')
            np.savetxt(Path(directory) / 'catalog_offset.txt',
                       [[1, halo, 3*i, 3*i+3] for i, halo in enumerate(halos)], fmt='%d')
            params = dict(seed=7, verbose=False, t_lag=.01, redshift_snap=np.array([2., 1.]),
                          resultspath=directory+'/', allcat_name='catalog.txt', cosmo=Cosmo(),
                          base_tree='', base_halo='', h100=1, log_Mhmin=8,
                          collisionless_only=True, frac_rs=1, assign_at_peaks=peaks,
                          rmin_peaks_finding=.1, peak_radius=.15, evenly_distribute=False)
            with patch.object(assignment.loader, 'load_merger_tree', load_tree), \
                 patch.object(assignment.loader, 'load_halo', load_halo), \
                 patch.object(np.random, 'permutation', side_effect=AssertionError('global RNG used')):
                assignment.assign(params)
            result = np.loadtxt(Path(directory) / 'catalog_gcid.txt', ndmin=2, dtype=np.int64)
            return {halo: result[3*i:3*i+3].tolist() for i, halo in enumerate(halos)}


if __name__ == '__main__':
    unittest.main()
