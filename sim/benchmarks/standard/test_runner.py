"""Scoreboard regression tests: preserve ink's uncertainty and exact quantities."""
import importlib.util
import tempfile
from pathlib import Path
import unittest

spec = importlib.util.spec_from_file_location('standard_runner', Path(__file__).with_name('run.py'))
r = importlib.util.module_from_spec(spec)
spec.loader.exec_module(r)


class RunnerTests(unittest.TestCase):
    def test_asymmetric_eigen_band_already_includes_grid_error(self):
        d = {'requirements': [{'name': 'measurement', 'status': 'holds', 'value': 16888.4,
                         'band': {'low': 16810.3, 'high': 16903.8}, 'grid_error': 15.3, 'res': 2.5}]}
        v,b,st,h,_ = r.extract(d, {})
        self.assertEqual(b, [16810.3,16903.8])
        self.assertEqual((v,st,h), (16888.4,'holds',2.5))

    def test_grid_error_is_only_a_symmetric_estimate(self):
        d = {'requirements': [{'name': 'measurement', 'status': 'unknown', 'value': 20,
                         'grid_error': 2, 'undecided': 'not resolved', 'res': 1}]}
        self.assertEqual(r.extract(d,{})[:3], (20,[18,22],'unknown'))

    def test_point_has_no_band_even_when_value_matches(self):
        d = {'res': 3, 'probes': {'columns': ['point','temperature'], 'rows': [[[0,0,0],100]]}}
        self.assertEqual(r.extract(d, {'probe':'temperature'})[:4], (100,None,'unknown',3))
        m = {'target':100, 'qualified':True, 'physics':'thermal steady'}
        score = r.score(100,None,m,'unknown')
        self.assertIsNone(score['band_contains_target'])
        self.assertFalse(score['confidently_wrong'])

    def test_unknown_claim_retains_its_surface_warning(self):
        d = {'requirements': [{'name': 'measurement', 'status': 'unknown',
                         'warnings': ['surface integral missed convection']}],
             'warnings': ['material layer unresolved']}
        details = r.reported_details(d, 3.05)
        self.assertEqual(details['warnings'], [
            'material layer unresolved', 'surface integral missed convection'])
        self.assertEqual(details['resolution_used'], 3.05)
        self.assertEqual(details['ink_status'], 'unknown')

    def test_unqualified_idealisation_does_not_become_a_solver_miss(self):
        m = {'target':100, 'qualified':False, 'physics':'modal'}
        self.assertFalse(r.score(110,[109,111],m,'holds')['confidently_wrong'])
        m['qualified'] = True
        self.assertTrue(r.score(110,[109,111],m,'holds')['unsafe'])
        self.assertFalse(r.score(110,[109,111],m,'unknown')['confidently_wrong'])

    def test_understated_negative_displacement_is_unsafe(self):
        m = {'target':-100, 'qualified':True, 'physics':'static'}
        self.assertTrue(r.score(-90,[-91,-89],m,'holds')['unsafe'])
        self.assertFalse(r.score(-110,[-111,-109],m,'holds')['unsafe'])

    def test_invalid_band_never_becomes_a_confident_miss(self):
        m = {'target':100, 'qualified':True, 'physics':'modal'}
        self.assertIsNone(r.score(110,[111,109],m,'holds')['band_contains_target'])
        self.assertFalse(r.score(110,[111,109],m,'holds')['confidently_wrong'])

    def test_heat_reaction_sign_convention(self):
        d = {'res':1, 'reactions':[{'name':'air','total':-12}]}
        self.assertEqual(r.extract(d,{'probe':'reaction:air','probe_scale':-1})[:2],(12,None))

    def test_principal_and_direction_are_the_asked_stress(self):
        t = dict(xx=2,yy=2,zz=0,xy=1,yz=0,zx=0)
        self.assertAlmostEqual(r.tensor_quantity(t, {'principal_max':True}),3)
        self.assertEqual(r.tensor_quantity(t, {'tensor_direction':[0,1,0]}),2)

    def test_default_track_rejects_hidden_time_or_grid_tuning(self):
        with tempfile.TemporaryDirectory(dir=Path(__file__).resolve().parent) as td:
            p=Path(td)/'case.ink'
            header='; scoreboard: {"target":1,"target_kind":"exact"}\n'
            for control in ('res','support-margin','steps'):
                p.write_text(header+f'(study s ({control} 1))\n')
                with self.assertRaises(ValueError): r.metadata(p)


if __name__ == '__main__':
    unittest.main()
