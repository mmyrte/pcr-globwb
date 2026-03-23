import shutil
import tempfile
import unittest
from pathlib import Path
from pcrglobwb.configuration import Configuration

class TestConfiguration(unittest.TestCase):
    def setUp(self):
        # Create a temp output directory
        self.tmp_output = Path(tempfile.mkdtemp())
        # Copy the comprehensive ini file to a temp location and patch outputDir
        ini_sample = Path(__file__).parent / 'data' / 'test_config_comprehensive.ini'
        with ini_sample.open('r') as f:
            ini_lines = f.readlines()
        patched_ini = []
        for line in ini_lines:
            if line.strip().startswith('outputDir'):
                patched_ini.append(f'outputDir    = {self.tmp_output}\n')
            else:
                patched_ini.append(line)
        self.tmp_ini = self.tmp_output / 'test_config_comprehensive_patched.ini'
        with self.tmp_ini.open('w') as f:
            f.writelines(patched_ini)

    def tearDown(self):
        # Remove the temp output directory and all its contents
        shutil.rmtree(self.tmp_output, ignore_errors=True)

    def test_configuration_parsing_and_directories(self):
        config = Configuration(str(self.tmp_ini))
        # Check that some expected sections and options are present
        self.assertTrue(hasattr(config, 'globalOptions'))
        self.assertIn('outputDir', config.globalOptions)
        self.assertEqual(config.globalOptions['outputDir'], str(self.tmp_output))
        self.assertTrue(hasattr(config, 'meteoOptions'))
        self.assertIn('precipitationNC', config.meteoOptions)
        # Check that output directories were created
        self.assertTrue(self.tmp_output.is_dir())
        self.assertTrue((self.tmp_output / 'tmp').is_dir())
        self.assertTrue((self.tmp_output / 'netcdf').is_dir())
        self.assertTrue((self.tmp_output / 'scripts').is_dir())
        self.assertTrue((self.tmp_output / 'log').is_dir())
        self.assertTrue((self.tmp_output / 'states').is_dir())
        self.assertTrue((self.tmp_output / 'maps').is_dir())

if __name__ == '__main__':
    unittest.main()
