import argparse
import io
import json
import os
import subprocess
import unittest

def parse_args():
    """Helper function to obtain the binary directory path of Axom from CLI"""
    parser = argparse.ArgumentParser(description="Unit test arguments")
    parser.add_argument("-bd", "--binary-dir", type=str,
                        help="Path to the binary directory for Axom")
    # Add other arguments as needed
    return parser.parse_args()


class TestFortranExampleIntegration(unittest.TestCase):
    
    @classmethod
    def setUpClass(cls):
        """
        Obtain the binary directory from the CLI and compile the sina fortran
        example needed for these tests if necessary.
        """
        cwd = os.getcwd()

        args = parse_args()
        cls.binary_dir = args.binary_dir
        if cls.binary_dir is None:
            # Assume we're at /path/to/build_dir/axom/sina/tests so move up to build_dir
            cls.binary_dir = f"{cwd}/../../../"
        
        os.chdir(cls.binary_dir)

        if not os.path.exists(f"{cls.binary_dir}/examples/sina_fortran_ex"):
            subprocess.run(["make", "sina_fortran_ex"])

        os.chdir(cwd)


    def setUp(self):
        """ Invoke example Fortran application to dump a sina file """
        subprocess.run([f"{self.binary_dir}/examples/sina_fortran_ex"])
        self.dump_file = "sina_dump.json"
        
    def tearDown(self):
        """ Clean up output directory after each test. """
        os.remove(self.dump_file)
    
    def test_file_validity(self):
        """ Make sure the files we're importing follow the Sina schema. """
        try:
            import jsonschema
            schema_file = os.path.join(f"{self.binary_dir}/tests/sina_schema.json")
            with io.open(schema_file, "r", encoding="utf-8") as schema:
                schema = json.load(schema)
                with io.open(self.dump_file, "r", encoding="utf-8") as loaded_test:
                    file_json = json.load(loaded_test)
                    jsonschema.validate(file_json, schema)
        except ModuleNotFoundError:
            print("jsonschema module not found. Skipping test_file_validity.")
            pass
    
                    
    def test_validate_contents_of_record(self):
        """ Ensure that the  record written out matches what we expect """
        with open(self.dump_file, "r", encoding="utf-8") as loaded_test:
            rec = json.load(loaded_test)

        record = rec["records"][0]
        
        # Test the metadata in the record
        self.assertEqual("my_rec_id", record["id"])
        self.assertEqual("my_type", record["type"])
        
        # Test the files
        self.assertEqual(list(record["files"].keys()), ["/path/to/my/file/my_other_file.txt", "/path/to/my/file/my_file.txt"])
        self.assertEqual(record["files"]["/path/to/my/file/my_other_file.txt"]["mimetype"], "png")
        self.assertEqual(record["files"]["/path/to/my/file/my_file.txt"]["mimetype"], "txt")
        
        # Test the signed variants 
        self.assertEqual("A", record["data"]["char"]["value"])
        self.assertEqual(10, record["data"]["int"]["value"])
        self.assertEqual(0, record["data"]["logical"]["value"])
        self.assertEqual(1000000000.0, record["data"]["long"]["value"])
        self.assertEqual(1.23456704616547, record["data"]["real"]["value"])
        self.assertEqual(0.810000002384186, record["data"]["double"]["value"])
        
        # Test the unsigned variants
        self.assertEqual("A", record["data"]["u_char"]["value"])
        self.assertEqual("kg", record["data"]["u_char"]["units"])
        self.assertEqual(10, record["data"]["u_int"]["value"])
        self.assertEqual("kg", record["data"]["u_int"]["units"])
        self.assertEqual(1.0, record["data"]["u_logical"]["value"])
        self.assertEqual("kg", record["data"]["u_logical"]["units"])
        self.assertEqual(1000000000.0, record["data"]["u_long"]["value"])
        self.assertEqual("kg", record["data"]["u_long"]["units"])
        self.assertEqual(1.23456704616547, record["data"]["u_real"]["value"])
        self.assertEqual("kg", record["data"]["u_real"]["units"])
        self.assertEqual(0.810000002384186, record["data"]["u_double"]["value"])
        self.assertEqual("kg", record["data"]["u_double"]["units"])
        self.assertEqual(0.810000002384186, record["data"]["u_double_w_tag"]["value"])
        self.assertEqual("kg", record["data"]["u_double_w_tag"]["units"])
        self.assertEqual(["new_fancy_tag"], record["data"]["u_double_w_tag"]["tags"])
        
        # Test the curves
        nums = range(1, 21)
        real_arr = [i for i in nums]
        double_arr = [i*2 for i in nums]
        int_arr = [i*3 for i in nums]
        long_arr = [i*4 for i in nums]
        curveset = "my_curveset"
        for kind, loc in (("indep", "independent"), ("dep", "dependent")):
            for val_type, target in (("real", real_arr), ("double", double_arr), ("int", int_arr), ("long", long_arr)):
                name = "my_{}_curve_{}".format(kind, val_type)
                self.assertEqual(target, record["curve_sets"][curveset][loc][name]["value"])
        double_2_name = "my_dep_curve_double_2"
        self.assertEqual(double_arr, record["curve_sets"][curveset]["dependent"][double_2_name]["value"])


if __name__ == "__main__":
    # Doing the below instead of unittest.main() so that we can print to stdout
    suite = unittest.TestLoader().loadTestsFromTestCase(TestFortranExampleIntegration)
    runner = unittest.TextTestRunner(buffer=False)
    runner.run(suite)