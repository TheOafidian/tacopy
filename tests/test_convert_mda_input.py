import tempfile
import unittest
import os
import tacopy as tp
import pandas as pd
from difflib import unified_diff



class TestConvertMDAMethods(unittest.TestCase):

    def test_convert_mda_formula(self):
        TESTFML = ["~location","~method"]
        tmpfile = tempfile.NamedTemporaryFile(delete=False)
        tmpfile.close()
        tp.write_formula(TESTFML, tmpfile.name)

        with open("tests/data/mda_input/fmls.txt", "r") as f:
            expected_lines = f.readlines()

        with open(tmpfile.name, "r") as f:
            actual_lines = f.readlines()

        os.unlink(tmpfile.name)
        diff = list(unified_diff(expected_lines, actual_lines))
        assert diff == [], "Unexpected file contents:\n" + "".join(diff)


    def test_convert_mda_inputs(self):
        ta = tp.read_tacos("tests/data/urt")
        
        tmpfile = tempfile.NamedTemporaryFile(delete=False)
        tmpfile.close()

        tmpfile2 = tempfile.NamedTemporaryFile(delete=False)
        tmpfile2.close()

        tp.write_mda_inputs(ta, tmpfile.name, tmpfile2.name)
        
        with open("tests/data/mda_input/counts.tsv", "r") as f:
            expected_lines_count = f.readlines()

        with open("tests/data/mda_input/meta.tsv", "r") as f:
            expected_lines_meta = f.readlines()

        with open(tmpfile.name, "r") as f:
            actual_lines_count = f.readlines()

        with open(tmpfile2.name, "r") as f:
            actual_lines_meta = f.readlines()

        os.unlink(tmpfile.name)
        os.unlink(tmpfile2.name)
        diff = list(unified_diff(expected_lines_count, actual_lines_count))
        assert diff == [], "Unexpected file contents:\n" + "".join(diff)

        diff = list(unified_diff(expected_lines_meta, actual_lines_meta))
        assert diff == [], "Unexpected file contents:\n" + "".join(diff)


if __name__ == '__main__':
    unittest.main()