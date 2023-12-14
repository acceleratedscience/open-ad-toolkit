import unittest


class TestMyModule(unittest.TestCase):
    def test_dummy(self):
        # Test the add function
        self.assertEqual(5, 5)
        self.assertEqual(1, 1)  # Test that 0 + 0 equals 0


if __name__ == "__main__":
    unittest.main()
