"""
Module containing unittest classes to assess that the functions in the
`glow.geometry_layouts.fillable_layouts` module have a valid implementation.
"""
import unittest

from glow.geometry_layouts.cells import CartesianCell
from glow.geometry_layouts.fillable_layouts import find_region_path_in_tree, \
    follow_path_to_node, get_region_in_tree, print_subtree
from glow.geometry_layouts.geometries import Circle, Rectangle
from glow.geometry_layouts.layouts import Region
from glow.interface.geom_interface import ShapeType
from glow.support.types import PropertyType
from glow.support.utility import are_same_shapes
from tests.unittest.support_funcs import capture_output


class TestFindRegionPathInTree(unittest.TestCase):
    """
    Test case for verifying the implementation and behaviour of the function
    `find_region_path_in_tree` in the `fillable_layouts.py` module.
    """
    def test_find_region_in_single_layer(self) -> None:
        """
        Method that verifies finding a region path in a simple hierarchical
        tree with a single layer containing regions.
        """
        # Create a simple tree with one layer containing regions
        region1 = Region(Circle(), properties={PropertyType.MATERIAL: "MAT1"})
        region2 = Region(
            Circle(radius=2), properties={PropertyType.MATERIAL: "MAT2"}
        )
        layers = [[region1, region2]]

        # Find the path to the first region
        path = find_region_path_in_tree(layers, region1.geom_obj)
        self.assertIsNotNone(path)
        self.assertEqual(path, [(0, 0)])

        # Find the path to the second region
        path = find_region_path_in_tree(layers, region2.geom_obj)
        self.assertIsNotNone(path)
        self.assertEqual(path, [(0, 1)])

    def test_find_region_in_nested_fillable(self) -> None:
        """
        Method that tests finding a region path in a hierarchical tree
        containing nested `Fillable` objects.
        """
        # Create nested fillable structure
        inner_cell = CartesianCell()
        inner_region = Region(
            Circle(), properties={PropertyType.MATERIAL: "MAT1"}
        )
        inner_cell.add(inner_region)
        inner_cell.update_hierarchical_structure()

        outer_cell = CartesianCell(width_height=(4, 4))
        outer_cell.add(inner_cell)
        outer_cell.update_hierarchical_structure()

        # Find the path to the nested region
        path = find_region_path_in_tree(
            outer_cell.layers, inner_region.geom_obj
        )
        self.assertIsNotNone(path)
        # Path should have two elements: one for outer layer, one for inner
        self.assertEqual(len(path), 2)
        self.assertEqual(path, [(1, 0), (1, 0)])

    def test_find_region_not_in_tree(self) -> None:
        """
        Method that tests that `None` is returned when searching for a region
        not present in the tree.
        """
        # Create a tree with one region
        region1 = Region(Circle(), properties={PropertyType.MATERIAL: "MAT1"})
        layers = [[region1]]

        # Create a different region not in the tree
        region_not_in_tree = Region(
            Circle(radius=5), properties={PropertyType.MATERIAL: "MAT2"}
        )

        # Try to find the region not in the tree
        path = find_region_path_in_tree(layers, region_not_in_tree.geom_obj)
        self.assertIsNone(path)

    def test_find_region_in_multiple_layers(self) -> None:
        """
        Method that tests finding regions in a tree with multiple layers.
        """
        # Create a tree with two layers
        region1 = Region(Circle(), properties={PropertyType.MATERIAL: "MAT1"})
        region2 = Region(
            Circle(radius=1.5), properties={PropertyType.MATERIAL: "MAT2"}
        )
        region3 = Region(
            Circle(radius=2.5), properties={PropertyType.MATERIAL: "MAT3"}
        )
        layers = [[region1, region2], [region3]]

        # Find regions in different layers
        path1 = find_region_path_in_tree(layers, region1.geom_obj)
        self.assertEqual(path1, [(0, 0)])
        path2 = find_region_path_in_tree(layers, region2.geom_obj)
        self.assertEqual(path2, [(0, 1)])
        path3 = find_region_path_in_tree(layers, region3.geom_obj)
        self.assertEqual(path3, [(1, 0)])

    def test_find_region_empty_tree(self) -> None:
        """
        Method that tests finding a region in an empty hierarchical tree.
        """
        # Create a region to search for
        region = Region(Circle(), properties={PropertyType.MATERIAL: "MAT1"})

        # Try to find the region in the empty tree
        path = find_region_path_in_tree([[]], region.geom_obj)
        self.assertIsNone(path)

    def test_initial_path_parameter(self) -> None:
        """
        Method that tests the behaviour of the function when providing the
        initial `path` parameter.
        """
        # Create nested fillable structure
        inner_cell = CartesianCell()
        inner_region = Region(
            Circle(), properties={PropertyType.MATERIAL: "MAT1"}
        )
        inner_cell.add(inner_region)
        inner_cell.update_hierarchical_structure()

        outer_cell = CartesianCell(width_height=(4, 4))
        outer_cell.add(inner_cell)
        outer_cell.update_hierarchical_structure()

        # Find with an initial path
        initial_path = [(0, 0)]
        path = find_region_path_in_tree(
            outer_cell.layers, inner_region.geom_obj, path=initial_path
        )
        self.assertIsNotNone(path)
        # The returned path should include the initial path
        self.assertEqual(path[0], (0, 0))


class TestFollowPathToNode(unittest.TestCase):
    """
    Test case for verifying the implementation and behaviour of the function
    `follow_path_to_node` in the `fillable_layouts.py` module.

    Attributes
    ----------
    fillable : CartesianCell
        The `Fillable` subclass object used in the tests.
    """
    def setUp(self) -> None:
        self.fillable: CartesianCell = CartesianCell(width_height=(2, 2))

    def test_follow_path_to_node_single_level(self) -> None:
        """
        Method that tests the `follow_path_to_node` function when following
        a path of a single level in the hierarchy.
        """
        # Add a region to the fillable
        region = Region(Circle(), properties={PropertyType.MATERIAL: "MAT"})
        self.fillable.add(region)
        # Build a path to the region
        path = [(1, 0)]
        # Follow the path to get the node
        result = follow_path_to_node(self.fillable, path)
        # Verify the result is the region
        self.assertIsInstance(result, Region)
        self.assertTrue(are_same_shapes(result, region, ShapeType.FACE))

    def test_follow_path_to_node_multi_level(self) -> None:
        """
        Method that tests the `follow_path_to_node` function when following
        a path with multiple levels in the hierarchy.
        """
        # Create a nested fillable structure
        nested_cell = CartesianCell()
        region = Region(Circle(), properties={PropertyType.MATERIAL: "MAT"})
        nested_cell.add(region)
        # Add the nested cell to the main fillable
        self.fillable.add(nested_cell)
        # Build a path to the region in the nested structure
        path = [(1, 0), (1, 0)]

        # Follow the path to get the node
        result = follow_path_to_node(self.fillable, path)
        # Verify the result is the region
        self.assertIsInstance(result, Region)
        self.assertTrue(are_same_shapes(result, region, ShapeType.FACE))

    def test_follow_path_to_node_empty_path(self) -> None:
        """
        Method that tests the `follow_path_to_node` function when following
        an empty path, which returns the root node itself.
        """
        # Follow the path to get the node
        result = follow_path_to_node(self.fillable, [])
        # Verify the result is the root fillable itself
        self.assertIs(result, self.fillable)

    def test_follow_path_to_node_invalid_layer_index(self) -> None:
        """
        Method that tests the `follow_path_to_node` function when following
        a path with an invalid layer index.
        """
        # Add a region to the fillable
        region = Region(Circle(), properties={PropertyType.MATERIAL: "MAT"})
        self.fillable.add(region)
        # Verify the exception is raised when a path with an invalid layer
        # index is provided
        with self.assertRaises(IndexError):
            follow_path_to_node(self.fillable, [(5, 0)])

class TestGetRegionInTree(unittest.TestCase):
    """
    Test case for verifying the implementation and behaviour of the function
    `get_region_in_tree` in the `fillable_layouts.py` module.
    """
    def test_simple_tree(self) -> None:
        """
        Method that tests `get_region_in_tree` with a simple tree containing
        a single layer with three regions.
        """
        # Create multiple regions in the same layer
        region1 = Region(
            Circle(radius=1), properties={PropertyType.MATERIAL: "MAT1"}
        )
        region2 = Region(
            Circle(radius=2), properties={PropertyType.MATERIAL: "MAT2"}
        )
        region3 = Region(
            Circle(radius=3), properties={PropertyType.MATERIAL: "MAT3"}
        )
        tree = [[region1, region2, region3]]

        # Get the second region by its shape
        found_region = get_region_in_tree(tree, region2.geom_obj)

        # Verify the correct region is found
        self.assertIsNotNone(found_region)
        self.assertIs(found_region, region2)

    def test_multi_layer_tree(self) -> None:
        """
        Method that tests `get_region_in_tree` with a tree containing
        multiple layers with multiple regions each.
        """
        # Create regions across multiple layers
        region1 = Region(
            Circle(radius=1), properties={PropertyType.MATERIAL: "MAT1"}
        )
        region2 = Region(
            Circle(radius=2), properties={PropertyType.MATERIAL: "MAT2"}
        )
        region3 = Region(
            Circle(radius=3), properties={PropertyType.MATERIAL: "MAT3"}
        )
        region4 = Region(
            Circle(radius=4), properties={PropertyType.MATERIAL: "MAT4"}
        )
        tree = [[region1, region2], [region3, region4]]

        # Get a region from the second layer
        found_region = get_region_in_tree(tree, region4.geom_obj)

        # Verify the correct region is found
        self.assertIsNotNone(found_region)
        self.assertIs(found_region, region4)

    def test_nested_fillable_tree(self) -> None:
        """
        Method that tests `get_region_in_tree` with a nested tree containing
        `Fillable` objects.
        """
        # Create nested fillable with region
        nested_region = Region(Circle(radius=0.5))
        nested_fillable = CartesianCell()
        nested_fillable.add(nested_region)
        nested_fillable.update_hierarchical_structure()
        # Create root fillable with nested fillable
        root_region = Region(Circle((2.5, 2.5, 0.0), radius=0.2))
        root_fillable = CartesianCell(width_height=(3, 3))
        root_fillable.add(root_region)
        root_fillable.add(nested_fillable)
        root_fillable.update_hierarchical_structure()

        # Get the region in the nested hierarchy
        found_region = get_region_in_tree(
            root_fillable.layers, nested_region.geom_obj
        )
        # Verify the nested region is found by checking if they are the same
        # shape
        self.assertIsNotNone(found_region)
        self.assertTrue(
            are_same_shapes(found_region, nested_region, ShapeType.FACE)
        )

        # Get the root region
        found_region = get_region_in_tree(
            root_fillable.layers, root_region.geom_obj
        )
        # Verify the root region is found by checking if they are the same
        # shape
        self.assertIsNotNone(found_region)
        self.assertTrue(
            are_same_shapes(found_region, root_region, ShapeType.FACE)
        )

    def test_no_matching_region(self) -> None:
        """
        Method that tests `get_region_in_tree` when no region matches the
        given shape, and `None` is returned.
        """
        # Create a tree
        region1 = Region(
            Circle(radius=1.0), properties={PropertyType.MATERIAL: "MAT1"}
        )
        region2 = Region(
            Circle(radius=2.0), properties={PropertyType.MATERIAL: "MAT2"}
        )
        tree = [[region1, region2]]
        # Create a shape that is not in the tree
        non_matching_shape = Circle((3.0, 3.0, 0.0), 5.0)

        # Attempt to find the non-matching region
        found_region = get_region_in_tree(tree, non_matching_shape)
        # Verify 'None' is returned
        self.assertIsNone(found_region)

    def test_empty_tree(self) -> None:
        """
        Method that tests `get_region_in_tree` with an empty tree or with
        empty layers.
        """
        # Create an empty tree
        tree = []
        region = Region(Circle(), properties={PropertyType.MATERIAL: "MAT1"})

        # Attempt to find a region in an empty tree
        found_region = get_region_in_tree(tree, region.geom_obj)
        # Verify None is returned
        self.assertIsNone(found_region)

        # Create a tree with empty layer
        region = Region(Circle(), properties={PropertyType.MATERIAL: "MAT1"})
        tree = [[], [region], []]
        # Get the region
        found_region = get_region_in_tree(tree, region.geom_obj)

        # Verify the correct region is found
        self.assertIsNotNone(found_region)
        self.assertIs(found_region, region)

    def test_point_inside_shape_matching(self) -> None:
        """
        Method that tests `get_region_in_tree` when a vertex built within
        the input shape matches a region.
        """
        # Create a region with a specific geometry
        region = Region(Circle(radius=2.0, center=(0.0, 0.0, 0.0)))
        tree = [[region]]
        # Create a shape that is half the region
        shape = region - Rectangle((1.0, 0.0, 0.0), 2, 1)

        # Get the region using point-inside-shape logic
        found_region = get_region_in_tree(tree, shape)

        # Verify the region is found based on point inclusion
        self.assertIsNotNone(found_region)
        self.assertIs(found_region, region)

class TestPrintSubtree(unittest.TestCase):
    """
    Test case for verifying the implementation and behaviour of the function
    `print_subtree` in the `fillable_layouts.py` module.

    Attributes
    ----------
    fillable : CartesianCell
        The `Fillable` subclass object used in the tests.
    """
    def setUp(self) -> None:
        self.fillable: CartesianCell = CartesianCell(width_height=(2, 2))

    def test_print_subtree(self) -> None:
        """
        Method that tests the `print_subtree` function for printing the tree
        structure from a root `Fillable` to a target `Region` or `Fillable`.
        """
        # Create a nested structure with regions
        cell = CartesianCell()
        region = Region(
            Circle(), "TestRegion", {PropertyType.MATERIAL: "MAT1"}
        )
        cell.add(region)
        self.fillable.add(cell)
        # Update the hierarchical structure
        self.fillable.update_hierarchical_structure()

        # Get the path to the region
        path = find_region_path_in_tree(self.fillable.layers, region.geom_obj)
        self.assertIsNotNone(path)

        # Capture the output
        captured = capture_output(
            print_subtree, self.fillable, region, path
        )
        # Verify the root is printed
        self.assertIn(f"{self.fillable.name} (root)", captured)
        # Verify the target region name is in the output
        self.assertIn(f"'{region.name}'", captured)
        # Verify the access code is printed
        self.assertIn("Access Code:", captured)
        # Verify layer indices are in the access code
        for layer in path:
            self.assertIn(f"layers[{layer[0]}][{layer[1]}]", captured)

    def test_print_subtree_nested_fillables(self) -> None:
        """
        Method that tests the `print_subtree` function with deeply nested
        Fillable objects.
        """
        # Create a nested structure
        cell1 = CartesianCell(width_height=(0.5, 0.5))
        cell2 = CartesianCell(width_height=(0.25, 0.25))
        region = Region(Circle(radius=0.1))
        cell2.add(region)
        cell1.add(cell2)
        self.fillable.add(cell1)
        self.fillable.update_hierarchical_structure()

        # Get the path to the region
        path = find_region_path_in_tree(self.fillable.layers, region.geom_obj)
        self.assertIsNotNone(path)

        # Capture the output
        captured = capture_output(print_subtree, self.fillable, region, path)
        # Verify the entire tree path is printed
        self.assertIn(f"{self.fillable.name} (root)", captured)
        self.assertIn(f"{cell1.name}", captured)
        self.assertIn(f"{cell2.name}", captured)
        self.assertIn(f"'{region.name}'", captured)
        # Verify multiple layers are shown in the access code
        layer_count = captured.count("layers[")
        self.assertGreaterEqual(layer_count, len(path))
        # Verify layer indices are in the access code
        for layer in path:
            self.assertIn(f"layers[{layer[0]}][{layer[1]}]", captured)


if __name__ == "__main__":
    unittest.main()
