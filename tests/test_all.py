"""
test_all: A script for testing the redvals module

Copyright (c) 2025 Haig Bishop

MIT License (see LICENSE file for details)
"""

import pytest
import os
import sys

# Add the parent directory to the path so we can import redvals.py
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from redvals import RedTree, NodeInfo

# Sample paths for testing
BAC_TREE_PATH = "trees/bac120_r220.tree"
ARC_TREE_PATH = "trees/ar53_r220.tree"
BAC_DECORATED_TREE_PATH = "decorated_trees/bac120_r220_decorated.pkl"
ARC_DECORATED_TREE_PATH = "decorated_trees/ar53_r220_decorated.pkl"


class TestRedTree:
    @pytest.fixture
    def skip_if_files_not_exist(self):
        """Skip tests if required tree files don't exist"""
        if not os.path.exists(BAC_TREE_PATH) or not os.path.exists(ARC_TREE_PATH):
            pytest.skip("Tree files not found for testing")
    
    @pytest.fixture
    def skip_if_decorated_files_not_exist(self):
        """Skip tests if decorated tree files don't exist"""
        if not os.path.exists(BAC_DECORATED_TREE_PATH) or not os.path.exists(ARC_DECORATED_TREE_PATH):
            pytest.skip("Decorated tree files not found for testing")
    
    def test_init_from_tree_files(self, skip_if_files_not_exist):
        """Test initialization from tree files"""
        red_trees = RedTree(BAC_TREE_PATH, ARC_TREE_PATH)
        assert red_trees is not None
        assert hasattr(red_trees, 'bac_tree')
        assert hasattr(red_trees, 'arc_tree')
        assert not red_trees.is_decorated()
    
    def test_init_from_decorated_trees(self, skip_if_decorated_files_not_exist):
        """Test initialization from decorated tree files"""
        red_trees = RedTree(BAC_DECORATED_TREE_PATH, ARC_DECORATED_TREE_PATH)
        assert red_trees is not None
        assert hasattr(red_trees, 'bac_tree')
        assert hasattr(red_trees, 'arc_tree')
        assert red_trees.is_decorated()
    
    def test_node_info(self, skip_if_decorated_files_not_exist):
        """Test node info retrieval"""
        red_trees = RedTree(BAC_DECORATED_TREE_PATH, ARC_DECORATED_TREE_PATH)
        # Test with a bacterial node ID (testing with a known leaf node ID)
        node_info = red_trees.get_node_info("bac00000001")
        assert isinstance(node_info, NodeInfo)
        assert node_info.domain == "bac"
        assert node_info.is_terminal == True
        assert 0 <= node_info.red_value <= 1
    
    def test_get_node(self, skip_if_decorated_files_not_exist):
        """Test node retrieval"""
        red_trees = RedTree(BAC_DECORATED_TREE_PATH, ARC_DECORATED_TREE_PATH)
        node = red_trees.get_node("bac00000001")
        assert hasattr(node, 'redvals_id')
        assert node.redvals_id == "bac00000001"
    
    def test_dist_between_nodes(self, skip_if_decorated_files_not_exist):
        """Test distance calculation between nodes"""
        red_trees = RedTree(BAC_DECORATED_TREE_PATH, ARC_DECORATED_TREE_PATH)
        # Test distance between two bacterial nodes
        distance, mrca_id = red_trees.dist_between_nodes("bac00000001", "bac00000002")
        assert isinstance(distance, float)
        assert isinstance(mrca_id, str)
        assert distance >= 0
        assert mrca_id.startswith("bac")
        
        # Test with invalid input
        with pytest.raises(ValueError):
            # Nodes from different domains should raise an error
            red_trees.dist_between_nodes("bac00000001", "arc00000001")


class TestNodeInfo:
    def test_node_info_creation(self):
        """Test NodeInfo class creation and attributes"""
        node_info = NodeInfo(
            domain="bac",
            red_value=0.95,
            gtdb_id="GB_GCA_123456789.1",
            redvals_id="bac00012345",
            is_terminal=True
        )
        
        assert node_info.domain == "bac"
        assert node_info.red_value == 0.95
        assert node_info.gtdb_id == "GB_GCA_123456789.1"
        assert node_info.redvals_id == "bac00012345"
        assert node_info.is_terminal == True
        assert node_info.red_distance == pytest.approx(0.1)
