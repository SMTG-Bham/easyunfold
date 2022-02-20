"""
Pytest fixtures
"""

import pytest
from pathlib import Path


@pytest.fixture
def datapath():
    """Return the path relative to the test_data folder"""

    def _inner(rel):

        return Path(__file__).parent / 'test_data' / rel

    return _inner
