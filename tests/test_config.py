from hatchet import config


def test_config_reference():
    assert config.paths.reference is not None
