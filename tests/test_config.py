from hatchet import config


def test_config_hg19():
    assert config.paths.hg19 is not None
