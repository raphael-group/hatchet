import multiprocessing as mp
from hatchet import config


def test_config_hg19():
    assert config.paths.hg19 is not None


def test_cpus():
    assert mp.cpu_count() == 12
