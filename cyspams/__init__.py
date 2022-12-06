try:
    from importlib.metadata import version
except ImportError:
    from importlib_metadata import version
__version__ = version('cyspams')

from pathlib import Path

def get_include():
    include_dirs = []
    dir_path = Path(__file__).parent.resolve().joinpath('include')
    include_dirs.append(str(dir_path.joinpath('nnls')))
    include_dirs.append(str(dir_path.joinpath('spams')))
    include_dirs.append(str(dir_path.joinpath('spams', 'decomp')))
    include_dirs.append(str(dir_path.joinpath('spams', 'dictLearn')))
    include_dirs.append(str(dir_path.joinpath('spams', 'linalg')))
    include_dirs.append(str(dir_path.joinpath('spams', 'prox')))
    return include_dirs
