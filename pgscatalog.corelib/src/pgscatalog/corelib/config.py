import importlib.metadata
import pathlib

_package_version = f"{importlib.metadata.version('pgscatalog.corelib')}"
_package_string = f"pgscatalog.corelib/{_package_version}"

API_HEADER = {"user-agent": _package_string}
ROOT_DIR = pathlib.Path(__file__).resolve().parent.parent.parent.parent
