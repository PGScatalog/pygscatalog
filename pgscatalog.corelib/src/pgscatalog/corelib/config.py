import importlib.metadata

_package_version = f"{importlib.metadata.version('pgscatalog.corelib')}"
_package_string = f"pgscatalog.corelib/{_package_version}"

API_HEADER = {"user-agent": _package_string}
