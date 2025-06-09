import pathlib

import nox

nox.options.error_on_external_run = True
nox.options.default_venv_backend = "uv"


# test every version of python we support!
@nox.session(python=["3.10", "3.11", "3.12"])
def tests(session):
    """Run pytest for all supported python versions

    Test subpackages with:

        $ nox -s test -- pgscatalog.core

    If no subpackage is set, test pgscatalog.utils (
    """
    if session.posargs:
        package = session.posargs[0]
    else:
        package = "pgscatalog.utils"

    if package != "pgscatalog.utils":
        package_path = str(pathlib.Path("packages") / package)
        config_file = str(
            pathlib.Path("packages") / session.posargs[0] / "pyproject.toml"
        )
    else:
        package_path = "."
        config_file = "pyproject.toml"

    session.run_install(
        "uv",
        "sync",
        "--group",
        "test",
        "--package",
        package,
        env={"UV_PROJECT_ENVIRONMENT": session.virtualenv.location},
    )

    # run with --exact to make sure only the subpackage dependencies are loaded
    session.run(
        "uv",
        "run",
        "--group",
        "test",
        "--exact",
        "--package",
        package,
        "pytest",
        "-c",
        config_file,
        package_path,
        external=True,
        env={"UV_PROJECT_ENVIRONMENT": session.virtualenv.location},
    )


@nox.session
def lint(session):
    """Run linting checks"""
    if session.posargs:
        package_path = str(pathlib.Path("packages") / session.posargs[0])
        config_file = str(
            pathlib.Path("packages") / session.posargs[0] / "pyproject.toml"
        )
    else:
        package_path = "src"
        config_file = "pyproject.toml"

    session.run(
        "uv",
        "sync",
        "--group",
        "lint",
        env={"UV_PROJECT_ENVIRONMENT": session.virtualenv.location},
    )

    session.run("ruff", "check", package_path)
    session.run(
        "mypy", package_path, "--config-file", config_file, "--warn-unused-configs"
    )


@nox.session(default=False)
def coverage(session):
    """Generate coverage report"""
    if session.posargs:
        package = session.posargs[0]
        package_path = str(pathlib.Path("packages") / package)
        config_file = str(
            pathlib.Path("packages") / session.posargs[0] / "pyproject.toml"
        )
    else:
        package = "pgscatalog.utils"
        package_path = "."
        config_file = "pyproject.toml"

    session.run_install(
        "uv",
        "sync",
        "--group",
        "test",
        "--package",
        package,
        env={"UV_PROJECT_ENVIRONMENT": session.virtualenv.location},
    )

    # make sure environment is clean
    session.run("coverage", "erase")

    # run coverage <options> -m pytest <path>
    session.run(
        "uv",
        "run",
        "--group",
        "test",
        "--exact",
        "--package",
        package,
        "coverage",
        "run",
        "--source",
        package_path,
        "--rc",
        config_file,
        "-m",
        "pytest",
        package_path,
        external=True,
        env={"UV_PROJECT_ENVIRONMENT": session.virtualenv.location},
    )

    session.run(
        "coverage",
        "combine",
        env={"UV_PROJECT_ENVIRONMENT": session.virtualenv.location},
    )

    session.run("coverage", "xml")


@nox.session(default=False, python="3.12")
def build(session):
    """
    Build source and wheel distributions ready for publishing

    Build subpackages with:

        $ nox -s build -- pgscatalog.core

    If no subpackage is set, builds pgscatalog.utils
    """
    if session.posargs:
        package = session.posargs[0]
    else:
        package = "pgscatalog.utils"

    session.run(
        "uv",
        "sync",
        env={"UV_PROJECT_ENVIRONMENT": session.virtualenv.location},
    )
    session.run("uv", "build", "--package", package)


# nox cookbook: https://nox.thea.codes/en/stable/cookbook.html
# It's a good idea to keep your dev session out of the default list
# so it's not run twice accidentally
@nox.session(default=False)
def dev(session: nox.Session) -> None:
    """
    Set up a python development environment for the project at ".venv".

    Subpackages are configured in the same venv using a uv workspace
    """
    session.run("uv", "venv")
    session.run("uv", "sync", "--all-groups")
