import nox

nox.options.error_on_external_run = True
nox.options.default_venv_backend = "uv"


# test every version of python we support!
@nox.session
@nox.parametrize("python", ["3.12", "3.11", "3.10"])
def tests(session):
    """Run pytest for all supported python versions"""
    session.run_install(
        "uv",
        "sync",
        env={"UV_PROJECT_ENVIRONMENT": session.virtualenv.location},
    )
    session.run("pytest", "--ignore", "noxfile.py")


@nox.session
def coverage(session):
    """Run pytest and output a coverage report"""
    session.run_install(
        "uv",
        "sync",
        env={"UV_PROJECT_ENVIRONMENT": session.virtualenv.location},
    )
    session.run("pytest", "--ignore", "noxfile.py", "--cov", "--cov-report", "xml")


@nox.session
def lint(session):
    """Run linting checks"""
    # https://nox.thea.codes/en/stable/cookbook.html#using-a-lockfile
    # note: uv is set up to install the test and lint dependency groups
    session.run_install(
        "uv",
        "sync",
        env={"UV_PROJECT_ENVIRONMENT": session.virtualenv.location},
    )
    session.run("ruff", "check")
    session.run("mypy", ".")


# nox cookbook: https://nox.thea.codes/en/stable/cookbook.html
# It's a good idea to keep your dev session out of the default list
# so it's not run twice accidentally
@nox.session(default=False)
def dev(session: nox.Session) -> None:
    """
    Set up a python development environment for the project at ".venv".
    """
    session.run("uv", "venv")
    session.run("uv", "sync")
