import nox

nox.options.error_on_external_run = True
nox.options.default_venv_backend = "uv"


# nox cookbook: https://nox.thea.codes/en/stable/cookbook.html
# It's a good idea to keep your dev session out of the default list
# so it's not run twice accidentally
@nox.session(default=False)
def dev(session: nox.Session) -> None:
    """
    Set up a python development environment for the project at ".venv".
    """
    session.run("uv", "venv")
    session.run(
        "uv",
        "sync",
    )


@nox.session(default=False)
def build(session):
    """
    Build source and wheel distributions ready for publishing
    """
    session.run(
        "uv",
        "sync",
        env={"UV_PROJECT_ENVIRONMENT": session.virtualenv.location},
    )
    session.run("uv", "build")
