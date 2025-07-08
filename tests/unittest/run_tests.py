from subprocess import CalledProcessError, CompletedProcess, run, PIPE

command = [
    "python", "-m", "unittest", "discover", "-s", "tests/unittest/", "-v"]


def print_command_result(result: CompletedProcess[str]) -> None:
    """
    Function that prints out the result of the executed command passed
    as input. If the command return code is negative, the 'CalledProcessError'
    is caught and a 'RuntimeError' is raised providing the 'stderr' output.
    If the exit code is zero, the 'stdout' output is printed to the console,
    if any.

    Parameters
    ----------
    result : CompletedProcess[str]
        Result of the process that has finished running.
    """
    try:
        result.check_returncode()
    except CalledProcessError as e:
        if result.stderr:
            raise RuntimeError(result.stderr)
        raise RuntimeError(e)
    if result.stdout:
        print(result.stdout)
    if result.stderr:
        print(result.stderr)

print_command_result(run(command, stderr=PIPE, stdout=PIPE, text=True))
