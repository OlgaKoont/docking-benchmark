"""Environment management utilities."""

import subprocess
from pathlib import Path
from typing import Optional, List


def run_in_env(
    command: List[str],
    env_name: Optional[str] = None,
    cwd: Optional[Path] = None,
    check: bool = True,
    capture_output: bool = False,
    timeout: Optional[float] = None
) -> subprocess.CompletedProcess:
    """
    Run command in specified conda environment.
    
    Args:
        command: Command to run.
        env_name: Conda environment name. If None, runs in current environment.
        cwd: Working directory.
        check: Whether to check return code.
        capture_output: Whether to capture output.
        timeout: Maximum time in seconds to wait for command. If None, no timeout.
    
    Returns:
        CompletedProcess result.
    """
    import os
    
    # Prepare environment variables
    env = os.environ.copy()
    
    # If conda environment is specified, set up LD_LIBRARY_PATH for CUDA libraries
    if env_name and env_name != "":
        # Get conda prefix for the environment
        try:
            conda_prefix_result = subprocess.run(
                ["conda", "run", "-n", env_name, "python", "-c", "import os; print(os.environ.get('CONDA_PREFIX', ''))"],
                capture_output=True,
                text=True,
                check=True
            )
            conda_prefix = conda_prefix_result.stdout.strip()
            if conda_prefix:
                # For gnina, use system CUDA libraries to avoid version conflicts
                # Conda CUDA libraries have incompatible symbol versions
                current_ld_path = env.get("LD_LIBRARY_PATH", "")
                system_cuda = "/usr/lib/x86_64-linux-gnu"
                # Use system CUDA libraries first, skip conda CUDA libraries for compatibility
                env["LD_LIBRARY_PATH"] = f"{system_cuda}:{current_ld_path}"
        except Exception:
            pass  # If we can't get conda prefix, continue without modifying LD_LIBRARY_PATH
        
        # Use bash to export LD_LIBRARY_PATH and run command
        # conda run doesn't pass env vars directly, so we use bash -c
        ld_path = env.get('LD_LIBRARY_PATH', '')
        if ld_path:
            # Escape single quotes in LD_LIBRARY_PATH for bash
            ld_path_escaped = ld_path.replace("'", "'\"'\"'")
            # Escape command arguments
            cmd_parts = []
            for arg in command:
                arg_escaped = arg.replace("'", "'\"'\"'")
                cmd_parts.append(f"'{arg_escaped}'")
            cmd_str = " ".join(cmd_parts)
            full_command = [
                "conda", "run", "-n", env_name, "--no-capture-output",
                "bash", "-c",
                f"export LD_LIBRARY_PATH='{ld_path_escaped}' && {cmd_str}"
            ]
        else:
            full_command = ["conda", "run", "-n", env_name, "--no-capture-output"] + command
    else:
        full_command = command
    
    return subprocess.run(
        full_command,
        cwd=cwd,
        check=check,
        capture_output=capture_output,
        text=True,
        timeout=timeout,
        env=env
    )


def check_env_exists(env_name: str) -> bool:
    """Check if conda environment exists."""
    try:
        result = subprocess.run(
            ["conda", "env", "list"],
            capture_output=True,
            text=True,
            check=True
        )
        return env_name in result.stdout
    except (subprocess.CalledProcessError, FileNotFoundError):
        return False


def get_python_in_env(env_name: Optional[str] = None) -> str:
    """Get Python executable path in conda environment."""
    if env_name and env_name != "":
        try:
            result = subprocess.run(
                ["conda", "run", "-n", env_name, "which", "python"],
                capture_output=True,
                text=True,
                check=True
            )
            return result.stdout.strip()
        except subprocess.CalledProcessError:
            return "python3"
    else:
        return "python3"

