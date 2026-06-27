from pathlib import Path

ROOT = Path("data")

def ensure_workspace(name):
    path = ROOT / name
    path.mkdir(parents=True, exist_ok=True)
    return path
def ensure_project(workspace, project):
    path = ROOT / workspace / project / "files"
    path.mkdir(parents=True, exist_ok=True)
    return path