import json
from pathlib import Path

def load_categories(path):
    file = Path(path / "categories.json")
    if not file.exists():
        return{}
    return json.load(open(file))

def save_categories(path, data):
    file = Path(path) / "categories.json"
    with open(file, "w") as f:
        json.dump(data, f, indent=4)