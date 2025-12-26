import sys
import yaml

if len(sys.argv) < 2:
    print("Usage: python export_paths.py <config_file>")
    sys.exit(1)

config_file = sys.argv[1]

with open(config_file) as f:
    config = yaml.safe_load(f)

print(f"logs/{config['output_prefix']}")
