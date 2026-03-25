import sys

from config_compat import load_config_with_compat

if len(sys.argv) < 2:
    print("Usage: python export_paths.py <config_file>")
    sys.exit(1)

config_file = sys.argv[1]
config = load_config_with_compat(config_file)

print(f"logs/{config['output_prefix']}")
