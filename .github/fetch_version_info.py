#!/usr/bin/env -S uv run --script

# /// script
# dependencies = [
#   "requests<3",
# ]
# ///

"""
Fetch version information for nf-core pipelines from GitHub.

## What it does

1. Loads the pipeline list from https://nf-co.re/pipelines.json
2. For each pipeline and release tag:
   - Fetches `nextflow.config` and extracts the `nextflowVersion` value
   - Fetches `.nf-core.yml` and extracts the `nf_core_version` value
   - Stores the `published_at` timestamp from pipelines.json
3. Caches results in `.github/version_info.json` to avoid re-fetching
4. Saves progress after each pipeline to prevent data loss

## Usage

From the root of the repository:

```bash
uv run .github/fetch_version_info.py
```

## Output

The script creates/updates `.github/version_info.json` with a structure like:

```json
{
  "pipelines": {
    "abotyper": {
      "dev": {
        "nextflow_version": "!>=24.04.2",
        "nf_core_version": "3.3.1",
        "published_at": "2025-04-24T05:26:38Z"
      },
      "1.0.0": {
        "nextflow_version": "!>=23.10.0",
        "nf_core_version": "3.0.0",
        "published_at": "2024-12-05T10:30:00Z"
      }
    }
  }
}
```
"""

import json
import re
import requests
from pathlib import Path
from typing import Dict, Any, Optional
import sys


def load_pipelines_json() -> Dict[str, Any]:
    """Load the pipelines.json from nf-co.re"""
    url = "https://nf-co.re/pipelines.json"
    print(f"Fetching pipelines from {url}")

    try:
        response = requests.get(url, timeout=30)
        response.raise_for_status()
        return response.json()
    except Exception as e:
        print(f"Error fetching pipelines.json: {e}")
        sys.exit(1)


def load_cache() -> Dict[str, Any]:
    """Load the local cache file if it exists"""
    cache_path = Path(".github/version_info.json")

    if cache_path.exists():
        try:
            with open(cache_path, 'r') as f:
                return json.load(f)
        except Exception as e:
            print(f"Error loading cache: {e}")
            return {}

    return {}


def save_cache(cache: Dict[str, Any]) -> None:
    """Save the cache to disk"""
    cache_path = Path(".github/version_info.json")
    cache_path.parent.mkdir(parents=True, exist_ok=True)

    try:
        with open(cache_path, 'w') as f:
            json.dump(cache, f, indent=2, sort_keys=True)
    except Exception as e:
        print(f"Error saving cache: {e}")


def fetch_github_file(repo_name: str, tag: str, filename: str) -> Optional[str]:
    """Fetch a file from GitHub for a specific repo and tag"""
    # Use raw.githubusercontent.com for direct file access
    url = f"https://raw.githubusercontent.com/{repo_name}/{tag}/{filename}"

    try:
        response = requests.get(url, timeout=10)
        if response.status_code == 200:
            return response.text
        elif response.status_code == 404:
            print(f"  File {filename} not found for {repo_name}@{tag}")
            return None
        else:
            print(f"  Error fetching {filename}: HTTP {response.status_code}")
            return None
    except Exception as e:
        print(f"  Error fetching {filename}: {e}")
        return None


def extract_nextflow_version(content: str) -> Optional[str]:
    """Extract nextflowVersion from nextflow.config content"""
    if not content:
        return None

    # Look for patterns like: nextflowVersion = '!>=24.04.2'
    # Can have optional spaces and different quote types
    pattern = r'nextflowVersion\s*=\s*[\'"]([^\'\"]+)[\'"]'

    match = re.search(pattern, content)
    if match:
        return match.group(1)

    return None


def extract_nf_core_version(content: str) -> Optional[str]:
    """Extract nf_core_version from .nf-core.yml content"""
    if not content:
        return None

    # Look for patterns like: nf_core_version: 3.3.1
    # Can have optional spaces and quotes
    pattern = r'nf_core_version:\s*[\'"]?([^\s\'\"]+)[\'"]?'

    match = re.search(pattern, content)
    if match:
        return match.group(1)

    return None


def main():
    """Main function"""
    print("Starting version info fetch...")

    # Load pipelines data
    pipelines_data = load_pipelines_json()

    # Load cache
    cache = load_cache()

    # Ensure cache has proper structure
    if 'pipelines' not in cache:
        cache['pipelines'] = {}

    # Process each workflow
    for workflow in pipelines_data.get('remote_workflows', []):
        pipeline_name = workflow['name']
        full_name = workflow['full_name']

        print(f"\nProcessing {pipeline_name}...")

        # Ensure pipeline exists in cache
        if pipeline_name not in cache['pipelines']:
            cache['pipelines'][pipeline_name] = {}

        # Get list of releases
        releases = workflow.get('releases', [])

        # Ensure 'dev' is included
        release_tags = [r['tag_name'] for r in releases]
        if 'dev' not in release_tags:
            releases.append({'tag_name': 'dev'})

        # Process each release
        for release in releases:
            tag = release['tag_name']

            # Check if we already have this info cached
            if tag in cache['pipelines'][pipeline_name]:
                print(f"  {tag}: Already cached, skipping")
                continue

            print(f"  {tag}: Fetching version info...")

            # Initialize entry for this release
            version_info = {
                'nextflow_version': None,
                'nf_core_version': None,
                'published_at': release.get('published_at', None)
            }

            # Fetch nextflow.config
            config_content = fetch_github_file(full_name, tag, 'nextflow.config')
            if config_content:
                version_info['nextflow_version'] = extract_nextflow_version(config_content)

            # Fetch .nf-core.yml
            nf_core_content = fetch_github_file(full_name, tag, '.nf-core.yml')
            if nf_core_content:
                version_info['nf_core_version'] = extract_nf_core_version(nf_core_content)

            # Store in cache
            cache['pipelines'][pipeline_name][tag] = version_info

            # Save cache after each iteration
            save_cache(cache)

            print(f"    Nextflow version: {version_info['nextflow_version'] or 'Not found'}")
            print(f"    nf-core version: {version_info['nf_core_version'] or 'Not found'}")
            print(f"    Published at: {version_info['published_at'] or 'Not found'}")

    print(f"\nCompleted! Processed {len(pipelines_data.get('remote_workflows', []))} pipelines")
    print("Results saved to .github/version_info.json")


if __name__ == "__main__":
    main()
