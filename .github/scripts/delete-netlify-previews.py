#!/usr/bin/env python3
# /// script
# requires-python = ">=3.12"
# dependencies = ["requests"]
# ///
"""Delete all Netlify deploy previews associated with a merged PR."""

import os
import sys

import requests

NETLIFY_API = "https://api.netlify.com/api/v1"
NETLIFY_TOKEN = os.environ["NETLIFY_TOKEN"]
PR_NUMBER = int(os.environ["PR_NUMBER"])
PER_PAGE = 100

# Netlify sites that deploy previews for the nf-core/website repo
SITE_NAMES = ["nf-co.re", "nf-core-docs"]


def netlify_get(url, params=None):
    r = requests.get(
        url,
        headers={"Authorization": f"Bearer {NETLIFY_TOKEN}"},
        params=params,
        timeout=30,
    )
    r.raise_for_status()
    return r.json()


def netlify_delete(url):
    r = requests.delete(
        url,
        headers={"Authorization": f"Bearer {NETLIFY_TOKEN}"},
        timeout=30,
    )
    r.raise_for_status()


def resolve_site_ids():
    """Map site names to Netlify site IDs."""
    all_sites = netlify_get(
        f"{NETLIFY_API}/sites", {"per_page": PER_PAGE, "filter": "all"}
    )
    site_ids = {}
    for name in SITE_NAMES:
        match = next(
            (
                s
                for s in all_sites
                if s.get("name") == name or s.get("subdomain") == name
            ),
            None,
        )
        if match:
            site_ids[name] = match["id"]
            print(f"  {name} -> {match['id']}")
        else:
            print(f"  {name} -> NOT FOUND, skipping")
    return site_ids


def get_pr_deploys(site_id):
    """Return all deploy-preview deploy IDs for the given PR on a site."""
    deploy_ids = []
    page = 1
    while True:
        deploys = netlify_get(
            f"{NETLIFY_API}/sites/{site_id}/deploys",
            {"per_page": PER_PAGE, "page": page},
        )
        for d in deploys:
            if d.get("context") == "deploy-preview" and d.get("review_id") == PR_NUMBER:
                deploy_ids.append(d["id"])
        if len(deploys) < PER_PAGE:
            break
        page += 1
    return deploy_ids


def main():
    print(f"Deleting Netlify deploy previews for PR #{PR_NUMBER}\n")

    print("Resolving site IDs...")
    site_ids = resolve_site_ids()
    if not site_ids:
        sys.exit("No Netlify sites found")

    total_deleted = 0
    for site_name, site_id in site_ids.items():
        print(f"\n{site_name}:")
        deploy_ids = get_pr_deploys(site_id)
        if not deploy_ids:
            print("  No deploy previews found")
            continue
        print(f"  Found {len(deploy_ids)} deploy preview(s)")
        for deploy_id in deploy_ids:
            print(f"  Deleting {deploy_id}...")
            netlify_delete(f"{NETLIFY_API}/deploys/{deploy_id}")
            total_deleted += 1

    print(f"\nDone. Deleted {total_deleted} deploy preview(s).")


if __name__ == "__main__":
    main()
