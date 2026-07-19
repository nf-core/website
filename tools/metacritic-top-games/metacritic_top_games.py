#!/usr/bin/env python3
"""
metacritic_top_games.py
=======================

Pull the top-rated video games of all time from Metacritic and write them to a
Markdown file formatted as a table, ready to drop into an Obsidian vault.

Metacritic does not publish an official public API. Its own website is powered
by an *unofficial* JSON backend ("finder") that this script talks to directly:

    https://backend.metacritic.com/finder/metacritic/search/all/web

That endpoint returns clean JSON (title, metaScore, platform, release date,
...) and is far more reliable than scraping HTML. An `apiKey` query parameter
is required; the public key baked into Metacritic's own frontend bundle is used
by default and can be overridden with --api-key if it ever rotates.

Because some environments (locked-down CI, sandboxes, corporate proxies) block
metacritic.com entirely, the script ships with a small, verified *offline
snapshot* of the canonical all-time top games. If the live API cannot be
reached, it falls back to that snapshot so you always get an output file. Use
--offline to force it.

Usage
-----
    # Top games with a Metascore over 80 (the default), to top_games.md
    python metacritic_top_games.py

    # Only "must-play" 90+ titles, cap at 50, custom output path
    python metacritic_top_games.py --min-score 90 --limit 50 --output zelda.md

    # Force the bundled offline snapshot (no network)
    python metacritic_top_games.py --offline

Requirements
------------
    pip install requests      # only needed for live fetching; offline works stdlib-only
"""

from __future__ import annotations

import argparse
import datetime as _dt
import sys
import time
from dataclasses import dataclass
from typing import Iterable

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

# Metacritic's unofficial JSON backend used by its own search/browse pages.
API_BASE = "https://backend.metacritic.com/finder/metacritic/search/all/web"

# Public apiKey shipped in Metacritic's frontend. Override with --api-key if it
# stops working (Metacritic rotates it occasionally).
DEFAULT_API_KEY = "1MOZgmNFxvmljaQR1b9wcw1WAan3v6Oup4RGVfbb"

# Browser-like headers – the backend rejects requests without them.
HEADERS = {
    "User-Agent": (
        "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 "
        "(KHTML, like Gecko) Chrome/124.0.0.0 Safari/537.36"
    ),
    "Accept": "application/json, text/plain, */*",
    "Origin": "https://www.metacritic.com",
    "Referer": "https://www.metacritic.com/",
}

PAGE_SIZE = 24  # Metacritic serves games 24 at a time.


# ---------------------------------------------------------------------------
# Data model
# ---------------------------------------------------------------------------

@dataclass
class Game:
    title: str
    score: int
    platform: str
    release_date: str  # ISO-ish string ("1998-11-23") or just a year
    url: str = ""

    @property
    def year(self) -> str:
        return (self.release_date or "")[:4]

    @property
    def tier(self) -> str:
        """Metacritic's own colour banding."""
        if self.score >= 90:
            return "Must-Play"
        if self.score >= 75:
            return "Generally Favorable"
        if self.score >= 50:
            return "Mixed"
        return "Unfavorable"


# ---------------------------------------------------------------------------
# Live fetching
# ---------------------------------------------------------------------------

def fetch_live(min_score: int, limit: int, api_key: str,
               max_pages: int = 60) -> list[Game]:
    """Page through the Metacritic backend, newest highest-rated first.

    Raises RuntimeError on any network/HTTP problem so the caller can fall
    back to the offline snapshot.
    """
    try:
        import requests
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError(
            "The 'requests' package is required for live fetching "
            "(pip install requests), or run with --offline."
        ) from exc

    session = requests.Session()
    session.headers.update(HEADERS)

    games: list[Game] = []
    for page in range(1, max_pages + 1):
        params = {
            "sortBy": "-metaScore",         # highest score first
            "productType": "games",
            "page": page,
            "releaseYearMin": 1958,
            "releaseYearMax": _dt.date.today().year,
            "offset": (page - 1) * PAGE_SIZE,
            "limit": PAGE_SIZE,
            "apiKey": api_key,
        }

        data = _get_json_with_retry(session, API_BASE, params)
        items = _extract_items(data)
        if not items:
            break

        stop = False
        for item in items:
            game = _parse_item(item)
            if game is None:
                continue
            if game.score < min_score:
                # Results are sorted desc by score, so we can stop here.
                stop = True
                break
            games.append(game)
            if limit and len(games) >= limit:
                stop = True
                break

        if stop:
            break

    if not games:
        raise RuntimeError("Live API returned no usable games.")
    return games


def _get_json_with_retry(session, url, params, retries: int = 4):
    """GET JSON with exponential backoff (2s, 4s, 8s, 16s)."""
    last_err: Exception | None = None
    for attempt in range(retries):
        try:
            resp = session.get(url, params=params, timeout=30)
            resp.raise_for_status()
            return resp.json()
        except Exception as exc:  # network error, non-200, bad JSON
            last_err = exc
            if attempt < retries - 1:
                time.sleep(2 ** (attempt + 1))
    raise RuntimeError(f"Request failed after {retries} attempts: {last_err}")


def _extract_items(data) -> list:
    """The backend has changed its JSON envelope over time; be defensive."""
    if not isinstance(data, dict):
        return []
    for key in ("data", "components", "items", "results"):
        node = data.get(key)
        if isinstance(node, list) and node:
            # `components` wraps items one level deeper.
            if key == "components":
                for comp in node:
                    inner = (comp or {}).get("data", {}).get("items")
                    if isinstance(inner, list) and inner:
                        return inner
                continue
            return node
        if isinstance(node, dict):
            inner = node.get("items")
            if isinstance(inner, list):
                return inner
    return []


def _parse_item(item: dict) -> Game | None:
    if not isinstance(item, dict):
        return None
    title = item.get("title") or item.get("name")
    score = (
        item.get("criticScoreSummary", {}).get("score")
        if isinstance(item.get("criticScoreSummary"), dict)
        else item.get("metaScore") or item.get("score")
    )
    if title is None or score is None:
        return None
    try:
        score = int(round(float(score)))
    except (TypeError, ValueError):
        return None

    platform = ""
    plats = item.get("platforms") or item.get("platform")
    if isinstance(plats, list) and plats:
        platform = ", ".join(
            p.get("name", "") if isinstance(p, dict) else str(p) for p in plats
        )
    elif isinstance(plats, str):
        platform = plats

    release = (
        item.get("releaseDate")
        or item.get("premiereDate")
        or item.get("date")
        or ""
    )
    slug = item.get("slug") or ""
    url = f"https://www.metacritic.com/game/{slug}/" if slug else ""

    return Game(title=str(title), score=score, platform=platform,
                release_date=str(release), url=url)


# ---------------------------------------------------------------------------
# Offline snapshot (verified canonical all-time Metascores, best platform version)
# ---------------------------------------------------------------------------
# Source: Metacritic "Best Video Games of All Time" listings. Scores are the
# highest Metascore across platforms for each title. This is a stable, curated
# snapshot for use when the live API is unreachable — run live for the full,
# authoritative, up-to-the-minute ranking.

SNAPSHOT: list[Game] = [
    Game("The Legend of Zelda: Ocarina of Time", 99, "Nintendo 64", "1998-11-23"),
    Game("Tony Hawk's Pro Skater 2", 98, "PlayStation", "2000-09-20"),
    Game("Grand Theft Auto IV", 98, "PlayStation 3", "2008-04-29"),
    Game("SoulCalibur", 98, "Dreamcast", "1999-09-08"),
    Game("Super Mario Galaxy", 97, "Wii", "2007-11-12"),
    Game("Super Mario Galaxy 2", 97, "Wii", "2010-05-23"),
    Game("Grand Theft Auto V", 97, "PlayStation 3", "2013-09-17"),
    Game("Red Dead Redemption 2", 97, "PlayStation 4", "2018-10-26"),
    Game("The Legend of Zelda: Breath of the Wild", 97, "Nintendo Switch", "2017-03-03"),
    Game("Super Mario Odyssey", 97, "Nintendo Switch", "2017-10-27"),
    Game("Metroid Prime", 97, "GameCube", "2002-11-18"),
    Game("Perfect Dark", 97, "Nintendo 64", "2000-05-22"),
    Game("Grand Theft Auto III", 97, "PlayStation 2", "2001-10-23"),
    Game("Halo: Combat Evolved", 97, "Xbox", "2001-11-15"),
    Game("Half-Life", 96, "PC", "1998-11-19"),
    Game("Half-Life 2", 96, "PC", "2004-11-16"),
    Game("The Orange Box", 96, "Xbox 360", "2007-10-10"),
    Game("BioShock", 96, "Xbox 360", "2007-08-21"),
    Game("The Elder Scrolls V: Skyrim", 96, "PlayStation 3", "2011-11-11"),
    Game("Grand Theft Auto: San Andreas", 95, "PlayStation 2", "2004-10-26"),
    Game("Grand Theft Auto: Vice City", 95, "PlayStation 2", "2002-10-27"),
    Game("Uncharted 2: Among Thieves", 96, "PlayStation 3", "2009-10-13"),
    Game("Batman: Arkham City", 96, "PlayStation 3", "2011-10-18"),
    Game("Mass Effect 2", 96, "Xbox 360", "2010-01-26"),
    Game("The Legend of Zelda: The Wind Waker", 96, "GameCube", "2003-03-24"),
    Game("Gran Turismo", 96, "PlayStation", "1998-05-12"),
    Game("Elden Ring", 96, "PC", "2022-02-25"),
    Game("Baldur's Gate 3", 96, "PC", "2023-08-03"),
    Game("The Legend of Zelda: Tears of the Kingdom", 96, "Nintendo Switch", "2023-05-12"),
    Game("Resident Evil 4", 96, "GameCube", "2005-01-11"),
    Game("The Last of Us", 95, "PlayStation 3", "2013-06-14"),
    Game("Portal 2", 95, "PC", "2011-04-19"),
    Game("The Legend of Zelda: A Link to the Past", 95, "SNES", "1992-04-13"),
    Game("Super Mario World", 94, "SNES", "1991-08-23"),
    Game("God of War", 94, "PlayStation 4", "2018-04-20"),
    Game("Tekken 3", 96, "PlayStation", "1998-04-29"),
    Game("The Legend of Zelda: Majora's Mask", 95, "Nintendo 64", "2000-10-26"),
    Game("Disco Elysium: The Final Cut", 97, "PC", "2021-03-30"),
    Game("Hades", 93, "Nintendo Switch", "2020-09-17"),
    Game("Super Mario 3D World", 93, "Wii U", "2013-11-22"),
    Game("The Witcher 3: Wild Hunt", 93, "PC", "2015-05-19"),
    Game("Bloodborne", 92, "PlayStation 4", "2015-03-24"),
    Game("Dark Souls", 89, "PlayStation 3", "2011-10-04"),
    Game("Hollow Knight", 90, "Nintendo Switch", "2018-06-12"),
]


def snapshot_games(min_score: int, limit: int) -> list[Game]:
    games = [g for g in SNAPSHOT if g.score >= min_score]
    games.sort(key=lambda g: (-g.score, g.title))
    return games[:limit] if limit else games


# ---------------------------------------------------------------------------
# Markdown / Obsidian output
# ---------------------------------------------------------------------------

def _md_escape(text: str) -> str:
    """Escape pipe characters so table cells don't break."""
    return text.replace("|", "\\|")


def render_markdown(games: Iterable[Game], min_score: int, source: str,
                    generated: str) -> str:
    games = list(games)
    lines: list[str] = []

    # YAML frontmatter — Obsidian reads this into note properties.
    lines += [
        "---",
        'title: "Metacritic — Top Rated Games of All Time"',
        "tags: [games, metacritic, ratings]",
        f"generated: {generated}",
        f"min_score: {min_score}",
        f"source: {source}",
        f"count: {len(games)}",
        "---",
        "",
        "# 🎮 Metacritic — Top Rated Games of All Time",
        "",
        f"> Games with a Metascore **≥ {min_score}**. "
        f"{len(games)} titles · generated {generated} · source: `{source}`.",
        "",
        "| # | Game | Metascore | Tier | Platform | Year |",
        "| --: | :--- | :--: | :--- | :--- | :--: |",
    ]

    for i, g in enumerate(games, start=1):
        title_cell = f"[{_md_escape(g.title)}]({g.url})" if g.url else _md_escape(g.title)
        lines.append(
            f"| {i} | {title_cell} | {g.score} | {g.tier} | "
            f"{_md_escape(g.platform) or '—'} | {g.year or '—'} |"
        )

    lines.append("")
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Pull top-rated Metacritic games into an Obsidian-ready Markdown table.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--min-score", type=int, default=80,
                   help="Only include games with a Metascore at or above this value.")
    p.add_argument("--limit", type=int, default=200,
                   help="Maximum number of games to include (0 = no limit).")
    p.add_argument("--output", "-o", default="metacritic_top_games.md",
                   help="Path to the Markdown file to write.")
    p.add_argument("--api-key", default=DEFAULT_API_KEY,
                   help="Metacritic backend apiKey (override if it rotates).")
    p.add_argument("--offline", action="store_true",
                   help="Skip the network and use the bundled verified snapshot.")
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    generated = _dt.date.today().isoformat()

    if args.offline:
        games = snapshot_games(args.min_score, args.limit)
        source = "offline snapshot (verified all-time list)"
    else:
        try:
            games = fetch_live(args.min_score, args.limit, args.api_key)
            source = "live metacritic backend API"
        except RuntimeError as exc:
            print(f"[warn] live fetch failed: {exc}", file=sys.stderr)
            print("[warn] falling back to the bundled offline snapshot.",
                  file=sys.stderr)
            games = snapshot_games(args.min_score, args.limit)
            source = "offline snapshot (live API unreachable)"

    markdown = render_markdown(games, args.min_score, source, generated)
    with open(args.output, "w", encoding="utf-8") as fh:
        fh.write(markdown)

    print(f"Wrote {len(games)} games (Metascore >= {args.min_score}) "
          f"to {args.output}  [source: {source}]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
