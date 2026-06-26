#! /usr/bin/env node
// Auto-post newly published blog posts and advisories, and upcoming events,
// to Bluesky and Mastodon. Driven by the /social-feed.json endpoint (built from
// Astro content collections) and a small JSON log that records what has already
// been posted to each network, so nothing is posted twice.
//
// Usage:
//   node bin/social-post.js [--dry-run] [--seed] [--source <url>] [--log <path>]
//
//   --dry-run   Print what would be posted; no network calls, no log changes.
//   --seed      Mark everything currently in the feed as already posted and exit.
//               Run once at rollout so the back-catalogue isn't posted.
//
// Env: BSKY_IDENTIFIER (default "nf-co.re"), BSKY_PASSWORD,
//      MASTODON_URL (default "https://mstdn.science"), MASTODON_ACCESS_TOKEN.

import { readFileSync, writeFileSync } from "fs";
import { AtpAgent, RichText } from "@atproto/api";
import { createRestAPIClient } from "masto";

const args = process.argv.slice(2);
const flag = (name) => args.includes(`--${name}`);
const opt = (name, fallback) => {
    const i = args.indexOf(`--${name}`);
    return i !== -1 && args[i + 1] ? args[i + 1] : fallback;
};

const DRY_RUN = flag("dry-run");
const SEED = flag("seed");
const SOURCE = opt("source", "https://nf-co.re/social-feed.json");
const LOG_PATH = opt("log", "social-posts.json");
const HASHTAGS = "#nfcore #nextflow #bioinformatics #openscience";
const BLUESKY_MAX = 290; // grapheme budget, leaving headroom under the 300 limit
const PRUNE_DAYS = 90; // forget log entries this old — far outside the feed window, so they can't re-post

const ymd = (date) => date.toISOString().slice(0, 10);

// --- decide what is due to post --------------------------------------------

const today = ymd(new Date());
const tomorrow = ymd(new Date(Date.now() + 864e5));

// Returns the lead-in line for an item, or null if it is not due to post.
const lead = (item) => {
    switch (item.type) {
        case "blog":
            return `New blog post: ${item.title}`;
        case "advisory":
            return `New advisory: ${item.title}`;
        case "event":
            if (item.startDate === tomorrow) return `Happening tomorrow: ${item.title}`;
            if (item.startDate === today) return `Happening today: ${item.title}`;
            return null; // events only post the day before / day of
        default:
            return null;
    }
};

// --- message composition ----------------------------------------------------

const blueskyText = (item, leadLine) => {
    // The link is carried by the embed card, so it is left out of the text.
    let text = `${leadLine}\n\n${item.subtitle}`;
    if ([...text].length > BLUESKY_MAX) {
        const room = BLUESKY_MAX - [...leadLine].length - 3; // "\n\n" + ellipsis
        text = `${leadLine}\n\n${[...item.subtitle].slice(0, Math.max(0, room)).join("")}…`;
    }
    return text;
};

const mastodonText = (item, leadLine) => `${leadLine}\n\n${item.subtitle}\n\n${item.url}\n\n${HASHTAGS}`;

// --- network clients (lazy singletons) -------------------------------------

let bsky;
async function blueskyAgent() {
    if (bsky) return bsky;
    bsky = new AtpAgent({ service: "https://bsky.social" });
    await bsky.login({ identifier: process.env.BSKY_IDENTIFIER || "nf-co.re", password: process.env.BSKY_PASSWORD });
    return bsky;
}

async function postBluesky(item, leadLine) {
    const agent = await blueskyAgent();
    const external = { uri: item.url, title: item.title, description: item.subtitle };
    if (item.image) {
        try {
            const res = await fetch(item.image);
            if (res.ok) {
                const blob = new Uint8Array(await res.arrayBuffer());
                const encoding = res.headers.get("content-type") || "image/jpeg";
                const upload = await agent.uploadBlob(blob, { encoding });
                external.thumb = upload.data.blob;
            }
        } catch (e) {
            console.warn(`  ⚠ could not attach Bluesky thumbnail: ${e.message}`);
        }
    }
    const rt = new RichText({ text: blueskyText(item, leadLine) });
    await rt.detectFacets(agent);
    await agent.post({
        text: rt.text,
        facets: rt.facets,
        embed: { $type: "app.bsky.embed.external", external },
        createdAt: new Date().toISOString(),
    });
}

let mastodon;
async function postMastodon(item, leadLine) {
    mastodon ||= createRestAPIClient({
        url: process.env.MASTODON_URL || "https://mstdn.science",
        accessToken: process.env.MASTODON_ACCESS_TOKEN,
    });
    await mastodon.v1.statuses.create({ status: mastodonText(item, leadLine), visibility: "public" });
}

const NETWORKS = {
    bluesky: { post: postBluesky, enabled: () => !!process.env.BSKY_PASSWORD },
    mastodon: { post: postMastodon, enabled: () => !!process.env.MASTODON_ACCESS_TOKEN },
};

// --- main -------------------------------------------------------------------

async function main() {
    let log = {};
    try {
        log = JSON.parse(readFileSync(LOG_PATH, "utf8"));
    } catch (e) {
        if (e.code !== "ENOENT") throw e;
    }
    const save = () => {
        const cutoff = Date.now() - PRUNE_DAYS * 864e5;
        for (const id of Object.keys(log)) {
            const newest = Math.max(...Object.values(log[id]).map((t) => new Date(t).getTime()));
            if (newest < cutoff) delete log[id];
        }
        writeFileSync(LOG_PATH, JSON.stringify(log, null, 2) + "\n");
    };

    const res = await fetch(SOURCE);
    if (!res.ok) throw new Error(`Failed to fetch ${SOURCE}: ${res.status}`);
    const items = await res.json();

    if (SEED) {
        const now = new Date().toISOString();
        for (const item of items) {
            log[item.id] ||= {};
            for (const name of Object.keys(NETWORKS)) log[item.id][name] ||= now;
        }
        save();
        console.log(`Seeded ${items.length} item(s) into ${LOG_PATH} — nothing posted.`);
        return;
    }

    let posted = 0;
    for (const item of items) {
        const leadLine = lead(item);
        if (!leadLine) continue;

        log[item.id] ||= {};
        for (const [name, network] of Object.entries(NETWORKS)) {
            if (log[item.id][name]) continue; // already posted here

            if (DRY_RUN) {
                console.log(`[dry-run] ${name} ← ${leadLine}\n           ${item.url}`);
                continue;
            }
            if (!network.enabled()) continue; // network not configured
            try {
                await network.post(item, leadLine);
                log[item.id][name] = new Date().toISOString();
                save();
                posted++;
                console.log(`✓ ${name}: ${leadLine}`);
            } catch (e) {
                console.error(`✗ ${name} failed for ${item.id}: ${e.message}`);
            }
        }
    }
    console.log(DRY_RUN ? "Dry run complete." : `Done — ${posted} post(s) sent.`);
}

main().catch((e) => {
    console.error(e);
    process.exit(1);
});
