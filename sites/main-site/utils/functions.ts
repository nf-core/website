import type { SidebarEntry } from './types.ts';
import type { CollectionEntry } from 'astro:content';
import { NextflowVersions } from '../src/components/store.ts';
import octokit from '../src/components/octokit';
import type { OctokitResponse } from '@octokit/types';
import type { NextflowVersion } from '../src/components/store.ts';

export const createLinkOrGroup = (
    id: string,
    part: string,
    path: string,
    isLastPart: boolean,
    currentUrl: string,
    weight?: number,
) => {
    let entry: SidebarEntry = {
        id: id,
        label: part,
        href: path,
        isCurrent: currentUrl.replace('.html', '') === path,
        type: isLastPart ? 'link' : 'group',
        weight: weight,
    };

    if (entry.type === 'group') {
        entry.entries = [];
        entry.collapsed = true;
    }

    return entry;
};

export const flattenSidebar = (sidebar: SidebarEntry[]): SidebarEntry[] => {
    if (!sidebar) {
        return [];
    }
    return sidebar.flatMap((entry) =>
        entry.type === 'group' && entry.entries ? flattenSidebar(entry.entries) : entry,
    );
};

export const findCurrentGroup = (sections: SidebarEntry[]) => {
    let currentGroup: SidebarEntry[] = [];
    sections.forEach((section) => {
        if (section.isCurrent) {
            currentGroup = section.type === 'group' ? [section] : sections;
            return currentGroup;
        } else if (section.type === 'group' && section.entries) {
            section.entries.forEach((entry) => {
                if (entry.isCurrent) {
                    currentGroup = section.entries;
                } else if (entry.type === 'group' && entry.entries) {
                    const group = findCurrentGroup(entry.entries);
                    if (group.length > 0) {
                        currentGroup = group;
                    }
                }
            });
        }
    });
    return currentGroup;
};

export const findByProperty = (obj: object, predicate: (obj: object) => boolean) => {
    if (predicate(obj)) return obj;
    for (const n of Object.values(obj)
        .filter(Boolean)
        .filter((v) => typeof v === 'object')) {
        const found = findByProperty(n, predicate);
        if (found) return n;
    }
};
const sortEntries = (entries: SidebarEntry[]) => {
    entries.sort((a, b) => {
        const weightA = a.weight || 1000000;
        const weightB = b.weight || 1000000;

        if (weightA !== weightB) {
            return weightA - weightB;
        } else {
            return a.label.localeCompare(b.label);
        }
    });
    entries.forEach((entry) => {
        if (entry.type === 'group' && entry.entries) {
            sortEntries(entry.entries);
        }
    });
};

export const addEntriesToSection = (sections, docs: CollectionEntry<'docs'>[], url: string, url_prefix: string = '/docs/') => {
    docs.sort((a, b) => a.id.localeCompare(b.id));

    docs.forEach((doc) => {
        const parts = doc.id.replace(/\.[^/.]+$/, '').split('/');
        let currentLevel = sections;

        parts.forEach((part, i) => {
            let label = part.replaceAll('_', ' ')?.replace(/(^)\S/g, (match) => match.toUpperCase());
            // replace nf-core-tools with nf-core/tools
            if (label === 'Nf-core-tools') {
                label = 'nf-core/tools';
            }
            const existingEntry = currentLevel.find((entry) => entry?.label === label);
            const lastPart = i === parts.length - 1;
            const secondToLastPart = i === parts.length - 2;

            if (existingEntry && existingEntry.type === 'link') {
                existingEntry.type = 'group';
                existingEntry.collapsed = true;
                existingEntry.entries = [];
            }

            if (existingEntry) {
                currentLevel = existingEntry.entries;
                if (secondToLastPart && doc.data.parentWeight) {
                    existingEntry.weight = doc.data.parentWeight;
                }
            } else {
                const newEntry = createLinkOrGroup(
                    parts.slice(0, i + 1).join('_'),
                    lastPart ? doc.data.shortTitle || label : label,
                    lastPart ? url_prefix + doc.id.replace(/\.[^/.]+$/, '') : '', // add href to group if they have an index file
                    lastPart,
                    url,
                    secondToLastPart && doc.data.parentWeight ? doc.data.parentWeight : doc.data.weight,
                );
                currentLevel.push(newEntry);

                if (newEntry.type === 'group' && newEntry.entries) {
                    currentLevel = newEntry.entries;
                }
            }
        });
    });

    sortEntries(sections);
};
export const getSectionParents = (sections: SidebarEntry[], currentSection: SidebarEntry) => {
    let parents: SidebarEntry[] = [];
    findByProperty(sections, (entry) => {
        if (entry.type === 'group' && entry.entries) {
            const found = entry.entries.find((subEntry) => subEntry.id === currentSection.id);
            if (found) {
                parents.push(entry);
                parents = [getSectionParents(sections, entry).flat(), ...parents].flat();
            }
        }
    });

    return parents;
};

export const sanitizeNfCoreLabels = (label: string) => {
    return (
        (label.startsWith('Nf-') ? 'nf-' : '') +
        label
            .substring(label.startsWith('Nf-') ? 3 : 0)
            .replace('Nf-core-tools', 'nf-core/tools')
            .replaceAll('nf-core-tools', 'nf-core/tools')
            .split('nf-')
            .map((part) => part.replaceAll(/-/g, ' '))
            .join(' nf-')
    );
};

// Helper function to check if cache is older than 24 hours
const isCacheExpired = (lastUpdated: number): boolean => {
    return Date.now() - lastUpdated > 24 * 60 * 60 * 1000;
};

/**
 * Synchronously get the cached Nextflow versions.
 * This is safe to use in components and during build time.
 * Returns empty array if no versions are cached yet.
 */
export function getCachedNextflowVersions(): NextflowVersion[] {
    const cached = NextflowVersions.get();
    return cached.versions;
}

/**
 * Asynchronously fetch Nextflow versions from GitHub and update the cache.
 * During build time, this will populate the cache for synchronous access.
 */
export async function getNextflowVersions(renew = false): Promise<NextflowVersion[]> {
    const cached = NextflowVersions.get();
    const cacheExpired = isCacheExpired(cached.lastUpdated);

    // Use cache if it exists, is not expired, and renewal is not requested
    if (cached.versions.length > 0 && !cacheExpired && !renew) {
        console.log("Using cached Nextflow versions");
        return cached.versions;
    }

    // If cache is expired or renewal requested, fetch from GitHub
    console.log("Fetching Nextflow versions from GitHub");
    const versions: {}[] = [];
    let page = 1;
    let releases: OctokitResponse<any>;
    do {
        releases = await octokit.request('GET /repos/{owner}/{repo}/releases', {
            owner: 'nextflow-io',
            repo: 'nextflow',
            page,
        });
        versions.push(...releases.data);
        page++;
    } while (releases.headers.link?.includes('rel="next"'));

    const formattedVersions = versions.map((version: any) => ({
        version: version['tag_name'],
        isEdge: version['prerelease'],
        downloadUrl: version['assets'][0] && version['assets'][0]['browser_download_url'],
        downloadUrlAll: version['assets'][1] && version['assets'][1]['browser_download_url'],
        published_at: version['published_at'],
    })).sort((a: any, b: any) => {
        const aVersion = a['version'].replace('-edge', '').replace('v', '');
        const bVersion = b['version'].replace('-edge', '').replace('v', '');
        const aVersionSplit = aVersion.split('.');
        const bVersionSplit = bVersion.split('.');
        for (let i = 0; i < 3; i++) {
            if (parseInt(aVersionSplit[i]) > parseInt(bVersionSplit[i])) {
                return -1;
            } else if (parseInt(aVersionSplit[i]) < parseInt(bVersionSplit[i])) {
                return 1;
            }
        }
        return 0;
    });

    // Update the store with versions and current timestamp
    NextflowVersions.set({
        versions: formattedVersions,
        lastUpdated: Date.now()
    });

    return formattedVersions;
}

