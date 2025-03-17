import type { SidebarEntry } from '@utils/types';
import type { CollectionEntry } from 'astro:content';

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

