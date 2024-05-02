import type { SidebarEntry } from '@utils/types';
import type { CollectionEntry } from 'astro:content';

export const createLinkOrGroup = (
    part: string,
    path: string,
    isLastPart: boolean,
    currentUrl: string,
    weight?: number,
) => {
    let entry: SidebarEntry = {
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
    for (let n of Object.values(obj)
        .filter(Boolean)
        .filter((v) => typeof v === 'object')) {
        let found = findByProperty(n, predicate);
        if (found) return n;
    }
};

export const addEntriesToSection = (sections, docs: CollectionEntry<'docs'>[], url: string) => {
    docs.sort((a, b) => a.slug.localeCompare(b.slug));

    docs.forEach((doc) => {
        const parts = doc.slug.split('/');
        let currentLevel = sections;

        parts.forEach((part, i) => {
            part = part.replaceAll('_', ' ')?.replace(/(^)\S/g, (match) => match.toUpperCase());
            // replace nf-core-tools with nf-core/tools
            if (part === 'Nf-core-tools') {
                part = 'nf-core/tools';
            }
            const existingEntry = currentLevel.find((entry) => entry?.label === part);
            const lastPart = i === parts.length - 1;
            const secondToLastPart = i === parts.length - 2;

            if (existingEntry && existingEntry.type === 'link') {
                existingEntry.type = 'group';
                existingEntry.collapsed = true;
                existingEntry.entries = [];
            }

            if (existingEntry && /index\.(md|mdx)$/.test(doc.id) && lastPart) {
                existingEntry.href = '/docs/' + doc.slug;
            }

            if (existingEntry) {
                currentLevel = existingEntry.entries;
                if (secondToLastPart && doc.data.parentWeight) {
                    existingEntry.weight = doc.data.parentWeight;
                }
            } else {
                const newEntry = createLinkOrGroup(
                    lastPart ? doc.data.shortTitle || part : part,
                    lastPart ? '/docs/' + doc.slug : /index\.(md|mdx)$/.test(doc.id) ? '/docs/' + doc.slug : '', // add href to group if they have an index file
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
};

export const sanitizeNfCoreLabels = (label: string) =>
    (label.startsWith('Nf-') ? 'nf-' : '') +
    label
        .substring(label.startsWith('Nf-') ? 3 : 0)
        .replace('Nf-core-tools', 'nf-core/tools')
        .replaceAll('nf-core-tools', 'nf-core/tools')
        .split('nf-')
        .map((part) => part.replaceAll(/-/g, ' '))
        .join(' nf-');
