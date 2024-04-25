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

export const addEntriesToSection = (sections, docs: CollectionEntry<'docs'>[], url: str) => {
    docs.sort((a, b) => a.slug.localeCompare(b.slug));

    docs.forEach((doc) => {
        const parts = doc.slug.split('/');
        let currentLevel = sections;

        parts.forEach((part, i) => {
            part = part.replaceAll('_', ' ')?.replace(/(^)\S/g, (match) => match.toUpperCase());
            // replace nf-core-tools with nf-core/tools
            if (part === 'nf-ore-tools') {
                part = 'nf-core/tools';
            }
            const existingEntry = currentLevel.find((entry) => entry?.label === part);
            const lastPart = i === parts.length - 1;
            const secondToLastPart = i === parts.length - 2;
            if (existingEntry) {
                if (existingEntry.type === 'group') {
                    // workaround for index files in nested events, where another element could already have created the group
                } else {
                    existingEntry.type = 'group';
                    existingEntry.collapsed = true;
                    existingEntry.entries = [];
                }
                if (/index\.(md|mdx)$/.test(doc.id) && lastPart) {
                    existingEntry.href = '/docs/' + doc.slug;
                    existingEntry.isCurrent = url === '/docs/' + doc.slug;
                    existingEntry.collapsed = url !== '/docs/' + doc.slug;
                }

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
        .split('-nf-')
        .map((part) => part.replaceAll(/-/g, ' '))
        .join(' nf-');
