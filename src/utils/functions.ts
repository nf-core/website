import type { SidebarEntry } from '@utils/types';
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
