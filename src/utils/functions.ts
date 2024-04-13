import type { SidebarEntry } from '@utils/types';
export const createLinkOrGroup = (part, path, isLastPart, currentUrl) => {
    let entry: SidebarEntry = {
        label: part,
        href: path,
        isCurrent: currentUrl.replace('.html', '') === path,
        type: isLastPart ? 'link' : 'group',
    };

    if (entry.type === 'group') {
        entry.entries = [];
        entry.collapsed = true; // You can set this dynamically as well
    }

    return entry;
};

const flattenSidebar = (sidebar: SidebarEntry[]): SidebarEntry[] => {
    if (!sidebar) {
        return [];
    }
    return sidebar.flatMap((entry) =>
        entry.type === 'group' && entry.entries ? flattenSidebar(entry.entries) : entry,
    );
};
