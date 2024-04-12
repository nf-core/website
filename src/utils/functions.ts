import type { SidebarEntry } from '@utils/types';
export const createLinkOrGroup = (part, path, isLastPart, currentUrl) => {
    let entry: SidebarEntry = {
        label: part,
        href: path,
        isCurrent: currentUrl === path,
        type: isLastPart ? 'link' : 'group',
    };

    if (entry.type === 'group') {
        entry.entries = [];
        entry.collapsed = true; // You can set this dynamically as well
    }

    return entry;
};
