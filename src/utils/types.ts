export interface Link {
    type: 'link';
    label: string;
    href: string;
    isCurrent: boolean;
}

interface Group {
    type: 'group';
    label: string;
    entries: (Link | Group)[];
    collapsed: boolean;
    href?: string;
}

export type SidebarEntry = Link | Group;
