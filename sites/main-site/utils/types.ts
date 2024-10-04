export interface Link {
    type: 'link';
    id: string;
    label: string;
    href: string;
    isCurrent: boolean;
    weight?: number;
    icon?: string;
}

interface Group {
    type: 'group';
    id: string;
    label: string;
    entries: (Link | Group)[];
    collapsed: boolean;
    href?: string;
    isCurrent?: boolean;
    weight?: number;
    icon?: string;
}

export type SidebarEntry = Link | Group;
