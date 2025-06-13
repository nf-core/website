import { persistentAtom } from '@nanostores/persistent';
import { atom } from 'nanostores';

export interface FilterItem {
    name: string;
    displayName?: string;
    class?: string;
    icon?: string;
    count?: number;
}

export const CurrentFilter = atom<{ name: string }[]>([{ name: '' }]);
export const Filters = atom<FilterItem[]>([]);
export const SortBy = atom<string>('');
// add persistentatom for display style to be either table or grid

export const DisplayStyle = persistentAtom<string>('DisplayStyle', 'grid', {
    encode(value) {
        return JSON.stringify(value);
    },
    decode(value) {
        return JSON.parse(value);
    },
});

export const Checkboxes = persistentAtom<Array<{id: string, checked: boolean}>>('Checkboxes', [], {
    encode(value) {
        return JSON.stringify(value);
    },
    decode(value) {
        try {
            return JSON.parse(value);
        } catch {
            return [];
        }
    },
});
export const currentHeading = atom('');
export const currentPage = atom<number>(1);
export const CurrentTab = atom<string>('');
export const SearchQuery = atom<string>('');
export const showHidden = atom<boolean>(false);
export const showHelp = atom<boolean>(false);

export interface NextflowVersion {
    version: string;
    isEdge: boolean;
    downloadUrl: string;
    downloadUrlAll: string;
    published_at: string;
}

interface CachedVersions {
    versions: NextflowVersion[];
    lastUpdated: number; // timestamp in milliseconds
}

export const NextflowVersions = atom<CachedVersions>({ versions: [], lastUpdated: 0 });

// Helper function to check if cache is expired
export const isCacheExpired = (lastUpdated: number): boolean => {
    // Cache expiration time (24 hours in milliseconds)
    return Date.now() - lastUpdated > 24 * 60 * 60 * 1000;
};
