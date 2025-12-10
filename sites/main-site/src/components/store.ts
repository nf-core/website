import { atom } from "nanostores";
import { persistentAtom } from "@nanostores/persistent";

export interface FilterItem {
    name: string;
    displayName?: string;
    class?: string;
    icon?: string;
    count?: number;
}

export const CurrentFilter = atom<{ name: string }[]>([{ name: "" }]);
export const Filters = atom<FilterItem[]>([]);
export const SortBy = atom<string>("");

// Use regular atom during SSR, persistentAtom only in browser
const isBrowser = typeof window !== 'undefined';

export const DisplayStyle = isBrowser
    ? persistentAtom<string>('DisplayStyle', 'grid')
    : atom<string>('grid');

export const Checkboxes = isBrowser
    ? persistentAtom<Array<{ id: string; checked: boolean }>>("Checkboxes", [], {
        decode: (value) => {
            try {
                return JSON.parse(value);
            } catch {
                return [];
            }
        },
    })
    : atom<Array<{ id: string; checked: boolean }>>([]);
export const currentHeading = atom("");
export const currentPage = atom<number>(1);
export const CurrentTab = atom<string>("");
export const SearchQuery = atom<string>("");
export const showHidden = atom<boolean>(false);
export const showHelp = atom<boolean>(false);
