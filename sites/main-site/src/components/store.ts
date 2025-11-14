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

export const DisplayStyle = persistentAtom<string>('DisplayStyle', 'grid', {
    encode: JSON.stringify,
    decode(value) {
        try {
            return JSON.parse(value);
        } catch {
            return 'grid';
        }
    },
});

export const Checkboxes = persistentAtom<Array<{ id: string; checked: boolean }>>("Checkboxes", [], {
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
export const currentHeading = atom("");
export const currentPage = atom<number>(1);
export const CurrentTab = atom<string>("");
export const SearchQuery = atom<string>("");
export const showHidden = atom<boolean>(false);
export const showHelp = atom<boolean>(false);
