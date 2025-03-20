import { persistentAtom } from '@nanostores/persistent';
import { atom } from 'nanostores';

export const CurrentFilter = atom([{ name: '' }]);
export const Filters = atom([]);
export const SortBy = atom('');
// add persistentatom for display style to be either table or grid

export const DisplayStyle = persistentAtom('DisplayStyle', 'grid', {
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
export const currentPage = atom(1);
export const CurrentTab = atom('');
export const SearchQuery = atom('');
export const showHidden = atom(false);
export const showHelp = atom(false);
