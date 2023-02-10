import { persistentAtom } from '@nanostores/persistent';
import { atom } from 'nanostores';


export const CurrentFilter = atom([]);
export const SortBy = atom('');
export const DisplayStyle = atom('');
export const SearchQuery = atom('');
export const showHidden = atom(false);
export const showHelp = atom(false);
export const currentHeading = atom('');
export const EventIsOngoing = persistentAtom('EventIsOngoing', false, {
  encode(value) {
    return JSON.stringify(value);
  },
  decode(value) {
    return JSON.parse(value);
  },
});