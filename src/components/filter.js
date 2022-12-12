import { atom } from 'nanostores';


export const CurrentFilter = atom([]);
export const SortBy = atom('');
export const DisplayStyle = atom('');

CurrentFilter.set(['released', 'under_development', 'archived']);
SortBy.set('alphabetical');
DisplayStyle.set('grid');
