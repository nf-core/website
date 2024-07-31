import Cache from 'file-system-cache';
import { existsSync, mkdirSync } from 'fs';

let cache;

//create .cache directory if it doesn't exist
if (!existsSync('./.cache')) {
  mkdirSync('./.cache');
}

// need this workaround because in ssr/dev mode vite doesn't use the default export of Cache
// but it does in build/static mode
if (Cache.default !== undefined) {
  cache = Cache.default({
    basePath: './.cache',
    ns: 'nf-core',
  });
} else {
  cache = Cache({
    basePath: './.cache',
    ns: 'nf-core',
  });
}

export default cache;
