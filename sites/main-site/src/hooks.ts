import { getNextflowVersions } from '../utils/functions';

// This hook runs during build time to populate the Nextflow versions cache
export async function onBuildStart() {
    console.log('Populating Nextflow versions cache during build...');
    try {
        await getNextflowVersions(true); // Force a fresh fetch during build
        console.log('Successfully populated Nextflow versions cache');
    } catch (error) {
        console.error('Failed to populate Nextflow versions cache:', error);
    }
}
