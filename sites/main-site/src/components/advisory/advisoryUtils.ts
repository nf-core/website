import { advisories_type_classes, advisories_type_icons } from "./advisoryTypes";

export function formatAdvisoryType(type: string): string {
    return type.split('_').map(word => word.charAt(0).toUpperCase() + word.slice(1)).join(' ');
}

export function formatAdvisoryCategory(category: string): string {
    return category.charAt(0).toUpperCase() + category.slice(1);
}

export function getAdvisoryTypeIcon(type: string): string {
    return advisories_type_icons[type] || 'fa-circle-info';
}

export function getAdvisoryTypeClass(type: string): string {
    return advisories_type_classes[type] || 'secondary';
}

export function getAdvisorySeverityIcon(severity: string): string {
    return advisories_type_icons[severity] || 'fa-circle-info';
}

export function getAdvisorySeverityClass(severity: string): string {
    return advisories_type_classes[severity] || 'secondary';
}

export function getAdvisoryCategoryIcon(category: string): string {
    switch(category) {
        case 'pipelines': return 'fa-project-diagram';
        case 'modules': return 'fa-object-group';
        case 'subworkflows': return 'fa-sitemap';
        case 'configuration': return 'fa-cogs';
        default: return 'fa-cogs';
    }
}

interface SoftwareDependency {
    name: string;
    versions?: string[];
}

type DependencyItem = string | SoftwareDependency;

export function formatNextflowVersions(versions: string[]): string {
    return versions.join(', ');
}

export function formatNextflowExecutors(executors: string[]): string {
    return executors.join(', ');
}

export function formatSoftwareDependency(dep: DependencyItem): string {
    if (typeof dep === 'string') {
        return dep;
    }
    return `${dep.name}${dep.versions ? ` (${dep.versions.join(', ')})` : ''}`;
}

export function formatSoftwareDependencies(dependencies: DependencyItem[] | string): string {
    if (typeof dependencies === 'string') {
        return dependencies;
    }
    return dependencies.map(formatSoftwareDependency).join(', ');
}

export interface AdvisoryMetadata {
    nextflowVersions?: string[] | null;
    nextflowExecutors?: string[] | null;
    softwareDependencies?: DependencyItem[] | string | null;
}

export interface MetadataItem {
    icon: string;
    label: string;
    value: string;
}

export function getAdvisoryMetadataItems(metadata: AdvisoryMetadata): MetadataItem[] {
    const items: MetadataItem[] = [];

    if (metadata.nextflowVersions?.length) {
        items.push({
            icon: 'fa-code-branch',
            label: 'Nextflow',
            value: formatNextflowVersions(metadata.nextflowVersions)
        });
    }

    if (metadata.nextflowExecutors?.length) {
        items.push({
            icon: 'fa-server',
            label: 'Executors',
            value: formatNextflowExecutors(metadata.nextflowExecutors)
        });
    }

    if (metadata.softwareDependencies) {
        items.push({
            icon: 'fa-object-group',
            label: 'Dependencies',
            value: formatSoftwareDependencies(metadata.softwareDependencies)
        });
    }

    return items;
}
