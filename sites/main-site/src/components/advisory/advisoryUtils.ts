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
