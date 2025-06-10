export const advisories_type_classes = {
    // Advisory types
    known_regression: "dark",
    incompatibility: "warning",
    security: "danger",
    performance: "success",
    data_corruption: "primary",
    scientific_advice: "info",
    other: "secondary",
    // Severity levels
    low: "primary",
    medium: "warning",
    high: "danger",
    critical: "danger",
} as const;

export const advisories_type_icons = {
    // Advisory types
    known_regression: "fa-solid fa-bug",
    incompatibility: "fa-solid fa-exclamation-triangle",
    security: "fa-solid fa-shield-alt",
    performance: "fa-solid fa-tachometer-alt",
    data_corruption: "fa-solid fa-database",
    scientific_advice: "fa-solid fa-flask",
    other: "fa-solid fa-info-circle",
    // Severity levels
    low: "fa-solid fa-arrow-down",
    medium: "fa-solid fa-minus",
    high: "fa-solid fa-arrow-up",
    critical: "fa-solid fa-exclamation-circle",
} as const;

// Only include actual advisory types, not severity levels
const advisory_types = [
    "known_regression",
    "incompatibility",
    "security",
    "performance",
    "data_corruption",
    "scientific_advice",
    "other",
] as const;

export const advisories_types = advisory_types.map((type) => {
    return {
        name: type,
        class: advisories_type_classes[type],
        icon: advisories_type_icons[type],
    };
});

export function formatAdvisoryType(type: string): string {
    return type
        .split("_")
        .map((word) => word.charAt(0).toUpperCase() + word.slice(1))
        .join(" ");
}

export function formatAdvisoryCategory(category: string): string {
    return category.charAt(0).toUpperCase() + category.slice(1);
}

export function getAdvisoryTypeIcon(type: string): string {
    return advisories_type_icons[type] || "fa-circle-info";
}

export function getAdvisoryTypeClass(type: string): string {
    return advisories_type_classes[type] || "secondary";
}

export function getAdvisorySeverityIcon(severity: string): string {
    return advisories_type_icons[severity] || "fa-circle-info";
}

export function getAdvisorySeverityClass(severity: string): string {
    return advisories_type_classes[severity] || "secondary";
}

export function getAdvisoryCategoryIcon(category: string): string {
    switch (category) {
        case "pipelines":
            return "fa-project-diagram";
        case "modules":
            return "fa-cube";
        case "subworkflows":
            return "fa-sitemap";
        case "configuration":
            return "fa-cogs";
        default:
            return "fa-cogs";
    }
}

interface SoftwareDependency {
    name: string;
    versions?: string[];
}

type DependencyItem = string | SoftwareDependency;

export function formatNextflowVersions(versions: string[]): string {
    return versions.join(", ");
}

export function formatNextflowExecutors(executors: string[]): string {
    return executors.join(", ");
}

export function formatSoftwareDependency(dep: DependencyItem): string {
    if (typeof dep === "string") {
        return dep;
    }
    return `${dep.name}${dep.versions ? ` (${dep.versions.join(", ")})` : ""}`;
}

export function formatSoftwareDependencies(dependencies: DependencyItem[] | string): string {
    if (typeof dependencies === "string") {
        return dependencies;
    }
    return dependencies.map(formatSoftwareDependency).join(", ");
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
            icon: "fa-code-branch",
            label: "Nextflow",
            value: formatNextflowVersions(metadata.nextflowVersions),
        });
    }

    if (metadata.nextflowExecutors?.length) {
        items.push({
            icon: "fa-server",
            label: "Executors",
            value: formatNextflowExecutors(metadata.nextflowExecutors),
        });
    }

    if (metadata.softwareDependencies) {
        items.push({
            icon: "fa-object-group",
            label: "Dependencies",
            value: formatSoftwareDependencies(metadata.softwareDependencies),
        });
    }

    return items;
}
