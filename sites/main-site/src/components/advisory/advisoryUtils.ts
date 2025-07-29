// Combined dictionary for all advisory classes (types and severities)
export const advisoryClasses = {
    // Advisory types
    known_regression: "info",
    incompatibility: "warning",
    security: "danger",
    performance: "success",
    data_corruption: "primary",
    other: "secondary",
    // Severity levels
    low: "primary",
    medium: "warning",
    high: "danger",
    critical: "danger",
} as const;

// Combined dictionary for all advisory icons (types and severities)
export const advisoryIcons = {
    // Advisory types
    known_regression: "fa-solid fa-bug",
    incompatibility: "fa-solid fa-exclamation-triangle",
    security: "fa-solid fa-shield-alt",
    performance: "fa-solid fa-tachometer-alt",
    data_corruption: "fa-solid fa-database",
    other: "fa-solid fa-info-circle",
    // Severity levels
    low: "fa-solid fa-arrow-down",
    medium: "fa-solid fa-minus",
    high: "fa-solid fa-arrow-up",
    critical: "fa-solid fa-exclamation-circle",
    // Categories
    pipelines: "fa-solid fa-project-diagram",
    modules: "fa-solid fa-cube",
    subworkflows: "fa-solid fa-sitemap",
    configuration: "fa-solid fa-cogs",
} as const;

export const advisoryTypes = [
    "known_regression",
    "incompatibility",
    "security",
    "performance",
    "data_corruption",
    "other",
].map((type) => {
    return {
        name: type,
        class: advisoryClasses[type],
        icon: advisoryIcons[type],
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
            icon: "logos:nextflow",
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
