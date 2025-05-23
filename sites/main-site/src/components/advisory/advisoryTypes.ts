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
    critical: "danger"
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
    critical: "fa-solid fa-exclamation-circle"
} as const;

// Only include actual advisory types, not severity levels
const advisory_types = [
    'known_regression',
    'incompatibility',
    'security',
    'performance',
    'data_corruption',
    'scientific_advice',
    'other'
] as const;

export const advisories_types = advisory_types.map((type) => {
    return {
        name: type,
        class: advisories_type_classes[type],
        icon: advisories_type_icons[type],
    };
});
