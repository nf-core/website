export const advisories_type_classes = {
    // Advisory types
    known_regression: "success",
    incompatibility: "warning",
    security: "danger",
    performance: "success",
    data_corruption: "primary",
    scientific_advice: "secondary",
    other: "secondary",
    // Severity levels
    low: "success",
    medium: "warning",
    high: "danger",
    critical: "danger"
} as const;

export const advisories_type_icons = {
    // Advisory types
    known_regression: "fa-bug",
    incompatibility: "fa-exclamation-triangle",
    security: "fa-shield-alt",
    performance: "fa-tachometer-alt",
    data_corruption: "fa-database",
    scientific_advice: "fa-flask",
    other: "fa-info-circle",
    // Severity levels
    low: "fa-arrow-down",
    medium: "fa-minus",
    high: "fa-arrow-up",
    critical: "fa-exclamation-circle"
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
