"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
const _ = require("lodash");
function filterOnly(configs, onlyString) {
    if (!onlyString) {
        return configs;
    }
    let onlyTargets = onlyString.split(",");
    if (_.includes(onlyTargets, "hosting")) {
        return configs;
    }
    onlyTargets = onlyTargets
        .filter((target) => target.startsWith("hosting:"))
        .map((target) => target.replace("hosting:", ""));
    return configs.filter((config) => _.includes(onlyTargets, config.target || config.site));
}
function normalizedHostingConfigs(options) {
    let configs = options.config.get("hosting");
    if (!configs) {
        return [];
    }
    else if (!_.isArray(configs)) {
        if (!configs.target && !configs.site) {
            configs.site = options.instance;
        }
        configs = [configs];
    }
    return filterOnly(configs, options.only);
}
exports.normalizedHostingConfigs = normalizedHostingConfigs;
