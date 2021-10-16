"use strict";
var __awaiter = (this && this.__awaiter) || function (thisArg, _arguments, P, generator) {
    function adopt(value) { return value instanceof P ? value : new P(function (resolve) { resolve(value); }); }
    return new (P || (P = Promise))(function (resolve, reject) {
        function fulfilled(value) { try { step(generator.next(value)); } catch (e) { reject(e); } }
        function rejected(value) { try { step(generator["throw"](value)); } catch (e) { reject(e); } }
        function step(result) { result.done ? resolve(result.value) : adopt(result.value).then(fulfilled, rejected); }
        step((generator = generator.apply(thisArg, _arguments || [])).next());
    });
};
Object.defineProperty(exports, "__esModule", { value: true });
const _ = require("lodash");
const opn = require("open");
const qs = require("querystring");
const command_1 = require("../command");
const error_1 = require("../error");
const gcp = require("../gcp");
const getProjectId = require("../getProjectId");
const logger = require("../logger");
const requirePermissions_1 = require("../requirePermissions");
module.exports = new command_1.Command("functions:log")
    .description("read logs from deployed functions")
    .option("--only <function_names>", 'only show logs of specified, comma-seperated functions (e.g. "funcA,funcB")')
    .option("-n, --lines <num_lines>", "specify number of log lines to fetch")
    .option("--open", "open logs page in web browser")
    .before(requirePermissions_1.requirePermissions, ["logging.logEntries.list", "logging.logs.list"])
    .action((options) => __awaiter(void 0, void 0, void 0, function* () {
    try {
        const projectId = getProjectId(options, false);
        let apiFilter = `resource.type="cloud_function"`;
        if (options.only) {
            const funcNames = options.only.split(",");
            const apiFuncFilters = _.map(funcNames, (funcName) => {
                return `resource.labels.function_name="${funcName}"`;
            });
            apiFilter += `\n(${apiFuncFilters.join(" OR ")})`;
        }
        if (options.open) {
            const url = `https://console.developers.google.com/logs/viewer?advancedFilter=${qs.escape(apiFilter)}&project=${projectId}`;
            opn(url);
            return;
        }
        const entries = yield gcp.cloudlogging.listEntries(projectId, apiFilter, options.lines || 35, "desc");
        for (let i = _.size(entries) - 1; i >= 0; i--) {
            const entry = entries[i];
            logger.info(entry.timestamp, _.get(entry, "severity", "?").substring(0, 1), _.get(entry, "resource.labels.function_name") + ":", _.get(entry, "textPayload", ""));
        }
        if (_.isEmpty(entries)) {
            logger.info("No log entries found.");
        }
        return entries;
    }
    catch (err) {
        throw new error_1.FirebaseError(`Failed to list log entries ${err.message}`, { exit: 1 });
    }
}));
