"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
const request = require("request");
const responseToError = require("../responseToError");
const command_1 = require("../command");
const error_1 = require("../error");
const requirePermissions_1 = require("../requirePermissions");
const utils = require("../utils");
const api = require("../api");
const requireInstance = require("../requireInstance");
const settings_1 = require("../database/settings");
exports.default = new command_1.Command("database:settings:get <path>")
    .description("read the realtime database setting at path")
    .option("--instance <instance>", "use the database <instance>.firebaseio.com (if omitted, uses default database instance)")
    .help(settings_1.HELP_TEXT)
    .before(requirePermissions_1.requirePermissions, ["firebasedatabase.instances.get"])
    .before(requireInstance)
    .action((path, options) => {
    if (!settings_1.DATABASE_SETTINGS.has(path)) {
        return utils.reject(settings_1.INVALID_PATH_ERROR, { exit: 1 });
    }
    return new Promise((resolve, reject) => {
        const reqOptions = {
            url: utils.addSubdomain(api.realtimeOrigin, options.instance) + "/.settings/" + path + ".json",
        };
        return api.addRequestHeaders(reqOptions).then((reqOptionsWithToken) => {
            request.get(reqOptionsWithToken, (err, res, body) => {
                if (err) {
                    return reject(new error_1.FirebaseError(`Unexpected error fetching configs at ${path}`, {
                        exit: 2,
                        original: err,
                    }));
                }
                else if (res.statusCode >= 400) {
                    return reject(responseToError(res, body));
                }
                utils.logSuccess(`For database instance ${options.instance}\n\t ${path} = ${body}`);
                resolve();
            });
        });
    });
});
